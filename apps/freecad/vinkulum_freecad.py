"""FreeCAD owns modelling and navigation; Vinkulum supplies captured mechanics.

SPDX-License-Identifier: Apache-2.0
"""

import bisect
import importlib.util
import json
import os
import time
import uuid
from pathlib import Path

import FreeCAD as App
import FreeCADGui as Gui
from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

ROOT = Path(__file__).resolve().parent
ICON = str(ROOT / "vinkulum.svg")
_spec = importlib.util.spec_from_file_location(
    "vinkulum_freecad_bridge", ROOT / "bridge.py"
)
bridge = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bridge)
_panel = None


def preferences():
    return App.ParamGet("User parameter:BaseApp/Preferences/Mod/Vinkulum")


def worker_python():
    return preferences().GetString(
        "WorkerPython", os.environ.get("VINKULUM_FREECAD_PYTHON", "")
    )


def calculation_root():
    default = str(Path(App.getUserAppDataDir()) / "Vinkulum" / "calculations")
    return Path(preferences().GetString("CalculationDirectory", default))


def selected_source():
    objects = Gui.Selection.getSelection()
    if len(objects) != 1:
        raise ValueError(
            "Sélectionnez une pièce ou un résultat Vinkulum dans l’arbre FreeCAD."
        )
    source = objects[0]
    if hasattr(source, "VinkulumRequest"):
        return source
    parent = source.getParentGeoFeatureGroup()
    if parent is not None and parent.isDerivedFrom("PartDesign::Body"):
        source = parent
    if not hasattr(source, "Shape") or len(source.Shape.Solids) != 1:
        raise ValueError("Sélectionnez une pièce comportant un seul solide.")
    return source


def number(value, low, high, decimals=6, suffix=""):
    widget = QDoubleSpinBox()
    widget.setRange(low, high)
    widget.setDecimals(decimals)
    widget.setValue(value)
    widget.setSuffix(suffix)
    widget.setKeyboardTracking(False)
    return widget


class AnalysisPanel:
    def __init__(self, selected):
        self.form = QWidget()
        self.form.setWindowTitle("Vinkulum · pièce rigide")
        self.form.setObjectName("vinkulum_analysis")
        self.source = (
            selected.VinkulumSource if hasattr(selected, "VinkulumSource") else selected
        )
        if self.source is None:
            raise ValueError("La pièce source n’existe plus.")
        self.document_name, self.source_name = (
            self.source.Document.Name,
            self.source.Name,
        )
        self._source_visibility = getattr(
            selected, "VinkulumSourceWasVisible", self.source.Visibility
        )
        self.job = self.record = self.result_object = None
        self._closing = self._closed = False
        self._origin = None
        layout = QVBoxLayout(self.form)
        layout.addWidget(QLabel(self.source.Label))
        layout.addWidget(QLabel("Pièce rigide · pivot fixé au monde"))
        inputs = QFormLayout()
        self.density = number(7800, 0.001, 1e8, 3, " kg/m³")
        inputs.addRow("Masse volumique", self.density)
        self.pivot = [number(0, -1e9, 1e9, suffix=" mm") for _ in range(3)]
        for axis, field in zip("XYZ", self.pivot):
            inputs.addRow(f"Pivot {axis} · repère global", field)
        self.axis = QComboBox()
        self.axis.addItems(["X", "Y", "Z"])
        self.axis.setCurrentIndex(1)
        inputs.addRow("Axe du pivot · repère global", self.axis)
        self.duration = number(2, 0.001, 10, suffix=" s")
        self.step = number(0.005, 0.000001, 10, suffix=" s")
        inputs.addRow("Durée", self.duration)
        inputs.addRow("Pas de calcul", self.step)
        layout.addLayout(inputs)
        self.run_button = QPushButton("Calculer")
        self.run_button.setObjectName("vinkulum_run")
        self.cancel_button = QPushButton("Annuler le calcul")
        self.cancel_button.setObjectName("vinkulum_cancel")
        self.play_button = QPushButton("Lire le mouvement")
        self.play_button.setObjectName("vinkulum_play")
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setObjectName("vinkulum_time")
        self.time_label = QLabel("t = 0 s")
        self.status = QLabel("Sélection capturée au lancement du calcul.")
        self.status.setWordWrap(True)
        for item in (
            self.run_button,
            self.cancel_button,
            self.play_button,
            self.slider,
            self.time_label,
            self.status,
        ):
            layout.addWidget(item)
        self.timer = QTimer(self.form)
        self.timer.setInterval(16)
        self.timer.timeout.connect(self.tick)
        self.run_button.clicked.connect(self.run)
        self.cancel_button.clicked.connect(self.cancel)
        self.play_button.clicked.connect(self.play)
        self.slider.valueChanged.connect(self.frame)
        App.addDocumentObserver(self)
        if hasattr(selected, "VinkulumRequest"):
            try:
                self.reopen(selected)
            except Exception:
                App.removeDocumentObserver(self)
                raise
        self.controls()

    def controls(self):
        busy = self.job is not None and self.job.running
        self.run_button.setEnabled(not busy and not self._closing)
        self.cancel_button.setEnabled(busy)
        self.play_button.setEnabled(not busy and self.record is not None)
        self.slider.setEnabled(not busy and self.record is not None)
        for field in (self.density, self.axis, self.duration, self.step, *self.pivot):
            field.setEnabled(not busy)

    def run(self):
        if self._closing or (self.job is not None and self.job.running):
            return
        self.stop()
        try:
            if self.job is not None:
                self.job.deleteLater()
                self.job = None
            interpreter = worker_python()
            if not interpreter or not Path(interpreter).is_file():
                raise ValueError(
                    "Configurez le moteur dans Vinkulum → Configurer le moteur…"
                )
            directory = calculation_root() / str(uuid.uuid4())
            axis = [float(i == self.axis.currentIndex()) for i in range(3)]
            request = bridge.capture(
                self.source,
                directory,
                density=self.density.value(),
                body_id=str(uuid.uuid4()),
                project_id=str(uuid.uuid4()),
                duration=self.duration.value(),
                step=self.step.value(),
                pivot_m=[p.value() * 0.001 for p in self.pivot],
                axis_world=axis,
            )
            self.job = bridge.Job(self.source, directory, request, interpreter)
            self.job.completed.connect(self.completed)
            self.job.failed.connect(self.status.setText)
            self.job.cancelled.connect(self.status.setText)
            self.job.ended.connect(self.ended)
            self.job.start()
            self.status.setText("Calcul en cours…")
        except Exception as error:
            self.status.setText(str(error))
        self.controls()

    def cancel(self):
        if self.job is not None:
            self.job.cancel("Calcul annulé ; le résultat précédent est conservé.")

    def ended(self):
        self.controls()
        if self._closing:
            QTimer.singleShot(0, self.reject)

    def completed(self, result):
        if self._closing:
            return
        doc = self.source.Document
        doc.openTransaction("Résultat Vinkulum")
        try:
            obj = doc.addObject("Part::Feature", "VinkulumResult")
            obj.Label = self.source.Label + " · résultat Vinkulum"
            obj.Shape = self.source.Shape.copy()
            captured_brep = obj.Shape.exportBrepToString()
            obj.setEditorMode("Shape", 1)
            obj.ViewObject.ShapeColor = (0.95, 0.57, 0.15)
            for kind, name, value in (
                ("App::PropertyLink", "VinkulumSource", self.source),
                (
                    "App::PropertyBool",
                    "VinkulumSourceWasVisible",
                    self._source_visibility,
                ),
                ("App::PropertyString", "VinkulumDirectory", str(self.job.directory)),
                (
                    "App::PropertyString",
                    "VinkulumRequest",
                    bridge.request_bytes(self.job.request).decode("ascii"),
                ),
                (
                    "App::PropertyString",
                    "VinkulumPlaybackHash",
                    bridge.sha(
                        (self.job.directory / "result/playback.json").read_bytes()
                    ),
                ),
                (
                    "App::PropertyPlacement",
                    "VinkulumInitialPlacement",
                    App.Placement(obj.Placement),
                ),
                ("App::PropertyString", "VinkulumCapturedBrep", captured_brep),
                (
                    "App::PropertyString",
                    "VinkulumShapeHash",
                    bridge.sha(captured_brep.encode("ascii")),
                ),
            ):
                obj.addProperty(kind, name, "Vinkulum")
                setattr(obj, name, value)
                obj.setEditorMode(
                    name, 1 if name in ("VinkulumSource", "VinkulumDirectory") else 2
                )
            doc.recompute()
            doc.commitTransaction()
        except Exception as error:
            doc.abortTransaction()
            self.status.setText(str(error))
            return
        if self.result_object is not None:
            self.result_object.Visibility = False
        self.result_object = obj
        self.record = (self.job.request, result)
        self._set_result()
        self.status.setText(
            f"Calcul terminé · {len(result['time_s'])} échantillons · résultat conservé."
        )

    def reopen(self, obj):
        import Part

        directory = Path(obj.VinkulumDirectory)
        if (
            bridge.sha((directory / "result/playback.json").read_bytes())
            != obj.VinkulumPlaybackHash
        ):
            raise ValueError("Le fichier de mouvement a changé depuis le calcul.")
        brep = obj.VinkulumCapturedBrep
        if bridge.sha(brep.encode("ascii")) != obj.VinkulumShapeHash:
            raise ValueError("La géométrie capturée du résultat a changé.")
        request = json.loads(obj.VinkulumRequest)
        result = bridge.read_playback(directory, request)
        captured_shape = Part.Shape()
        captured_shape.importBrepFromString(brep)
        if not captured_shape.isValid() or len(captured_shape.Solids) != 1:
            raise ValueError("La géométrie capturée n’est plus un solide valide.")
        if any(
            abs(a * 0.001 - b) > 1e-9
            for a, b in zip(
                captured_shape.Solids[0].CenterOfMass,
                request["properties_si"]["centre_m"],
            )
        ):
            raise ValueError(
                "Le centre de la géométrie capturée ne correspond pas au résultat."
            )
        # FreeCAD may normalize its shape serialization when saving FCStd.
        # Restore the recorded result geometry, independently of that cache.
        obj.Shape = captured_shape
        obj.VinkulumInitialPlacement = App.Placement(obj.Placement)
        self.record, self.result_object = (request, result), obj
        self.density.setValue(request["density_kg_m3"])
        self.duration.setValue(request["duration_s"])
        self.step.setValue(request["step_s"])
        for widget, value in zip(self.pivot, request["pivot_m"]):
            widget.setValue(value * 1000)
        self.axis.setCurrentIndex(
            max(range(3), key=lambda i: abs(request["axis_world"][i]))
        )
        self._set_result()
        self.status.setText("Résultat conservé · géométrie capturée lors du calcul.")

    def _set_result(self):
        self.slider.setRange(0, len(self.record[1]["time_s"]) - 1)
        self.slider.setValue(0)
        self.frame(0)
        self.controls()

    def frame(self, index):
        if self.record is None or self._closing:
            return
        request, result = self.record
        try:
            self.result_object.Visibility = True
            self.source.Visibility = False
            bridge.set_pose(
                self.result_object,
                self.result_object.VinkulumInitialPlacement,
                request["properties_si"]["centre_m"],
                result["position_m"][index],
                result["rotation"][index],
            )
            self.time_label.setText(f"t = {result['time_s'][index]:.6g} s")
        except (RuntimeError, ReferenceError) as error:
            self.stop()
            self.status.setText(str(error))

    def play(self):
        if self.timer.isActive():
            self.stop()
            return
        if self.record is None:
            return
        if self.slider.value() == self.slider.maximum():
            self.slider.setValue(0)
        self._origin = time.monotonic() - self.record[1]["time_s"][self.slider.value()]
        self.timer.start()
        self.play_button.setText("Pause")

    def tick(self):
        times = self.record[1]["time_s"]
        elapsed = time.monotonic() - self._origin
        index = min(len(times) - 1, bisect.bisect_right(times, elapsed) - 1)
        self.slider.setValue(max(0, index))
        if index == len(times) - 1:
            self.stop()

    def stop(self):
        self.timer.stop()
        self.play_button.setText("Lire le mouvement")

    def slotDeletedDocument(self, document):
        if document.Name == self.document_name:
            self._closing = True
            self.stop()
            self.cancel()
            QTimer.singleShot(0, self.reject)

    def slotDeletedObject(self, obj):
        if obj.Document.Name == self.document_name and (
            obj.Name == self.source_name or obj is self.result_object
        ):
            self._closing = True
            self.stop()
            self.cancel()
            QTimer.singleShot(0, self.reject)

    def getStandardButtons(self):
        return QDialogButtonBox.StandardButton.Close.value

    def isAllowedAlterDocument(self):
        return True

    def isAllowedAlterSelection(self):
        return True

    def isAllowedAlterView(self):
        return True

    def reject(self):
        global _panel
        if self._closed:
            return True
        self._closing = True
        self.stop()
        if self.job is not None and self.job.running:
            self.cancel()
            return False
        self._closed = True
        App.removeDocumentObserver(self)
        doc = App.listDocuments().get(self.document_name)
        if doc is not None:
            source = doc.getObject(self.source_name)
            if source is not None:
                source.Visibility = self._source_visibility
            if self.result_object is not None:
                try:
                    self.result_object.Placement = (
                        self.result_object.VinkulumInitialPlacement
                    )
                    self.result_object.Visibility = False
                except (ReferenceError, RuntimeError):
                    pass
        if self.job is not None:
            self.job.deleteLater()
        _panel = None
        Gui.Control.closeDialog()
        return True


def open_analysis():
    global _panel
    if _panel is not None:
        return
    if Gui.Control.activeDialog():
        raise ValueError(
            "Terminez la tâche FreeCAD ouverte avant de créer une analyse."
        )
    _panel = AnalysisPanel(selected_source())
    Gui.Control.showDialog(_panel)


def configure():
    dialog = QDialog(Gui.getMainWindow())
    dialog.setWindowTitle("Moteur Vinkulum")
    layout = QFormLayout(dialog)
    fields = []
    for label, value, directory in (
        ("Python du moteur", worker_python(), False),
        ("Dossier des calculs", str(calculation_root()), True),
    ):
        line = QLineEdit(value)
        row = QHBoxLayout()
        row.addWidget(line)
        choose = QPushButton("Parcourir…")

        def pick(checked=False, field=line, folder=directory):
            value = (
                QFileDialog.getExistingDirectory(
                    dialog, "Dossier des calculs", field.text()
                )
                if folder
                else QFileDialog.getOpenFileName(
                    dialog, "Python du moteur", field.text()
                )[0]
            )
            if value:
                field.setText(value)

        choose.clicked.connect(pick)
        row.addWidget(choose)
        layout.addRow(label, row)
        fields.append(line)
    buttons = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
    )
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addRow(buttons)
    if dialog.exec():
        for key, field in zip(("WorkerPython", "CalculationDirectory"), fields):
            preferences().SetString(key, field.text())


def example():
    import Part
    import Sketcher

    doc = App.newDocument("PenduleVinkulum")
    body = doc.addObject("PartDesign::Body", "Pendule")
    sketch = body.newObject("Sketcher::SketchObject", "Section")
    points = [(-10, -15), (10, -15), (10, 15), (-10, 15)]
    for a, b in zip(points, points[1:] + points[:1]):
        i = sketch.addGeometry(
            Part.LineSegment(App.Vector(*a, 0), App.Vector(*b, 0)), False
        )
        sketch.addConstraint(Sketcher.Constraint("Block", i))
    pad = body.newObject("PartDesign::Pad", "Longueur")
    pad.Profile, pad.Length = sketch, 800
    body.Placement = App.Placement(App.Vector(), App.Rotation(App.Vector(0, 1, 0), 160))
    doc.recompute()
    sketch.Visibility = False
    view = Gui.activeDocument().activeView()
    animated = view.isAnimationEnabled()
    view.setAnimationEnabled(False)
    try:
        view.viewFront()
    finally:
        view.setAnimationEnabled(animated)
    Gui.Selection.clearSelection()
    Gui.Selection.addSelection(body)

    def fit_example():
        if App.ActiveDocument is doc:
            # An explicit margin uses FreeCAD's immediate bounding-box fit;
            # default animated fitAll can retain a viewer across document close.
            Gui.activeDocument().activeView().fitAll(1.15)

    QTimer.singleShot(0, fit_example)


class Command:
    def __init__(self, title, callback, needs_document=False):
        self.title, self.callback, self.needs_document = title, callback, needs_document

    def GetResources(self):
        return {"MenuText": self.title, "ToolTip": self.title, "Pixmap": ICON}

    def IsActive(self):
        return not self.needs_document or App.ActiveDocument is not None

    def Activated(self):
        try:
            self.callback()
        except Exception as error:
            QMessageBox.warning(Gui.getMainWindow(), "Vinkulum", str(error))


def register_commands():
    commands = {
        "Vinkulum_Analysis": Command(
            "Analyser la pièce sélectionnée…", open_analysis, True
        ),
        "Vinkulum_Example": Command("Exemple : pendule paramétrique", example),
        "Vinkulum_Configure": Command("Configurer le moteur…", configure),
    }
    for name, command in commands.items():
        Gui.addCommand(name, command)
    return list(commands)
