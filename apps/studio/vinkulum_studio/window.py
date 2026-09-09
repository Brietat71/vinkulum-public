"""G0 desktop workflow, deliberately independent of the native solver state."""

from bisect import bisect_right
from dataclasses import asdict
import json
import math
import time
from importlib.metadata import version

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QHBoxLayout, QLabel, QMainWindow,
                               QPushButton, QSlider, QTextEdit, QVBoxLayout, QWidget)

from . import __version__
from .controller import Controller
from .model import Parameters, load_parameters, save_parameters
from .views import AngleView, PendulumView


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"Vinkulum Studio {__version__} · Pendule G0")
        self.resize(1100, 820)
        self.controller = Controller(self)
        self.result = None
        self._previous = False
        self._times = ()
        self.timer = QTimer(self)
        self.timer.setInterval(33)
        self.timer.timeout.connect(self._tick)
        root = QWidget()
        self.setCentralWidget(root)
        layout = QVBoxLayout(root)
        layout.setContentsMargins(24, 20, 24, 18)
        title = QLabel("VINKULUM  /  STUDIO")
        title.setStyleSheet("font-size: 23px; font-weight: 700; color: #17314a;")
        layout.addWidget(title)
        layout.addWidget(QLabel(f"Prototype G0 · Un pendule calculé par le noyau Vinkulum {version('vinkulum')}"))
        columns = QHBoxLayout()
        layout.addLayout(columns, 1)
        panel = QWidget()
        panel.setFixedWidth(290)
        left = QVBoxLayout(panel)
        left.setContentsMargins(0, 15, 18, 0)
        left.addWidget(QLabel("PARAMÈTRES DU PROCHAIN CALCUL"))
        form = QFormLayout()
        self.fields = {}
        specs = [("length", "Longueur [m]", .05, 20., 4, .1),
                 ("mass", "Masse [kg]", .01, 100., 4, .1),
                 ("angle_deg", "Angle initial [°]", -170., 170., 2, 5.),
                 ("duration", "Durée [s]", .001, 600., 3, .5),
                 ("step", "Pas [s]", .000001, 1., 6, .001)]
        defaults = Parameters()
        for name, label, low, high, decimals, increment in specs:
            field = QDoubleSpinBox()
            field.setRange(low, high)
            field.setDecimals(decimals)
            field.setSingleStep(increment)
            field.setKeyboardTracking(False)
            field.setValue(getattr(defaults, name))
            field.setAccessibleName(label)
            field.valueChanged.connect(self._edited)
            self.fields[name] = field
            form.addRow(label, field)
        left.addLayout(form)
        self.run_button = QPushButton("Lancer le calcul")
        self.run_button.setShortcut("Ctrl+Return")
        self.run_button.setStyleSheet("background: #187b8d; color: white; padding: 9px; font-weight: 600;")
        self.run_button.clicked.connect(self.run)
        self.stop_button = QPushButton("Arrêter")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop)
        left.addWidget(self.run_button)
        left.addWidget(self.stop_button)
        for label, slot in [("Enregistrer les paramètres…", self._save_dialog),
                            ("Ouvrir les paramètres…", self._load_dialog)]:
            button = QPushButton(label)
            button.clicked.connect(slot)
            left.addWidget(button)
        self.export_button = QPushButton("Exporter les échantillons CSV…")
        self.export_button.clicked.connect(self._export_dialog)
        self.export_button.setEnabled(False)
        left.addWidget(self.export_button)
        self.provenance_button = QPushButton("Consulter la provenance")
        self.provenance_button.clicked.connect(self._provenance)
        self.provenance_button.setEnabled(False)
        left.addWidget(self.provenance_button)
        note = QLabel("Les modifications préparent le prochain calcul.\n\n"
                      "Gravité : 9,80665 m/s²\nVitesse initiale : 0 rad/s\n"
                      "α-généralisé · ρ∞ = 0,9\nInertie : 10⁻⁸ kg·m²\n\n"
                      "Limite G0 : 20 000 pas.\nLe dessin utilise les échantillons calculés.")
        note.setWordWrap(True)
        left.addWidget(note)
        left.addStretch()
        columns.addWidget(panel)
        right = QVBoxLayout()
        columns.addLayout(right, 1)
        self.result_label = QLabel()
        self.result_label.setWordWrap(True)
        self.result_label.setMinimumHeight(58)
        right.addWidget(self.result_label)
        self.pendulum = PendulumView()
        self.curve = AngleView()
        right.addWidget(self.pendulum, 1)
        right.addWidget(self.curve, 1)
        playback = QHBoxLayout()
        self.play_button = QPushButton("Lecture")
        self.play_button.clicked.connect(self.toggle_play)
        self.play_button.setEnabled(False)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setAccessibleName("Échantillon temporel du résultat")
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self._frame)
        self.slider.sliderPressed.connect(self.pause)
        self.time_label = QLabel("t = —")
        playback.addWidget(self.play_button)
        playback.addWidget(self.slider, 1)
        playback.addWidget(self.time_label)
        right.addLayout(playback)
        self.status = QLabel("Prêt. L'exemple est chargé ; lancez le premier calcul.")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        self.status.setAccessibleName("État du calcul")
        self.status.setMinimumHeight(42)
        layout.addWidget(self.status)
        limitations = QLabel("Domaine : démonstration de pendule plan. Calcul terminé ≠ précision certifiée. "
                             "Erreur de trajectoire non évaluée ; aucune reprise de calcul interrompu.")
        limitations.setWordWrap(True)
        limitations.setStyleSheet("color: #5c6573; font-size: 11px;")
        layout.addWidget(limitations)
        self.controller.busy_changed.connect(self._busy)
        self.controller.completed.connect(self._completed)
        self.controller.problem.connect(self._problem)
        self._edited()

    def parameters(self):
        for field in self.fields.values():
            field.interpretText()
        return Parameters(**{key: field.value() for key, field in self.fields.items()})

    def set_parameters(self, parameters):
        for key, value in asdict(parameters).items():
            self.fields[key].setValue(value)

    def _edited(self):
        if self.result is None:
            length = self.fields["length"].value()
            angle = math.radians(self.fields["angle_deg"].value())
            self.pendulum.length = length
            self.pendulum.sample = (0., length * math.sin(angle), -length * math.cos(angle), angle)
            self.pendulum.update()
        self._result_caption()

    def _result_caption(self):
        if self.result is None:
            self.result_label.setText("Aperçu de l'état initial · Aucun résultat calculé.")
            return
        p = self.result.parameters
        draft_differs = any(self.fields[k].value() != v for k, v in asdict(p).items())
        state = "Résultat précédent conservé" if self._previous else "Résultat calculé"
        if draft_differs:
            state += " · paramètres différents du brouillon"
        self.result_label.setText(f"{state} · {self.result.run_id[:8]}\n"
                                  f"L = {p.length:g} m · m = {p.mass:g} kg · θ₀ = {p.angle_deg:g}° · "
                                  f"durée = {p.duration:g} s · pas = {p.step:g} s")

    def run(self):
        try:
            parameters = self.parameters()
            self.controller.start(parameters)
            self._previous = self.result is not None
            self._result_caption()
            self.status.setText("Calcul en cours dans un processus séparé. Les champs restent éditables.")
        except (ValueError, RuntimeError, OSError) as exc:
            self._problem(str(exc))

    def stop(self):
        self.controller.cancel()
        self.stop_button.setEnabled(False)
        self.status.setText("Arrêt demandé…")

    def _busy(self, busy):
        self.run_button.setEnabled(not busy)
        self.stop_button.setEnabled(busy)

    def _completed(self, result):
        self.pause()
        self.result = result
        self._previous = False
        self._times = tuple(row[0] for row in result.samples)
        self.pendulum.length = result.parameters.length
        self.curve.samples = result.samples
        self.slider.setRange(0, len(result.samples) - 1)
        self.slider.setValue(0)
        self._frame(0)
        for control in (self.play_button, self.slider, self.export_button, self.provenance_button):
            control.setEnabled(True)
        self._result_caption()
        self.status.setText(f"Calcul terminé · {len(result.samples)} échantillons. "
                            "Contrôles de cohérence réussis ; erreur numérique non évaluée.")

    def _problem(self, message):
        self._previous = self.result is not None
        self._result_caption()
        self.status.setText(message + (" Le résultat précédent reste affiché." if self.result else ""))

    def _frame(self, index):
        if self.result is None:
            return
        self.pendulum.sample = self.result.samples[index]
        self.curve.index = index
        self.pendulum.update()
        self.curve.update()
        self.time_label.setText(f"t = {self.result.samples[index][0]:.4f} s")
        if self.timer.isActive():
            self._play_origin = time.monotonic() - self.result.samples[index][0]

    def toggle_play(self):
        if self.timer.isActive():
            self.pause()
        elif self.result:
            if self.slider.value() == self.slider.maximum():
                self.slider.setValue(0)
            self._play_origin = time.monotonic() - self.result.samples[self.slider.value()][0]
            self.timer.start()
            self.play_button.setText("Pause")

    def pause(self):
        self.timer.stop()
        self.play_button.setText("Lecture")

    def _tick(self):
        elapsed = time.monotonic() - self._play_origin
        index = min(len(self._times) - 1, max(0, bisect_right(self._times, elapsed) - 1))
        origin = self._play_origin
        self.slider.setValue(index)
        self._play_origin = origin
        if elapsed >= self._times[-1]:
            self.pause()

    def save(self, path):
        save_parameters(path, self.parameters())
        self.status.setText("Paramètres du prochain calcul enregistrés.")

    def load(self, path):
        parameters = load_parameters(path)  # Fully validate before editing any widget.
        for key, value in asdict(parameters).items():
            if round(value, self.fields[key].decimals()) != value:
                raise ValueError(f"{key} : précision supérieure aux {self.fields[key].decimals()} "
                                 "décimales saisissables dans G0 ; chargement refusé sans arrondi.")
        self.set_parameters(parameters)
        self.status.setText("Paramètres chargés pour le prochain calcul.")

    def _file_action(self, save, caption, file_filter, action):
        picker = QFileDialog.getSaveFileName if save else QFileDialog.getOpenFileName
        path, _ = picker(self, caption, "", file_filter)
        if path:
            try:
                action(path)
            except (OSError, ValueError) as exc:
                self.status.setText(f"Opération impossible : {exc}")

    def _save_dialog(self):
        self._file_action(True, "Enregistrer les paramètres", "Paramètres (*.json)", self.save)

    def _load_dialog(self):
        self._file_action(False, "Ouvrir les paramètres", "Paramètres (*.json)", self.load)

    def _export_dialog(self):
        def export(path):
            self.result.export_csv(path)
            self.status.setText("Échantillons et provenance du résultat affiché exportés.")
        self._file_action(True, "Exporter le résultat affiché", "Échantillons (*.csv)", export)

    def _provenance(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Provenance du résultat affiché")
        dialog.resize(650, 520)
        layout = QVBoxLayout(dialog)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText(json.dumps({"run_id": self.result.run_id,
                                      "parameters": asdict(self.result.parameters),
                                      "manifest": json.loads(self.result.manifest_json)},
                                     ensure_ascii=False, indent=2))
        layout.addWidget(text)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec()

    def closeEvent(self, event):
        self.pause()
        self.controller.shutdown()
        event.accept()
