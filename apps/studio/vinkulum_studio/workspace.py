"""Dock workspace and discoverable commands for the mechanical editor."""

from PySide6.QtCore import QSize, Qt, QTimer
from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDockWidget,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QPushButton,
    QScrollArea,
    QSlider,
    QTabBar,
    QTableView,
    QTabWidget,
    QToolButton,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)

from .commands import CommandPalette, search_key
from .controls import PATHS, line_icon
from .document import Body
from .examples3d import EXAMPLES
from .run_archive import model_differences, series_keys
from .series_view import SampleTableModel, SeriesView
from .theme import apply_theme
from .viewport import Viewport


class Workspace:
    def _action(self, key, label, slot, shortcut=None, *, design=False, tip=None):
        action = QAction(label, self)
        action.setObjectName(key)
        if key in PATHS:
            action.setIcon(line_icon(key))
        action.setToolTip(tip or label)
        action.setStatusTip(tip or label)
        if shortcut:
            action.setShortcut(QKeySequence(shortcut))
        action.triggered.connect(slot)
        self.addAction(action)
        self.commands[key] = action
        if design:
            self._design_actions.append(action)
        return action

    def _dock(self, key, title, widget, area):
        dock = QDockWidget(title, self)
        dock.setObjectName(key)
        dock.setWidget(widget)
        self.addDockWidget(area, dock)
        self.docks[key] = dock
        return dock

    def _make_layout(self):
        self.commands, self.docks, self._design_actions = {}, {}, []
        self.setDockNestingEnabled(True)
        font = self.font()
        font.setPointSizeF(max(10.0, font.pointSizeF()))
        self.setFont(font)
        self.setMinimumSize(1000, 700)
        files = self.menuBar().addMenu("&File")
        edit = self.menuBar().addMenu("&Edit")
        create = self.menuBar().addMenu("&Create")
        view = self.menuBar().addMenu("&View")
        analysis = self.menuBar().addMenu("&Run")
        help_menu = self.menuBar().addMenu("&Help")
        for key, label, slot, shortcut in (
            ("new", "New project", self.new_project, QKeySequence.StandardKey.New),
            (
                "open",
                "Open project…",
                self.open_dialog,
                QKeySequence.StandardKey.Open,
            ),
            (
                "save",
                "Save project…",
                self.save_dialog,
                QKeySequence.StandardKey.Save,
            ),
            ("import", "Import G0 pendulum…", self.import_dialog, None),
        ):
            files.addAction(self._action(key, label, slot, shortcut))
        examples = files.addMenu("Open example")
        for i, name in enumerate(EXAMPLES):
            examples.addAction(
                self._action(
                    f"example_{i}",
                    name,
                    lambda checked=False, n=name: self.load_example(n),
                )
            )
        for key, label, slot, shortcut in (
            ("undo", "Undo", self.undo, QKeySequence.StandardKey.Undo),
            ("redo", "Redo", self.redo, QKeySequence.StandardKey.Redo),
            ("duplicate", "Duplicate selection", self.duplicate, "Ctrl+D"),
            ("delete", "Delete selection", self.delete, "Ctrl+Delete"),
        ):
            edit.addAction(self._action(key, label, slot, shortcut, design=True))
        for shape, label, shortcut in (
            ("box", "Box", "Alt+B"),
            ("cylinder", "Cylinder", "Alt+C"),
            ("sphere", "Sphere", "Alt+S"),
        ):
            create.addAction(
                self._action(
                    shape,
                    f"Add {label.lower()}" if shape != "cylinder" else "Add cylinder",
                    lambda checked=False, s=shape: self.add_body(s),
                    shortcut,
                    design=True,
                )
            )
        create.addSeparator()
        create.addAction(
            self._action("cad", "CAD design…", self.open_cad, "Alt+G", design=True)
        )
        create.addSeparator()
        create.addAction(
            self._action(
                "joint",
                "Add joint…",
                self.add_joint_dialog,
                "Alt+L",
                design=True,
            )
        )
        create.addAction(
            self._action("load", "Add load", self.add_load, "Alt+F", design=True)
        )
        toolbar = self.addToolBar("Main tools")
        toolbar.setObjectName("main_tools")
        toolbar.setMovable(False)
        brand = QLabel("Vinkulum  Studio")
        brand.setObjectName("brand")
        toolbar.addWidget(brand)
        toolbar.addSeparator()
        self.workspace_tabs = QTabBar()
        self.workspace_tabs.setAccessibleName("Workspace")
        for label in ("Model", "Simulate", "Inspect"):
            self.workspace_tabs.addTab(label)
        self.workspace_tabs.setTabEnabled(2, False)
        toolbar.addWidget(self.workspace_tabs)
        toolbar.addSeparator()
        toolbar.setIconSize(QSize(20, 20))
        toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        for key in ("open", "save", "undo", "redo"):
            toolbar.addAction(self.commands[key])
        toolbar.addSeparator()
        add = QToolButton()
        add.setText("+ Create")
        add.setMenu(create)
        add.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(add)
        cad_button = QToolButton()
        cad_button.setToolTip("Solid modelling · Alt+G")
        cad_button.setDefaultAction(self.commands["cad"])
        cad_button.setText("CAD")
        cad_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        toolbar.addWidget(cad_button)
        example_button = QToolButton()
        example_button.setText("Examples")
        example_button.setMenu(examples)
        example_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(example_button)
        stretch = QWidget()
        from PySide6.QtWidgets import QSizePolicy

        stretch.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        toolbar.addWidget(stretch)
        palette_action = self._action(
            "commands", "Find a command…", self.show_commands, "Ctrl+K"
        )
        palette_action.setIcon(line_icon("search"))
        toolbar.addAction(palette_action)
        help_menu.addAction(palette_action)
        self._action("context_tools", "Selection tools…", self.show_context_tools, "S")

        objects = QWidget()
        objects_layout = QVBoxLayout(objects)
        objects_layout.setContentsMargins(8, 8, 8, 6)
        objects_layout.setSpacing(6)
        self.project_label = QLabel()
        self.project_label.setWordWrap(True)
        objects_layout.addWidget(self.project_label)
        self.tree_search = QLineEdit()
        self.tree_search.setPlaceholderText("Filter objects…")
        self.tree_search.setAccessibleName("Filter objects by name, type or identity")
        self.tree_search.setClearButtonEnabled(True)
        self.tree_search.textChanged.connect(self._filter_tree)
        objects_layout.addWidget(self.tree_search)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["MECHANISM", "Type"])
        self.tree.setAccessibleName("Mechanism browser")
        self.tree.setUniformRowHeights(True)
        self.tree.setIndentation(16)
        self.tree.setHeaderHidden(True)
        self.tree.setColumnHidden(1, True)
        self.tree.setColumnWidth(0, 155)
        self.tree.currentItemChanged.connect(
            lambda current, previous: (
                self.select_object(current.data(0, Qt.ItemDataRole.UserRole))
                if current
                else None
            )
        )
        self.tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._tree_menu)
        objects_layout.addWidget(self.tree, 1)
        self.object_count = QLabel()
        self.object_count.setObjectName("muted")
        objects_layout.addWidget(self.object_count)
        self._dock(
            "objects", "Browser", objects, Qt.DockWidgetArea.LeftDockWidgetArea
        ).setMinimumWidth(220)

        inspector = QWidget()
        inspector_layout = QVBoxLayout(inspector)
        inspector_layout.setContentsMargins(8, 8, 8, 6)
        inspector_layout.setSpacing(4)
        self.inspector_title = QLabel("Project settings")
        self.inspector_title.setObjectName("inspector_title")
        inspector_layout.addWidget(self.inspector_title)
        self.inspector_hint = QLabel()
        self.inspector_hint.setWordWrap(True)
        self.inspector_hint.setObjectName("muted")
        inspector_layout.addWidget(self.inspector_hint)
        self.properties = QWidget()
        self.properties.setObjectName("property_fields")
        self.form = QFormLayout(self.properties)
        self.form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        self.form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapAllRows)
        self.form.setVerticalSpacing(4)
        self.form.setContentsMargins(0, 2, 0, 8)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.properties)
        inspector_layout.addWidget(scroll, 1)
        self._dock(
            "inspector",
            "Inspector · SI",
            inspector,
            Qt.DockWidgetArea.RightDockWidgetArea,
        ).setMinimumWidth(280)

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setCentralWidget(central)
        top = QHBoxLayout()
        top.setSpacing(4)
        layout.addLayout(top)
        self.mode = QComboBox()
        self.mode.setAccessibleName("Workspace")
        self.mode.addItems(["Design", "Calculation result"])
        self.mode.model().item(1).setEnabled(False)
        self.mode.currentIndexChanged.connect(self._mode_changed)
        top.addWidget(self.mode)
        self.scene_label = QLabel("3D scene · metres")
        self.scene_label.setObjectName("muted")
        top.addWidget(self.scene_label)
        top.addStretch()
        camera_menu = QMenu(self)
        for direction, label, shortcut in (
            ("iso", "Isometric view", "0"),
            ("front", "Front view", "1"),
            ("side", "Side view", "3"),
            ("top", "Top view", "7"),
        ):
            action = self._action(
                "camera_" + direction,
                label,
                lambda checked=False, d=direction: self.viewport.camera(d),
                shortcut,
            )
            camera_menu.addAction(action)
            view.addAction(action)
        fit = self._action(
            "fit",
            "Fit selection",
            lambda: self.viewport.camera("iso", selection=True),
            "F",
        )
        camera_menu.addAction(fit)
        view.addAction(fit)
        camera_button = QToolButton()
        camera_button.setText("Views")
        camera_button.setMenu(camera_menu)
        camera_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        top.addWidget(camera_button)
        self.viewport = Viewport()
        camera_menu.addSeparator()
        for key, label, parallel in (
            ("orthographic", "Orthographic projection", True),
            ("perspective", "Perspective projection", False),
        ):
            action = self._action(
                key,
                label,
                lambda checked=False, p=parallel: self.viewport.set_parallel_projection(
                    p
                ),
            )
            camera_menu.addAction(action)
            view.addAction(action)
        grid = self._action("grid", "XY reference grid", self.viewport.set_grid_visible)
        grid.setCheckable(True)
        grid.setChecked(True)
        grid.setIcon(line_icon("grid"))
        view.addAction(grid)
        grid_button = QToolButton()
        grid_button.setDefaultAction(grid)
        top.addWidget(grid_button)
        layout.addWidget(self.viewport, 1)
        self.navigation_hint = QLabel(
            "Orbit: drag · Pan: Shift + drag · Zoom: scroll · F: fit"
        )
        self.navigation_hint.setObjectName("muted")
        self.navigation_hint.setWordWrap(True)
        layout.addWidget(self.navigation_hint)
        self.transform_mode = QComboBox()
        self.transform_mode.setAccessibleName("3D manipulation mode")
        self.transform_mode.addItems(["Select", "Move", "Rotate", "Move and rotate"])
        self.transform_mode.setCurrentIndex(0)
        self.transform_mode.currentIndexChanged.connect(
            self.viewport.set_transform_mode
        )
        self.transform_mode.hide()
        self.transform_buttons = []
        from PySide6.QtWidgets import QButtonGroup

        self.transform_group = QButtonGroup(self)
        for i, (name, label) in enumerate(
            (("select", "Select"), ("move", "Move"), ("rotate", "Rotate"))
        ):
            button = QToolButton()
            button.setIcon(line_icon(name))
            button.setIconSize(QSize(20, 20))
            button.setToolTip(label)
            button.setAccessibleName(label)
            button.setCheckable(True)
            button.setChecked(i == 0)
            self.transform_group.addButton(button, i)
            self.transform_buttons.append(button)
            top.addWidget(button)
        self.transform_group.idClicked.connect(self.transform_mode.setCurrentIndex)
        self.viewport.set_transform_mode(0)
        for key, label, slot in (
            ("hide", "Toggle selection visibility", self.toggle_visibility),
            ("isolate", "Isolate selection", self.isolate_selection),
            ("show_all", "Show all objects", self.show_all_objects),
        ):
            view.addAction(self._action(key, label, slot))

        setup = QWidget()
        setup_layout = QVBoxLayout(setup)
        self.validation_label = QLabel()
        self.validation_label.setWordWrap(True)
        self.validation_label.setObjectName("badge")
        setup_layout.addWidget(self.validation_label)
        settings_form = QFormLayout()
        setup_layout.addLayout(settings_form)
        self.duration = QLineEdit("2.0")
        self.step = QLineEdit("0.005")
        self.duration.setAccessibleName("Calculation duration in seconds")
        self.step.setAccessibleName("Calculation time step in seconds")
        settings_form.addRow("Duration [s]", self.duration)
        settings_form.addRow("Time step [s]", self.step)
        for field in (self.duration, self.step):
            field.textEdited.connect(self._workspace_state)
        row = QHBoxLayout()
        setup_layout.addLayout(row)
        self.run_button = QPushButton("Run")
        self.run_button.setObjectName("primary")
        self.run_button.clicked.connect(self.run)
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop)
        self.stop_button.setEnabled(False)
        row.addWidget(self.run_button, 1)
        row.addWidget(self.stop_button)
        self.run_state = QLabel("Ready to run")
        self.run_state.setObjectName("muted")
        self.run_state.setWordWrap(True)
        setup_layout.addWidget(self.run_state)
        setup_layout.addStretch()
        self._dock(
            "analysis",
            "Run · dynamics",
            setup,
            Qt.DockWidgetArea.LeftDockWidgetArea,
        )
        analysis.addAction(self._action("run", "Run design", self.run, "Ctrl+Return"))
        analysis.addAction(self._action("stop", "Stop calculation", self.stop))

        results = QWidget()
        results_layout = QVBoxLayout(results)
        run_row = QHBoxLayout()
        results_layout.addLayout(run_row)
        self.run_combo = QComboBox()
        self.run_combo.setAccessibleName("Displayed run")
        self.run_combo.setMinimumContentsLength(15)
        self.run_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.run_combo.currentIndexChanged.connect(self.choose_run)
        run_row.addWidget(QLabel("Run"))
        run_row.addWidget(self.run_combo, 1)
        self.compare_combo = QComboBox()
        self.compare_combo.setAccessibleName("Reference run for comparison")
        self.compare_combo.setMinimumContentsLength(15)
        self.compare_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.compare_combo.addItem("No comparison", None)
        self.compare_combo.currentIndexChanged.connect(self.update_comparison)
        run_row.addWidget(QLabel("Compare with"))
        run_row.addWidget(self.compare_combo, 1)
        self.result_label = QLabel(
            "No result. Run a calculation to inspect the motion."
        )
        self.result_label.setWordWrap(True)
        self.result_label.setTextFormat(Qt.TextFormat.PlainText)
        results_layout.addWidget(self.result_label)
        self.comparison_label = QLabel()
        self.comparison_label.setObjectName("muted")
        self.comparison_label.setWordWrap(True)
        self.comparison_label.setTextFormat(Qt.TextFormat.PlainText)
        self.comparison_label.hide()
        results_layout.addWidget(self.comparison_label)
        tools = QHBoxLayout()
        results_layout.addLayout(tools)
        self.series_combo = QComboBox()
        self.series_combo.setMinimumContentsLength(16)
        self.series_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.series_combo.setAccessibleName("Displayed physical quantity")
        self.series_combo.currentIndexChanged.connect(self._select_series)
        tools.addWidget(self.series_combo, 1)
        self.export_button = QPushButton("Export CSV")
        self.export_button.clicked.connect(self.export_dialog)
        self.provenance_button = QPushButton("Provenance")
        self.provenance_button.clicked.connect(self.provenance)
        for button in (self.export_button, self.provenance_button):
            button.setEnabled(False)
            tools.addWidget(button)
        self.result_tabs = QTabWidget()
        self.curve = SeriesView()
        self.curve.sample_selected.connect(self.seek_sample)
        self.result_tabs.addTab(self.curve, "Curve")
        self.sample_table = QTableView()
        self.sample_table.setAccessibleName("Native samples of the displayed quantity")
        self.sample_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.sample_table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        self.sample_table.setAlternatingRowColors(True)
        self.sample_model = SampleTableModel(self)
        self.sample_table.setModel(self.sample_model)
        self.sample_table.horizontalHeader().setStretchLastSection(True)
        self.sample_table.selectionModel().currentRowChanged.connect(
            lambda current, previous: (
                self.seek_sample(current.row()) if current.isValid() else None
            )
        )
        self.result_tabs.addTab(self.sample_table, "Samples")
        results_layout.addWidget(self.result_tabs, 1)
        playback = QHBoxLayout()
        results_layout.addLayout(playback)
        self.play_button = QPushButton("Play")
        self.play_button.clicked.connect(self.toggle_play)
        self.play_button.setEnabled(False)
        playback.addWidget(self.play_button)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setAccessibleName("Time sample")
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self._frame)
        self.slider.sliderPressed.connect(self.pause)
        self.slider.actionTriggered.connect(lambda _: self.pause())
        playback.addWidget(self.slider, 1)
        self.speed = QComboBox()
        self.speed.setAccessibleName("Playback speed")
        for rate in (0.25, 0.5, 1.0, 2.0, 4.0):
            self.speed.addItem(f"×{rate:g}", rate)
        self.speed.setCurrentIndex(2)
        self.speed.currentIndexChanged.connect(self._speed_changed)
        playback.addWidget(self.speed)
        self.time_input = QDoubleSpinBox()
        self.time_input.setAccessibleName("Displayed time in seconds")
        self.time_input.setDecimals(6)
        self.time_input.setSuffix(" s")
        self.time_input.setKeyboardTracking(False)
        self.time_input.setEnabled(False)
        self.time_input.valueChanged.connect(self.seek_time)
        playback.addWidget(self.time_input)
        self.time_label = QLabel("t = —")
        self.time_label.setObjectName("muted")
        playback.addWidget(self.time_label)
        self._dock(
            "results",
            "Results · native samples",
            results,
            Qt.DockWidgetArea.BottomDockWidgetArea,
        )

        self.diagnostics = QTreeWidget()
        self.diagnostics.setHeaderLabels(["Object", "Diagnostic and required action"])
        self.diagnostics.setAccessibleName("Model diagnostics")
        self.diagnostics.itemActivated.connect(self._activate_diagnostic)
        self._dock(
            "diagnostics",
            "Diagnostics",
            self.diagnostics,
            Qt.DockWidgetArea.BottomDockWidgetArea,
        )
        self.tabifyDockWidget(self.docks["results"], self.docks["diagnostics"])
        self.docks["results"].raise_()
        panels = view.addMenu("Panels")
        for key, dock in self.docks.items():
            action = dock.toggleViewAction()
            panels.addAction(action)
            self.commands["panel_" + key] = action
        view.addAction(
            self._action("reset_layout", "Restore layout", self.reset_workspace)
        )
        themes = view.addMenu("Theme")
        self.theme_actions = QActionGroup(self)
        for name, label in (("dark", "Dark"), ("light", "Light")):
            action = self._action(
                "theme_" + name,
                f"Theme {label.lower()}",
                lambda checked=False, n=name: self.set_theme(n),
            )
            action.setCheckable(True)
            self.theme_actions.addAction(action)
            themes.addAction(action)
        help_menu.addAction(
            self._action(
                "navigation_help", "3D navigation help", self.navigation_help, "F1"
            )
        )
        self.status = QLabel("Ready. Select an object to edit.")
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        self.status.setWordWrap(True)
        self.statusBar().addWidget(self.status, 1)
        self.units_label = QLabel("SI · m · kg · s")
        self.statusBar().addPermanentWidget(self.units_label)
        self.docks["analysis"].hide()
        self.docks["results"].hide()
        self.docks["diagnostics"].hide()
        self._default_layout = self.saveState(3)
        self.workspace_tabs.currentChanged.connect(self.switch_workspace)
        self.set_theme("dark")
        self.resizeDocks(
            [self.docks["objects"], self.docks["inspector"]],
            [250, 300],
            Qt.Orientation.Horizontal,
        )
        if self.settings is not None:
            self.set_theme(self.settings.value("theme", "dark"))
            if self.settings.contains("geometry"):
                self.restoreGeometry(self.settings.value("geometry"))
            if self.settings.contains("layout_v3"):
                self.restoreState(self.settings.value("layout_v3"), 3)

    def switch_workspace(self, index):
        if index == 2:
            self.mode.setCurrentIndex(1)
            if self.mode.currentIndex() != 1:
                self.workspace_tabs.blockSignals(True)
                self.workspace_tabs.setCurrentIndex(0)
                self.workspace_tabs.blockSignals(False)
                return
        else:
            self.pause()
            self.mode.setCurrentIndex(0)
        self.viewport.setMinimumHeight(210 if index == 2 else 300)
        QTimer.singleShot(0, self.viewport.fit_scene)
        self.docks["analysis"].setVisible(index == 1)
        self.docks["results"].setVisible(index == 2)
        self.docks["diagnostics"].setVisible(index == 1)
        if index == 2:
            self.docks["results"].raise_()
        elif index == 1:
            self.docks["diagnostics"].raise_()
        self.resizeDocks([self.docks["results"]], [285], Qt.Orientation.Vertical)

    def show_commands(self):
        palette = CommandPalette(self.commands.values(), self)
        palette.exec()

    def show_context_tools(self):
        keys = (
            ("fit", "hide", "isolate", "duplicate", "delete", "joint", "load")
            if self.object(display=True)
            else ("box", "cylinder", "sphere", "show_all")
        )
        palette = CommandPalette([self.commands[key] for key in keys], self)
        palette.setWindowTitle("Selection tools")
        palette.exec()

    def set_theme(self, name):
        self.theme_name = apply_theme(self, name)
        self.commands["theme_" + self.theme_name].setChecked(True)
        self.curve.set_theme(self.theme_name)
        color = self.palette().text().color().name()
        for key, action in self.commands.items():
            if key in PATHS:
                action.setIcon(line_icon(key, color))
        self.commands["commands"].setIcon(line_icon("search", color))
        for g in range(self.tree.topLevelItemCount()):
            group = self.tree.topLevelItem(g)
            for i in range(group.childCount()):
                item = group.child(i)
                obj = self.object(item.data(0, Qt.ItemDataRole.UserRole), display=True)
                if obj is not None:
                    name = (
                        obj.shape
                        if isinstance(obj, Body)
                        else "joint"
                        if hasattr(obj, "kind")
                        else "load"
                    )
                    item.setIcon(0, line_icon(name, color))
        for button, name in zip(self.transform_buttons, ("select", "move", "rotate")):
            button.setIcon(line_icon(name, color))

    def reset_workspace(self):
        self.restoreState(self._default_layout, 3)
        self.switch_workspace(self.workspace_tabs.currentIndex())
        self.resizeDocks(
            [self.docks["objects"], self.docks["inspector"]],
            [250, 300],
            Qt.Orientation.Horizontal,
        )
        self.resizeDocks([self.docks["results"]], [285], Qt.Orientation.Vertical)

    def save_workspace(self):
        if self.settings is not None:
            self.settings.setValue("geometry", self.saveGeometry())
            self.settings.setValue("layout_v3", self.saveState(3))
            self.settings.setValue("theme", self.theme_name)
            self.settings.sync()

    def _filter_tree(self, *_):
        terms = search_key(self.tree_search.text()).split()
        visible = 0
        for g in range(self.tree.topLevelItemCount()):
            group = self.tree.topLevelItem(g)
            children_visible = 0
            for i in range(group.childCount()):
                item = group.child(i)
                text = search_key(
                    item.text(0)
                    + " "
                    + item.text(1)
                    + " "
                    + str(item.data(0, Qt.ItemDataRole.UserRole))
                )
                show = all(t in text for t in terms)
                item.setHidden(not show)
                children_visible += show
            group.setHidden(bool(terms) and not children_visible)
            visible += children_visible
        self.object_count.setText(
            f"{visible} object{'s' if visible != 1 else ''} shown in the browser"
        )

    def _tree_menu(self, point):
        item = self.tree.itemAt(point)
        if item:
            self.select_object(item.data(0, Qt.ItemDataRole.UserRole))
        menu = QMenu(self)
        for key in ("fit", "hide", "isolate", "show_all", "duplicate", "delete"):
            menu.addAction(self.commands[key])
        menu.exec(self.tree.viewport().mapToGlobal(point))

    def toggle_visibility(self):
        if self.selection:
            self.viewport.set_object_visible(
                self.selection, self.selection in self.viewport.hidden_objects
            )
            self._workspace_state()

    def isolate_selection(self):
        if self.selection:
            self.viewport.isolate(self.selection)
            self._workspace_state()

    def show_all_objects(self):
        self.viewport.show_all()
        self._workspace_state()

    def _activate_diagnostic(self, item, column):
        self.mode.setCurrentIndex(0)
        self.tree_search.clear()
        self.select_object(item.data(0, Qt.ItemDataRole.UserRole))
        self.docks["inspector"].show()
        self.docks["inspector"].raise_()

    def _workspace_state(self, *_):
        settings_pending = (self.duration.text(), self.step.text()) != (
            repr(self.project.duration),
            repr(self.project.step),
        )
        pending = self.dirty_fields or settings_pending
        issues = self.project.diagnostics()
        self.validation_label.setText(
            "Draft · unapplied properties"
            if pending
            else f"Invalid model · {len(issues)} diagnostics"
            if issues
            else "Inputs checked · ready to run"
        )
        self.project_label.setText(
            f"{self.project.name}\nRevision {self.project.revision}"
        )
        modified = self.project != self._saved or pending
        self.setWindowModified(modified)
        self.setWindowTitle(f"{self.project.name}[*] — Vinkulum Studio")
        result_mode = self.mode.currentIndex() == 1
        self.scene_label.setText(
            "Result · read only" if result_mode else "Design · metres"
        )
        self.transform_mode.setEnabled(not result_mode)
        for button in self.transform_buttons:
            button.setEnabled(not result_mode)
        self.workspace_tabs.setTabEnabled(2, self.result is not None)
        if result_mode and self.workspace_tabs.currentIndex() != 2:
            self.workspace_tabs.setCurrentIndex(2)
        elif not result_mode and self.workspace_tabs.currentIndex() == 2:
            self.workspace_tabs.setCurrentIndex(0)
        obj = self.object(display=True)
        self.inspector_title.setText(obj.name if obj else "Project settings")
        self.inspector_hint.setText(
            "Run snapshot · read only"
            if result_mode
            else "Pending property changes"
            if self.dirty_fields
            else "Document properties · SI units"
        )
        self.commands["undo"].setEnabled(not result_mode and bool(self.history.past))
        self.commands["redo"].setEnabled(not result_mode and bool(self.history.future))
        for key in ("duplicate", "delete"):
            self.commands[key].setEnabled(not result_mode and obj is not None)
        self.commands["load"].setEnabled(not result_mode and isinstance(obj, Body))
        self.commands["joint"].setEnabled(not result_mode and bool(self.project.bodies))
        for key in ("hide", "isolate", "fit"):
            self.commands[key].setEnabled(obj is not None)

    def _update_run_choices(self, selected_id):
        reference = self.compare_combo.currentData()
        for combo in (self.run_combo, self.compare_combo):
            combo.blockSignals(True)
            combo.clear()
        self.compare_combo.addItem("No comparison", None)
        for number, result in enumerate(self.run_archive.results, 1):
            label = f"{result.project.name} · r{result.project.revision} · {result.run_id[:8]}"
            self.run_combo.addItem(label, result.run_id)
            self.compare_combo.addItem(label, result.run_id)
        self.run_combo.setCurrentIndex(self.run_combo.findData(selected_id))
        self.compare_combo.setCurrentIndex(
            max(0, self.compare_combo.findData(reference))
        )
        for combo in (self.run_combo, self.compare_combo):
            combo.blockSignals(False)

    def choose_run(self, index):
        result = self.run_archive.find(self.run_combo.itemData(index))
        if result is not None:
            if not self.apply_properties():
                self.run_combo.blockSignals(True)
                self.run_combo.setCurrentIndex(
                    self.run_combo.findData(self.result.run_id)
                )
                self.run_combo.blockSignals(False)
                return
            self._display_result(result)
            self.mode.setCurrentIndex(1)

    def update_comparison(self, *_):
        self.curve.set_reference(None)
        self.comparison_label.hide()
        reference = self.run_archive.find(self.compare_combo.currentData())
        index = self.series_combo.currentIndex()
        if reference is None or self.result is None or index < 0:
            return
        self.comparison_label.show()
        current = self.result
        differences = model_differences(current.project, reference.project)
        import json

        current_manifest, reference_manifest = (
            json.loads(current.manifest_json),
            json.loads(reference.manifest_json),
        )
        for field, label in (
            ("method", "numerical method"),
            ("rho_infinity", "numerical dissipation"),
            ("kernel_sha256", "kernel binary"),
        ):
            if current_manifest.get(field) != reference_manifest.get(field):
                differences.append(label)
        description = (
            "Differences: " + "; ".join(differences)
            if differences
            else "Same model and numerical settings."
        )
        key = series_keys(current.project)[index]
        reference_keys = series_keys(reference.project)
        if key not in reference_keys:
            self.comparison_label.setText(
                description
                + " · Quantity missing from the reference: different object identity."
            )
            return
        label, unit, values = reference.series()[reference_keys.index(key)]
        if unit != self._series[index][1]:
            self.comparison_label.setText("Comparison unavailable: incompatible units.")
            return
        self.curve.set_reference((reference.time, values, reference.run_id[:8]))
        frame = (
            "World frame"
            if key[1] in {"position", "velocity"}
            else "Relative joint coordinate"
        )
        self.comparison_label.setText(
            description
            + f" · Solid line: {current.run_id[:8]}; dashed line: {reference.run_id[:8]}. {frame}; each run keeps its own timestamps."
        )

    def seek_sample(self, index):
        if self.result is not None and 0 <= index < len(self.result.time):
            self.pause()
            self.mode.setCurrentIndex(1)
            self.slider.setValue(index)

    def seek_time(self, value):
        if self.result is not None:
            import numpy as np

            right = min(
                int(np.searchsorted(self.result.time, value)), len(self.result.time) - 1
            )
            left = max(0, right - 1)
            index = (
                left
                if abs(self.result.time[left] - value)
                <= abs(self.result.time[right] - value)
                else right
            )
            self.seek_sample(index)
            self._frame(index)

    def _speed_changed(self, *_):
        if self.timer.isActive():
            self.pause()
            self.toggle_play()

    def navigation_help(self):
        from PySide6.QtWidgets import QMessageBox

        QMessageBox.information(
            self,
            "Navigation and editing",
            "Drag: orbit around the model.\nShift + drag: pan the camera. Scroll: zoom.\n0 / 1 / 3 / 7: isometric / front / side / top views. F: fit selection.\nChoose Select, Move or Rotate to set the manipulator.\nPositions are also editable in the inspector, in metres and degrees.\n\nMoving a body does not solve constraints. The grid is not a contact surface.\nOrange forces and purple moments: world directions, symbolic lengths.\nResults display native samples; their accuracy is not certified.",
        )
