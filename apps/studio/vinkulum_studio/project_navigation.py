"""One project window with persistent, lazily created analysis views.

The host changes presentation only. Each calculation keeps its captured inputs,
validation and close guards; switching selection never starts or resets a solver.
"""

from PySide6.QtCore import QSignalBlocker, Qt
from PySide6.QtWidgets import QDockWidget, QMenu, QStackedWidget, QTreeWidgetItem

from .theme import apply_theme


class ProjectNavigation:
    def initialize_project_pages(self, model_page):
        self.project_pages = QStackedWidget(self)
        self.project_pages.setObjectName("project_views")
        self.model_page = model_page
        self.project_pages.addWidget(model_page)
        self.setCentralWidget(self.project_pages)

    def active_viewport(self):
        entry = self.analysis_pages.get(self.active_analysis)
        return entry[0].viewport if entry else self.viewport

    def register_analysis(self, key, page, label, attribute=None):
        self.analysis_pages[key] = (page, label)
        self.project_pages.addWidget(page)
        apply_theme(page, self.theme_name)
        # Embedded panels stay attached to the project, including on macOS.
        for dock in page.findChildren(QDockWidget):
            dock.setFeatures(QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        page.closed.connect(lambda: self.remove_analysis(key, attribute))
        save = getattr(page, "save_action", None)
        if save is not None:
            save.changed.connect(self.sync_context_controls)
        for action_name in ("undo_action", "redo_action"):
            action = getattr(page, action_name, None)
            if action is not None:
                action.setShortcut("")
                action.changed.connect(self.sync_context_controls)
        self.sync_analysis_browser()

    def remove_analysis(self, key, attribute=None):
        entry = self.analysis_pages.pop(key, None)
        if attribute:
            setattr(self, attribute, None)
        if self.active_analysis == key:
            self.show_model_context()
        if entry:
            self.project_pages.removeWidget(entry[0])
        self.sync_analysis_browser()

    def activate_analysis(self, key):
        if key in ("native", "result"):
            return self.show_model_context(1 if key == "native" else 2)
        entry = self.analysis_pages.get(key)
        if entry is None:
            return False
        if self.active_analysis is None and not self.apply_properties():
            self.sync_analysis_browser()
            return False
        self.pause()
        self.active_analysis = key
        self.project_pages.setCurrentWidget(entry[0])
        self._workspace_state()
        return True

    def show_model_context(self, index=0):
        if index == 2 and self.result is None:
            return False
        if not self.apply_properties():
            self.sync_analysis_browser()
            return False
        self.active_analysis = None
        self.project_pages.setCurrentWidget(self.model_page)
        self.workspace_tabs.setCurrentIndex(index)
        self.switch_workspace(index)
        self._workspace_state()
        if index == 0:
            with QSignalBlocker(self.tree):
                self.tree.setCurrentItem(None)
                for g in range(self.tree.topLevelItemCount()):
                    group = self.tree.topLevelItem(g)
                    for i in range(group.childCount()):
                        item = group.child(i)
                        if item.data(0, Qt.ItemDataRole.UserRole) == self.selection:
                            self.tree.setCurrentItem(item)
        return True

    def activate_browser_item(self, item):
        identifier = item.data(0, Qt.ItemDataRole.UserRole) if item else None
        if isinstance(identifier, str) and identifier.startswith("analysis:"):
            self.activate_analysis(identifier.removeprefix("analysis:"))
        elif self.show_model_context():
            self.select_object(identifier)

    def sync_analysis_browser(self):
        with QSignalBlocker(self.tree):
            group = next(
                (
                    self.tree.topLevelItem(i)
                    for i in range(self.tree.topLevelItemCount())
                    if self.tree.topLevelItem(i).text(0) == "Analyses"
                ),
                None,
            )
            if group is None:
                group = QTreeWidgetItem(["Analyses"])
                self.tree.addTopLevelItem(group)
                group.setExpanded(True)
            entries = [("native", "Motion analysis")]
            if self.result is not None:
                entries.append(("result", "Motion results"))
            entries.extend(
                (key, value[1]) for key, value in self.analysis_pages.items()
            )
            wanted = {"analysis:" + key for key, _ in entries}
            for i in reversed(range(group.childCount())):
                if group.child(i).data(0, Qt.ItemDataRole.UserRole) not in wanted:
                    group.takeChild(i)
            for key, label in entries:
                identifier = "analysis:" + key
                item = next(
                    (
                        group.child(i)
                        for i in range(group.childCount())
                        if group.child(i).data(0, Qt.ItemDataRole.UserRole)
                        == identifier
                    ),
                    None,
                )
                if item is None:
                    item = QTreeWidgetItem([label, "Analysis"])
                    item.setData(0, Qt.ItemDataRole.UserRole, identifier)
                    group.addChild(item)
                item.setToolTip(
                    0, "Select to show parameters and results in this window"
                )
                current = self.active_analysis or {1: "native", 2: "result"}.get(
                    self.workspace_tabs.currentIndex()
                )
                if key == current:
                    self.tree.setCurrentItem(item)
        self._filter_tree()

    def sync_context_controls(self):
        if not hasattr(self, "status"):
            return
        embedded = self.active_analysis is not None
        index = self.workspace_tabs.currentIndex()
        context = (self.active_analysis, 0 if embedded else index)
        previous = getattr(self, "_presentation_context", (None, 0))
        if context != previous:
            keys = ("inspector", "analysis", "diagnostics", "results")
            if previous is not None and previous[0] is None:
                self._native_layouts[previous[1]] = {
                    key: not self.docks[key].isHidden() for key in keys
                }
                self._native_layout_state = self.saveState(3)
            defaults = {
                "inspector": index != 1,
                "analysis": index == 1,
                "diagnostics": index == 1,
                "results": index == 2,
            }
            visibility = self._native_layouts.get(index, defaults)
            for key in keys:
                self.docks[key].setVisible(not embedded and visibility[key])
            self._presentation_context = context
        self.statusBar().setVisible(not embedded)
        for key in ("inspector", "analysis", "diagnostics", "results"):
            self.commands["panel_" + key].setEnabled(not embedded)
        if embedded:
            for action in self._design_actions:
                action.setEnabled(False)
            for key in ("hide", "isolate", "show_all", "context_tools", "run", "stop"):
                self.commands[key].setEnabled(False)
            self.commands["fit"].setEnabled(True)
        else:
            for key in ("show_all", "context_tools"):
                self.commands[key].setEnabled(True)
        entry = self.analysis_pages.get(self.active_analysis)
        action = self.commands["save"]
        action.setText("Save analysis…" if entry else "Save project…")
        action.setToolTip(
            "Save the displayed analysis input to its own file"
            if entry
            else "Save the mechanism and motion settings"
        )
        if entry:
            page = entry[0]
            save = getattr(page, "save_action", None) or getattr(page, "save", None)
            action.setEnabled(save is not None and save.isEnabled())
            for key in ("undo", "redo"):
                edit = getattr(page, key + "_action", None)
                self.commands[key].setEnabled(edit is not None and edit.isEnabled())
        else:
            action.setEnabled(True)
        if not embedded:
            busy = self.controller.process is not None
            self.commands["run"].setEnabled(not busy and not self.project_diagnostics())
            self.commands["stop"].setEnabled(busy)

    def undo_current(self):
        if self.active_analysis is None:
            self.undo()
        else:
            action = getattr(
                self.analysis_pages[self.active_analysis][0], "undo_action", None
            )
            if action is not None:
                action.trigger()

    def redo_current(self):
        if self.active_analysis is None:
            self.redo()
        else:
            action = getattr(
                self.analysis_pages[self.active_analysis][0], "redo_action", None
            )
            if action is not None:
                action.trigger()

    def save_current(self):
        entry = self.analysis_pages.get(self.active_analysis)
        if entry:
            page = entry[0]
            action = getattr(page, "save_action", None)
            if action is not None:
                action.trigger()
            else:
                page.save.click()
        else:
            return self.save_dialog()

    def analysis_menu(self, key, position):
        menu = QMenu(self)
        menu.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        menu.addAction("Show analysis", lambda: self.activate_analysis(key))
        if key in self.analysis_pages:
            page = self.analysis_pages[key][0]
            menu.addSeparator()
            menu.addAction("Close analysis…", page.close)
        menu.popup(position)
