"""Keyboard-first command search over the same QActions used in menus."""

import unicodedata

from PySide6.QtCore import QEvent, Qt
from PySide6.QtWidgets import (
    QDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QVBoxLayout,
)


def search_key(text):
    return "".join(
        c
        for c in unicodedata.normalize("NFKD", text.casefold())
        if not unicodedata.combining(c)
    )


class CommandPalette(QDialog):
    def __init__(self, actions, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Find a command")
        self.resize(620, 440)
        self.actions = list(actions)
        layout = QVBoxLayout(self)
        self.query = QLineEdit()
        self.query.setPlaceholderText("What would you like to do?")
        self.query.setAccessibleName("Find a command")
        self.query.installEventFilter(self)
        layout.addWidget(self.query)
        self.results = QListWidget()
        self.results.setAccessibleName("Available commands")
        layout.addWidget(self.results)
        self.hint = QLabel("↑ ↓ Choose    Enter Run    Esc Close")
        self.hint.setObjectName("muted")
        layout.addWidget(self.hint)
        self.query.textChanged.connect(self._filter)
        self.query.returnPressed.connect(self._execute)
        self.results.itemActivated.connect(self._execute)
        self._filter("")

    def _filter(self, query):
        self.results.clear()
        terms = search_key(query).split()
        for action in self.actions:
            label = action.text().replace("&", "")
            if not all(
                term in search_key(label + " " + action.toolTip()) for term in terms
            ):
                continue
            shortcut = action.shortcut().toString(
                action.shortcut().SequenceFormat.NativeText
            )
            item = QListWidgetItem(label + (f"    {shortcut}" if shortcut else ""))
            item.setData(Qt.ItemDataRole.UserRole, action)
            item.setToolTip(action.toolTip())
            if not action.isEnabled():
                item.setText(item.text() + " — unavailable")
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self.results.addItem(item)
        for i in range(self.results.count()):
            if self.results.item(i).flags() & Qt.ItemFlag.ItemIsEnabled:
                self.results.setCurrentRow(i)
                break
        self.hint.setText(
            "↑ ↓ Choose    Enter Run    Esc Close"
            if self.results.count()
            else "No matching command. Try another word."
        )

    def _execute(self, *_):
        item = self.results.currentItem()
        if item is None:
            return
        action = item.data(Qt.ItemDataRole.UserRole)
        if action.isEnabled():
            self.accept()
            action.trigger()

    def eventFilter(self, obj, event):
        if obj is self.query and event.type() == QEvent.Type.KeyPress:
            if event.key() in (Qt.Key.Key_Down, Qt.Key.Key_Up):
                direction = 1 if event.key() == Qt.Key.Key_Down else -1
                row = self.results.currentRow() + direction
                while 0 <= row < self.results.count():
                    if self.results.item(row).flags() & Qt.ItemFlag.ItemIsEnabled:
                        self.results.setCurrentRow(row)
                        break
                    row += direction
                return True
        return super().eventFilter(obj, event)
