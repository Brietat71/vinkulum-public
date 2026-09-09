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
        self.setWindowTitle("Rechercher une commande")
        self.resize(620, 440)
        self.actions = list(actions)
        layout = QVBoxLayout(self)
        self.query = QLineEdit()
        self.query.setPlaceholderText("Que souhaitez-vous faire ?")
        self.query.setAccessibleName("Rechercher une commande")
        self.query.installEventFilter(self)
        layout.addWidget(self.query)
        self.results = QListWidget()
        self.results.setAccessibleName("Commandes disponibles")
        layout.addWidget(self.results)
        self.hint = QLabel("↑ ↓ Choisir    Entrée Exécuter    Échap Fermer")
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
                item.setText(item.text() + " — indisponible")
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self.results.addItem(item)
        for i in range(self.results.count()):
            if self.results.item(i).flags() & Qt.ItemFlag.ItemIsEnabled:
                self.results.setCurrentRow(i)
                break
        self.hint.setText(
            "↑ ↓ Choisir    Entrée Exécuter    Échap Fermer"
            if self.results.count()
            else "Aucune commande correspondante. Essayez un autre mot."
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
