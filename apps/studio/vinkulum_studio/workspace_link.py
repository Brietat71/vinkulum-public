"""Visible return path between workspaces, retaining their edited inputs."""

from PySide6.QtGui import QAction, QKeySequence


def add_workspace_return(window, toolbar):
    parent = window.parentWidget()
    if parent is None:
        return

    def return_to_parent():
        parent.show()
        parent.raise_()
        parent.activateWindow()
        window.hide()

    action = QAction("← Back", window)
    action.setObjectName("back_to_workspace")
    action.setToolTip("Return to the previous workspace; keep this study open")
    action.setShortcut(QKeySequence(QKeySequence.StandardKey.Back))
    action.triggered.connect(return_to_parent)
    window.addAction(action)
    toolbar.addAction(action)
    toolbar.addSeparator()
