"""Application colours shared by Qt controls and scientific views."""

from PySide6.QtGui import QColor, QPalette

THEMES = {
    "dark": dict(
        bg="#17191e",
        panel="#202329",
        field="#292d35",
        text="#eef0f6",
        muted="#a4abb9",
        border="#393e49",
        hover="#343d52",
        accent="#a2b6ff",
        accent_text="#17203b",
        grid="#393e49",
        curve="#a2b6ff",
        cursor="#ffca80",
    ),
    "light": dict(
        bg="#edf1f6",
        panel="#ffffff",
        field="#f6f8fb",
        text="#182c42",
        muted="#50637a",
        border="#b4c1cf",
        hover="#e0eaf4",
        accent="#086f65",
        accent_text="#ffffff",
        grid="#d5dfe9",
        curve="#086f65",
        cursor="#9a5100",
    ),
}


def apply_theme(window, name):
    name = name if name in THEMES else "dark"
    c = THEMES[name]
    palette = QPalette(window.palette())
    for role, key in (
        ("Window", "bg"),
        ("WindowText", "text"),
        ("Base", "field"),
        ("AlternateBase", "panel"),
        ("Text", "text"),
        ("Button", "panel"),
        ("ButtonText", "text"),
        ("Highlight", "accent"),
        ("HighlightedText", "accent_text"),
        ("ToolTipBase", "panel"),
        ("ToolTipText", "text"),
        ("PlaceholderText", "muted"),
    ):
        palette.setColor(getattr(QPalette.ColorRole, role), QColor(c[key]))
    for role in (
        QPalette.ColorRole.Text,
        QPalette.ColorRole.ButtonText,
        QPalette.ColorRole.WindowText,
    ):
        palette.setColor(QPalette.ColorGroup.Disabled, role, QColor(c["muted"]))
    window.setPalette(palette)
    window.setStyleSheet(
        """
        QWidget { color: %(text)s; }
        QLabel { background: transparent; }
        QWidget#property_fields { background: %(panel)s; }
        QMenuBar { background: %(panel)s; color: %(text)s; }
        QMenuBar::item:selected { background: %(hover)s; }
        QMainWindow, QDialog { background: %(bg)s; color: %(text)s; }
        QToolButton { border: 0; }
        QToolBar { background: %(panel)s; border: 0; spacing: 4px; padding: 6px 8px; }
        QToolBar::separator { background: %(border)s; width: 1px; margin: 5px; }
        QDockWidget { font-weight: 600; color: %(text)s; }
        QDockWidget::title { background: %(panel)s; padding: 8px; border-bottom: 1px solid %(border)s; }
        QMainWindow::separator { background: %(bg)s; width: 4px; height: 4px; }
        QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox, QTextEdit {
            background: %(field)s; color: %(text)s; border: 1px solid transparent;
            border-radius: 3px; padding: 4px 6px; min-height: 20px;
        }
        QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus, QSpinBox:focus,
        QTreeWidget:focus, QTableView:focus, QListWidget:focus {
            border: 1px solid %(accent)s;
        }
        QPushButton, QToolButton { color: %(text)s; background: %(panel)s;
            border: 1px solid transparent; border-radius: 3px; padding: 4px 8px; min-height: 20px; }
        QToolButton::menu-indicator { subcontrol-position: right center; right: 3px; }
        QToolButton:has-menu { padding-right: 16px; }
        QPushButton:hover, QToolButton:hover { background: %(hover)s; }
        QPushButton:focus, QToolButton:focus { border: 1px solid %(accent)s; }
        QPushButton:checked, QToolButton:checked { background: %(hover)s; border-color: %(accent)s; }
        QPushButton:disabled, QToolButton:disabled { color: %(muted)s; border-color: %(bg)s; }
        QPushButton#primary { background: %(accent)s; color: %(accent_text)s; font-weight: 600; }
        QPushButton#primary:disabled { background: %(hover)s; color: %(muted)s; }
        QTreeWidget, QListWidget, QTableView { background: %(panel)s; color: %(text)s;
            border: 1px solid transparent; border-radius: 4px; alternate-background-color: %(field)s; }
        QTreeWidget::item, QListWidget::item { padding: 4px 3px; }
        QTreeWidget::item:selected, QListWidget::item:selected { background: %(hover)s; color: %(text)s; }
        QHeaderView { background: %(panel)s; }
        QHeaderView::section { background: %(panel)s; color: %(muted)s; padding: 5px; border: 0; }
        QScrollArea { border: 0; }
        QTabWidget::pane { border: 1px solid %(border)s; }
        QTabBar::tab { color: %(muted)s; background: %(panel)s; padding: 6px 12px; }
        QTabBar::tab:selected { color: %(text)s; border-bottom: 2px solid %(accent)s; }
        QLabel#brand { color: %(text)s; font-weight: 600; }
        QLabel#section { color: %(muted)s; font-weight: 600; padding-top: 10px; padding-bottom: 3px; }
        QLabel#inspector_title { color: %(text)s; font-weight: 600; padding: 2px 0; }
        QLabel#muted { color: %(muted)s; }
        QLabel#badge { color: %(text)s; background: %(hover)s; border-radius: 5px; padding: 6px; }
        QWidget#vector_cell { background: %(field)s; border-radius: 3px; }
        QLineEdit#component { border: 1px solid transparent; padding: 4px 2px; }
        QLineEdit#component:focus { border: 1px solid %(accent)s; }
        QLabel#axis_label { color: %(muted)s; font-size: 10px; }
        QScrollBar:vertical { background: transparent; width: 8px; margin: 0; }
        QScrollBar::handle:vertical { background: %(border)s; min-height: 25px; border-radius: 4px; }
        QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
        QScrollBar:horizontal { background: transparent; height: 8px; margin: 0; }
        QScrollBar::handle:horizontal { background: %(border)s; min-width: 25px; border-radius: 4px; }
        QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
        QComboBox::drop-down { border: 0; width: 18px; }
        QTabWidget::pane { border: 0; }
        QStatusBar { background: %(panel)s; color: %(muted)s; }
        QMenu { background: %(panel)s; color: %(text)s; border: 1px solid %(border)s; padding: 5px; }
        QMenu::item { padding: 7px 24px; }
        QMenu::item:selected { background: %(hover)s; }
        QMenu::item:disabled { color: %(muted)s; }
        QToolTip { background: %(panel)s; color: %(text)s; border: 1px solid %(border)s; padding: 5px; }
    """
        % c
    )
    return name
