import sys
from PySide6.QtWidgets import QApplication
from .window import MainWindow


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("Vinkulum Studio")
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
