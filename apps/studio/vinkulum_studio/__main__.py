import sys


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--cad-worker":
        from .cad_worker import main as cad_main

        del sys.argv[1]
        return cad_main()
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        from .worker import main as worker_main

        del sys.argv[1]
        return worker_main()
    if len(sys.argv) == 3 and sys.argv[1] == "--bundle-check":
        from .bundle_check import main as check_main

        return check_main(sys.argv[2])
    from PySide6.QtWidgets import QApplication

    from .editor import EditorWindow

    app = QApplication(sys.argv)
    app.setApplicationName("Vinkulum Studio")
    app.setStyle("Fusion")
    from PySide6.QtCore import QSettings

    window = EditorWindow(QSettings("Vinkulum", "Studio"))
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
