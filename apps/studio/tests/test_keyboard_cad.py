"""Issue #2: real CAD editing with keyboard input, also runnable at 200% DPI."""

import importlib.util
import os
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import Qt, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog, QLineEdit


@unittest.skipUnless(
    os.environ.get("VINKULUM_KEYBOARD_TESTS") == "1"
    and importlib.util.find_spec("build123d"),
    "Desktop and CAD required",
)
class KeyboardCadRecipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def key(self, key, modifiers=Qt.KeyboardModifier.NoModifier):
        focus = self.app.focusWidget()
        self.assertIsNotNone(focus, "Keyboard focus was lost")
        QTest.keyClick(focus, key, modifiers)
        QTest.qWait(10)

    def enter(self, text):
        self.key(Qt.Key.Key_A, Qt.KeyboardModifier.ControlModifier)
        QTest.keyClicks(self.app.focusWidget(), text)

    def tab_to(self, target, backwards=False):
        sequence = []
        for _ in range(120):
            focus = self.app.focusWidget()
            if focus is target:
                self.assertTrue(target.isVisible())
                self.assertTrue(target.isEnabled())
                self.assertFalse(target.visibleRegion().isEmpty())
                print(
                    f"Tab route: {len(sequence)} → {target.accessibleName() or target.objectName() or type(target).__name__}",
                    flush=True,
                )
                return len(sequence)
            sequence.append(type(focus).__name__)
            self.key(Qt.Key.Key_Backtab if backwards else Qt.Key.Key_Tab)
        self.fail(
            f"Tab cannot reach {target.accessibleName() or target.objectName()}: {sequence}"
        )

    def visible_button_focus(self, button):
        self.assertIs(self.app.focusWidget(), button)
        focused = button.grab().toImage()
        self.key(Qt.Key.Key_Backtab)
        self.assertIsNot(self.app.focusWidget(), button)
        unfocused = button.grab().toImage()
        self.key(Qt.Key.Key_Tab)
        self.assertIs(self.app.focusWidget(), button)
        self.assertNotEqual(
            focused, unfocused, "Primary action has no visible focus indicator"
        )

    def modal(self, key, modifiers, interact):
        errors = []
        handled = []

        def handle():
            handled.append(True)
            try:
                interact(self.app.activeModalWidget())
            except Exception as error:
                errors.append(error)
                dialog = self.app.activeModalWidget()
                if dialog is not None:
                    dialog.reject()

        timer = QTimer()
        timer.setSingleShot(True)

        def timeout():
            errors.append(AssertionError("Modal keyboard recipe timed out"))
            dialog = self.app.activeModalWidget()
            if dialog is not None:
                dialog.reject()

        timer.timeout.connect(timeout)
        timer.start(30000)
        interaction = QTimer()
        interaction.setSingleShot(True)
        interaction.timeout.connect(handle)
        interaction.start(60)
        try:
            self.key(key, modifiers)
        finally:
            interaction.stop()
            timer.stop()
        if errors:
            raise errors[0]
        self.assertTrue(handled, "Shortcut did not open its modal dialog")

    def command(self, text):
        def choose(dialog):
            from vinkulum_studio.commands import CommandPalette

            self.assertIsInstance(dialog, CommandPalette)
            self.assertIs(self.app.focusWidget(), dialog.query)
            QTest.keyClicks(dialog.query, text)
            self.assertEqual(dialog.results.count(), 1)
            self.key(Qt.Key.Key_Return)

        self.modal(Qt.Key.Key_K, Qt.KeyboardModifier.ControlModifier, choose)

    def test_create_edit_undo_redo_save_at_normal_and_narrow_inspector(self):
        from vinkulum_studio.editor import EditorWindow
        from vinkulum_studio.cad_dialog import CadDialog
        from vinkulum_studio.document import load_project

        with tempfile.TemporaryDirectory() as tmp, patch("sys.excepthook") as unhandled:
            for width in (420, 280):
                with self.subTest(inspector=width):
                    window = EditorWindow()
                    window._discard_allowed = lambda: True
                    self.addCleanup(window.close)
                    window.new_project()
                    window.resize(1280, 860)
                    window.show()
                    window.activateWindow()
                    self.assertTrue(QTest.qWaitForWindowActive(window))
                    window.resizeDocks(
                        [window.docks["inspector"]], [width], Qt.Orientation.Horizontal
                    )
                    window.viewport.setFocus()  # initial desktop focus, before the recipe
                    initial_focus = self.app.focusWidget()

                    def cancel(dialog):
                        self.assertIsInstance(dialog, CadDialog)
                        self.key(Qt.Key.Key_Escape)

                    self.modal(Qt.Key.Key_G, Qt.KeyboardModifier.AltModifier, cancel)
                    QTest.qWait(150)
                    self.assertIs(self.app.focusWidget(), initial_focus)
                    self.assertFalse(window.project.bodies)

                    def create(dialog):
                        self.assertIsInstance(dialog, CadDialog)
                        self.assertIs(self.app.focusWidget(), dialog.operation)
                        self.key(Qt.Key.Key_Tab)
                        self.assertIs(self.app.focusWidget(), dialog.name)
                        self.key(Qt.Key.Key_Backtab)
                        self.assertIs(self.app.focusWidget(), dialog.operation)
                        self.key(Qt.Key.Key_Tab)
                        self.enter("Keyboard box 137")
                        for field, value in zip(dialog.dimensions, ("20", "30", "40")):
                            self.key(Qt.Key.Key_Tab)
                            self.assertIs(self.app.focusWidget(), field)
                            self.enter(value)
                        self.tab_to(dialog.apply_button)
                        self.visible_button_focus(dialog.apply_button)
                        self.key(Qt.Key.Key_Space)

                    self.modal(Qt.Key.Key_G, Qt.KeyboardModifier.AltModifier, create)
                    self.assertEqual(
                        len(window.project.bodies), 1, window.status.text()
                    )
                    created = window.project.bodies[0]
                    self.assertIsNotNone(created.cad)
                    self.assertEqual(created.name, "Keyboard box 137")
                    self.tab_to(window.fields["mass"])
                    self.enter("1.3267053771417465")
                    self.tab_to(window.fields["position"].components[0])
                    self.enter("1.234567890123456")
                    self.tab_to(window.apply_button)
                    self.visible_button_focus(window.apply_button)
                    self.key(Qt.Key.Key_Space)
                    edited = window.project.bodies[0]
                    self.assertEqual(edited.mass, 1.3267053771417465)
                    self.assertEqual(edited.position[0], 1.234567890123456)
                    self.assertEqual(edited.orientation, created.orientation)
                    for resized in (320, width):
                        window.resizeDocks(
                            [window.docks["inspector"]],
                            [resized],
                            Qt.Orientation.Horizontal,
                        )
                        QTest.qWait(30)
                        self.assertEqual(window.project.bodies[0], edited)
                        self.assertFalse(window.dirty_fields)
                    self.tab_to(window.fields["name"])
                    self.enter("Renamed box 137")
                    self.tab_to(window.apply_button)
                    self.visible_button_focus(window.apply_button)
                    self.key(Qt.Key.Key_Space)
                    renamed = window.project.bodies[0]
                    self.assertEqual(renamed, replace(edited, name="Renamed box 137"))
                    self.command("Undo")
                    self.assertEqual(window.project.bodies[0], edited)
                    self.command("Redo")
                    self.assertEqual(window.project.bodies[0], renamed)
                    path = Path(tmp) / f"keyboard-{width}.json"

                    def save(dialog):
                        self.assertIsInstance(dialog, QFileDialog)
                        self.assertIsInstance(self.app.focusWidget(), QLineEdit)
                        self.enter(str(path))
                        self.key(Qt.Key.Key_Return)

                    self.modal(Qt.Key.Key_S, Qt.KeyboardModifier.ControlModifier, save)
                    self.assertEqual(load_project(path), window.project)
                    self.command("Fit selection")
                    print(
                        f"keyboard CAD: DPI={window.devicePixelRatioF()} window={window.width()}x{window.height()} inspector={window.docks['inspector'].width()} requested={width} saved={path.name}",
                        flush=True,
                    )
                    out = os.environ.get("VINKULUM_KEYBOARD_SCREENSHOTS")
                    if out:
                        Path(out).mkdir(parents=True, exist_ok=True)
                        self.tab_to(window.fields["position"].components[0])
                        window.viewport.render()
                        QTest.qWait(80)
                        image = window.screen().grabWindow(window.winId())
                        self.assertFalse(image.isNull())
                        self.assertTrue(
                            image.save(str(Path(out) / f"keyboard-{width}.png"))
                        )
                    window.close()
            QTest.qWait(50)
            unhandled.assert_not_called()


if __name__ == "__main__":
    unittest.main()
