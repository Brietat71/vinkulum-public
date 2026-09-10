"""Actual sketch picking, decimal edits, conflicts, drawing, dragging and undo."""

import os
import unittest
from dataclasses import asdict, replace
from unittest.mock import patch

from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio.sketch import profile, solve
from vinkulum_studio.sketch_dialog import SketchDialog


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires a Qt desktop"
)
class SketchWorkspace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.guard = patch("sys.excepthook")
        self.errors = self.guard.start()
        self.addCleanup(self.guard.stop)
        self.addCleanup(self.errors.assert_not_called)
        self.window = SketchDialog()
        self.addCleanup(self.window.reject)
        self.window.show()
        self.window.activateWindow()
        QTest.qWait(40)

    def text(self, field, value):
        field.setFocus()
        QTest.keyClick(field, Qt.Key.Key_A, Qt.KeyboardModifier.ControlModifier)
        QTest.keyClicks(field, value)

    def choose(self, kind, key):
        self.assertTrue(self.window.select(((kind, key),)))
        QTest.qWait(10)

    def add(self, kind):
        w = self.window
        w.constraint_kind.setCurrentIndex(w.constraint_kind.findData(kind))
        QTest.mouseClick(w.add_constraint, Qt.MouseButton.LeftButton)

    def test_dimension_pick_edit_conflict_and_keyboard_undo_preserve_source(self):
        w = self.window
        original = w.profile
        before = asdict(original)
        width = next(c for c in original.constraints if c.name == "Width")
        rectangle = next(
            rect for rect, key in w.canvas._dimension_hits if key == width.id
        )
        QTest.mouseDClick(
            w.canvas, Qt.MouseButton.LeftButton, pos=rectangle.center().toPoint()
        )
        QTest.qWait(20)  # The newly created property field must become visible.
        self.assertEqual(w.selection, (("constraint", width.id),))
        self.assertTrue(w.value_fields[0].hasFocus())
        self.text(w.value_fields[0], "100")
        QTest.mouseClick(w.update_button, Qt.MouseButton.LeftButton)
        self.assertEqual(w.profile.points[1].xy_mm, (100.0, 0.0))
        self.assertTrue(w.accept_button.isEnabled())
        self.assertEqual(asdict(original), before)
        w.canvas.setFocus()
        QTest.keyClick(w.canvas, Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
        self.assertEqual(w.profile, original)
        QTest.keyClick(
            w.canvas,
            Qt.Key.Key_Z,
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier,
        )
        self.assertEqual(w.profile.points[1].xy_mm, (100.0, 0.0))
        self.choose("edge", w.profile.edges[0].id)
        self.add("distance_x")
        duplicate = w.selection[0][1]
        self.assertFalse(solve(w.profile).conflicts)
        self.assertEqual(solve(w.profile).redundant_equations, 1)
        self.text(w.value_fields[0], "120")
        QTest.mouseClick(w.update_button, Qt.MouseButton.LeftButton)
        self.assertIn(duplicate, solve(w.profile).conflicts)
        self.assertFalse(w.accept_button.isEnabled())
        self.assertIn("Conflicting", w.summary.text())
        QTest.mouseClick(w.remove_button, Qt.MouseButton.LeftButton)
        self.assertFalse(solve(w.profile).conflicts)
        self.assertEqual(w.profile.points[1].xy_mm, (100.0, 0.0))
        self.assertTrue(w.accept_button.isEnabled())

    def test_invalid_pending_dimension_keeps_text_and_selection_until_corrected(self):
        w = self.window
        width = next(c for c in w.profile.constraints if c.name == "Width")
        height = next(c for c in w.profile.constraints if c.name == "Height")
        self.choose("constraint", width.id)
        before = w.profile
        self.text(w.value_fields[0], "nan")
        item = w.items[("constraint", height.id)]
        QTest.mouseClick(
            w.tree.viewport(),
            Qt.MouseButton.LeftButton,
            pos=w.tree.visualItemRect(item).center(),
        )
        self.assertEqual(w.selection, (("constraint", width.id),))
        self.assertEqual(w.value_fields[0].text(), "nan")
        self.assertEqual(w.profile, before)
        self.assertIn("decimal", w.status.text())
        QTest.mouseClick(w.revert_button, Qt.MouseButton.LeftButton)
        self.assertEqual(w.value_fields[0].text(), "80")
        self.choose("constraint", height.id)

    def test_draw_close_constrain_drag_and_revert_one_transaction(self):
        w = self.window
        w.empty_action.trigger()
        w.canvas.setFocus()
        QTest.keyClick(w.canvas, Qt.Key.Key_L)
        for point in ((0, 0), (60, 0), (60, 40), (0, 40), (0, 0)):
            QTest.mouseClick(
                w.canvas,
                Qt.MouseButton.LeftButton,
                pos=w.canvas.screen(point).toPoint(),
            )
            QTest.qWait(10)
        self.assertEqual((len(w.profile.points), len(w.profile.edges)), (4, 4))
        self.assertTrue(w.accept_button.isEnabled())
        self.assertEqual(w.canvas.tool, "select")
        before = {p.id for p in w.profile.points}
        for i, edge in enumerate(w.profile.edges):
            if i == 0:
                QTest.mouseClick(
                    w.canvas,
                    Qt.MouseButton.LeftButton,
                    pos=w.canvas.screen((30, 0)).toPoint(),
                )
                self.assertEqual(w.selection, (("edge", edge.id),))
            else:
                self.choose("edge", edge.id)
            self.add("horizontal" if i % 2 == 0 else "vertical")
        self.choose("point", w.profile.points[0].id)
        self.add("fixed")
        self.assertEqual(solve(w.profile).degrees_of_freedom, 2)
        start, end = (
            w.canvas.screen((60, 40)).toPoint(),
            w.canvas.screen((80, 50)).toPoint(),
        )
        index = w.undo.index()
        saved = w.profile
        QTest.mousePress(w.canvas, Qt.MouseButton.LeftButton, pos=start)
        QTest.mouseMove(w.canvas, end, 30)
        QTest.mouseRelease(w.canvas, Qt.MouseButton.LeftButton, pos=end)
        self.assertEqual(w.undo.index(), index + 1)
        self.assertEqual(
            [p.xy_mm for p in w.profile.points], [(0, 0), (80, 0), (80, 50), (0, 50)]
        )
        self.assertEqual({p.id for p in w.profile.points}, before)
        self.assertEqual(w.profile.constraints, saved.constraints)
        w.canvas.setFocus()
        QTest.keyClick(w.canvas, Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
        self.assertEqual(w.profile, saved)

    def drag(self, start, end):
        canvas = self.window.canvas
        start, end = canvas.screen(start).toPoint(), canvas.screen(end).toPoint()
        QTest.mousePress(canvas, Qt.MouseButton.LeftButton, pos=start)
        QTest.mouseMove(canvas, end, 30)
        QTest.mouseRelease(canvas, Qt.MouseButton.LeftButton, pos=end)

    def test_drag_preserves_invalid_name_then_commits_corrected_fields_first(self):
        w = self.window
        w.preset(profile(((0, 0), (60, 0), (60, 40), (0, 40))))
        first = w.profile.points[0]
        self.choose("point", first.id)
        self.text(w.name_field, " ")
        before, selection, index = w.profile, w.selection, w.undo.index()
        self.drag((60, 40), (80, 50))
        self.assertEqual(w.profile, before)
        self.assertEqual(w.selection, selection)
        self.assertEqual(w.name_field.text(), " ")
        self.assertTrue(w._pending_fields)
        self.assertEqual(w.undo.index(), index)
        self.assertIsNone(w.canvas.preview)
        self.assertIsNone(w.canvas._drag)
        self.assertIn("names", w.status.text())

        self.text(w.name_field, "Origin")
        self.drag((60, 40), (80, 50))
        self.assertEqual(w.profile.points[0], replace(first, name="Origin"))
        self.assertEqual(w.profile.points[2].xy_mm, (80, 50))
        self.assertEqual(w.undo.index(), index + 2)
        w.undo.undo()  # Undo the gesture, retaining the preceding field edit.
        self.assertEqual(w.profile.points[2], before.points[2])
        self.assertEqual(w.profile.points[0].name, "Origin")
        w.undo.undo()
        self.assertEqual(w.profile, before)

    def test_drag_preserves_invalid_dimension_until_fields_are_reverted(self):
        w = self.window
        w.preset(profile(((0, 0), (60, 0), (60, 40), (0, 40))))
        self.choose("edge", w.profile.edges[0].id)
        self.add("distance_x")
        before, selection, index = w.profile, w.selection, w.undo.index()
        self.text(w.value_fields[0], "nan")
        self.drag((60, 40), (80, 50))
        self.assertEqual(w.profile, before)
        self.assertEqual(w.selection, selection)
        self.assertEqual(w.value_fields[0].text(), "nan")
        self.assertTrue(w._pending_fields)
        self.assertEqual(w.undo.index(), index)
        self.assertIn("decimal", w.status.text())
        QTest.mouseClick(w.revert_button, Qt.MouseButton.LeftButton)
        self.drag((60, 40), (80, 50))
        self.assertEqual(w.profile.points[2].xy_mm, (80, 50))
        self.assertEqual(w.profile.constraints, before.constraints)
        self.assertEqual(w.undo.index(), index + 1)

    def test_drag_of_current_point_also_respects_invalid_coordinate(self):
        w = self.window
        w.preset(profile(((0, 0), (60, 0), (60, 40), (0, 40))))
        point = w.profile.points[2]
        self.choose("point", point.id)
        self.text(w.value_fields[0], "nan")
        before, selection, index = w.profile, w.selection, w.undo.index()
        self.drag(point.xy_mm, (80, 50))
        self.assertEqual(w.profile, before)
        self.assertEqual(w.selection, selection)
        self.assertEqual(w.value_fields[0].text(), "nan")
        self.assertTrue(w._pending_fields)
        self.assertEqual(w.undo.index(), index)

    def test_deletion_requires_explicit_dependent_entities_and_cancel_keeps_input(self):
        w = self.window
        original = w.profile
        self.choose("point", original.points[0].id)
        QTest.mouseClick(w.remove_button, Qt.MouseButton.LeftButton)
        self.assertEqual(w.profile, original)
        self.assertIn("dependent", w.status.text())
        keys = tuple(w.items)
        self.assertTrue(w.select(keys))
        QTest.mouseClick(w.remove_button, Qt.MouseButton.LeftButton)
        self.assertFalse(w.profile.points or w.profile.edges or w.profile.constraints)
        self.assertFalse(w.accept_button.isEnabled())
        w.undo.undo()
        self.assertEqual(w.profile, original)
        w.reject()
        self.assertEqual(w.result(), SketchDialog.DialogCode.Rejected)

    def test_out_of_range_imported_dimensions_remain_inspectable_and_repairable(self):
        original = self.window.profile
        self.window.reject()
        origin = next(c for c in original.constraints if c.name == "Origin")
        imported = replace(
            original,
            constraints=tuple(
                replace(c, values_mm=("1000000", "0")) if c == origin else c
                for c in original.constraints
            ),
        )
        self.window = w = SketchDialog(imported)
        self.addCleanup(w.reject)
        w.show()
        w.activateWindow()
        QTest.qWait(40)
        self.assertFalse(w.accept_button.isEnabled())
        self.assertIn("unavailable", w.summary.text())
        point = w.profile.points[0]
        QTest.mouseClick(
            w.canvas,
            Qt.MouseButton.LeftButton,
            pos=w.canvas.screen(point.xy_mm).toPoint(),
        )
        self.assertEqual(w.selection, (("point", point.id),))
        self.assertTrue(all(not field.isEnabled() for field in w.value_fields))
        self.assertIn("exceed", w.status.text())
        self.choose("constraint", origin.id)
        self.text(w.value_fields[0], "0")
        QTest.mouseClick(w.update_button, Qt.MouseButton.LeftButton)
        self.assertTrue(w.accept_button.isEnabled())
        self.assertEqual(w.profile, original)

    def test_imported_profile_displays_solved_dimensions_without_mutating_input(self):
        original = self.window.profile
        self.window.reject()
        imported = replace(
            original,
            constraints=tuple(
                replace(c, values_mm=("100",)) if c.name == "Width" else c
                for c in original.constraints
            ),
        )
        before = asdict(imported)
        self.window = w = SketchDialog(imported)
        self.addCleanup(w.reject)
        w.show()
        w.activateWindow()
        QTest.qWait(40)
        point = w.profile.points[1]
        self.assertEqual(point.xy_mm, (100, 0))
        self.assertEqual(w.canvas.sketch, w.profile)
        self.choose("point", point.id)
        self.assertEqual(w.value_fields[0].value(), 100)
        self.assertEqual(asdict(imported), before)
        self.assertEqual(imported.points[1].xy_mm, (80, 0))
        self.assertTrue(w.accept_button.isEnabled())


if __name__ == "__main__":
    unittest.main()
