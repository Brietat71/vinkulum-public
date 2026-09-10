"""Run with VINKULUM_3D_TESTS=1 QT_QPA_PLATFORM=xcb under Xvfb or a real X server."""

import os
import unittest
from dataclasses import replace

import numpy as np
from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QPaintEvent
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio.document import Body, new_id, pendulum


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1",
    "VTK integration requires a graphics-enabled test run",
)
class ViewportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_fit_includes_anchors_respects_visibility_and_world_translation(self):
        from vinkulum_studio.viewport import Viewport

        view = Viewport()
        self.addCleanup(view.close)
        view.resize(900, 650)
        view.show()
        QTest.qWait(50)

        def visible(point):
            view.renderer.SetWorldPoint(*point, 1.0)
            view.renderer.WorldToDisplay()
            x, y, depth = view.renderer.GetDisplayPoint()
            width, height = view.view.GetRenderWindow().GetSize()
            self.assertTrue(
                5 < x < width - 5 and 5 < y < height - 5 and 0 < depth < 1,
                (point, x, y, depth),
            )

        project = pendulum()
        view.set_project(project, fit=True)
        visible(project.bodies[0].position)
        self.assertIsNone(project.joints[0].a)
        visible(project.joints[0].pa)
        extent = view._extent
        distance = view.renderer.GetActiveCamera().GetDistance()
        shift = np.array((1000.0, -2000.0, 3000.0))
        moved = replace(
            project,
            bodies=tuple(
                replace(body, position=tuple(np.array(body.position) + shift))
                for body in project.bodies
            ),
            joints=tuple(
                replace(
                    joint,
                    pa=tuple(np.array(joint.pa) + shift)
                    if joint.a is None
                    else joint.pa,
                    pb=tuple(np.array(joint.pb) + shift)
                    if joint.b is None
                    else joint.pb,
                )
                for joint in project.joints
            ),
        )
        view.set_project(moved, fit=True)
        self.assertAlmostEqual(view._extent / extent, 1.0, places=9)
        self.assertAlmostEqual(
            view.renderer.GetActiveCamera().GetDistance() / distance, 1.0, places=8
        )
        visible(moved.bodies[0].position)
        visible(moved.joints[0].pa)

        far = Body(new_id(), "Hidden outlier", position=(100.0, 0.0, 0.0))
        view.set_project(
            replace(project, joints=(), bodies=(*project.bodies, far)), fit=True
        )
        view.set_object_visible(far.id, False)
        view.fit_scene()
        np.testing.assert_allclose(
            view.renderer.GetActiveCamera().GetFocalPoint(),
            project.bodies[0].position,
            rtol=0,
            atol=1e-5,
        )

        joint = replace(project.joints[0], pa=(-2.0, 0.0, 0.0), pb=(2.0, 0.0, 0.0))
        view.set_project(replace(project, joints=(joint,)), fit=True)
        view.select(joint.id)
        view.fit_scene(selection=True)
        visible(joint.pa)
        body = project.bodies[0]
        visible(
            np.array(body.position)
            + np.array(body.orientation).reshape(3, 3) @ joint.pb
        )

    def test_scene_picking_manipulation_result_and_shutdown(self):
        from vinkulum_studio.viewport import Viewport, numpy_pose, vtk_matrix
        from vtkmodules.vtkCommonTransforms import vtkTransform

        project = pendulum()
        project = project.replace_object(
            Body(new_id(), "Box", position=(-0.4, 0.0, 0.0))
        )
        project = project.replace_object(
            Body(
                new_id(),
                "Cylinder",
                shape="cylinder",
                dimensions=(0.08, 0.4),
                position=(0.8, 0.0, 0.0),
            )
        )
        view = Viewport()
        try:
            view.resize(900, 650)
            view.show()
            view.set_project(project, fit=True)
            QTest.qWait(100)
            # Model Cocoa's feedback: an unsolicited paint must not cause yet
            # another VTK render; an explicit update must still repaint.
            view.view._coalesce_paints = True
            renders = []
            window = view.view.GetRenderWindow()
            observer = window.AddObserver("StartEvent", lambda *_: renders.append(1))
            view.view.update()
            QTest.qWait(30)
            self.assertTrue(renders)
            before = len(renders)
            for _ in range(20):
                self.app.sendEvent(view.view, QPaintEvent(view.view.rect()))
            self.assertEqual(len(renders), before)
            view.view.update()
            QTest.qWait(30)
            self.assertGreater(len(renders), before)
            window.RemoveObserver(observer)
            camera_before = view.renderer.GetActiveCamera().GetPosition()
            QTest.mousePress(view.view, Qt.MouseButton.LeftButton, pos=QPoint(40, 40))
            QTest.mouseMove(view.view, QPoint(110, 90), 20)
            QTest.mouseRelease(
                view.view, Qt.MouseButton.LeftButton, pos=QPoint(110, 90)
            )
            self.assertNotEqual(
                view.renderer.GetActiveCamera().GetPosition(), camera_before
            )
            view.camera("iso")
            body = project.bodies[1]
            selected = []
            view.selected.connect(selected.append)
            # Exercise the same picker called by the interactor mouse observer.
            view.renderer.SetWorldPoint(*body.position, 1.0)
            view.renderer.WorldToDisplay()
            x, y, _ = view.renderer.GetDisplayPoint()
            view.view.SetEventPosition(int(x), int(y))
            view._pick(None, None)
            self.assertEqual(selected[-1], body.id)
            self.assertTrue(view.box.GetEnabled())
            committed = []
            view.pose_committed.connect(lambda *args: committed.append(args))
            transform = vtkTransform()
            new_position = (-0.2, 0.1, 0.2)
            transform.SetMatrix(vtk_matrix(new_position, body.orientation))
            view.box_rep.SetTransform(transform)
            view.box.InvokeEvent("InteractionEvent")
            self.assertEqual(project.bodies[1].position, body.position)
            self.assertFalse(committed)
            view.box.InvokeEvent("EndInteractionEvent")
            self.assertEqual(committed[0][0], body.id)
            for actual, expected in zip(committed[0][1], new_position):
                self.assertAlmostEqual(actual, expected)
            actor = view._actors[body.id][0]
            view.set_editable(False)
            self.assertFalse(view.box.GetEnabled())
            poses = {b.id: (b.position, b.orientation) for b in project.bodies}
            poses[body.id] = ((0.3, 0.2, 0.1), body.orientation)
            view.set_poses(poses)
            self.assertIs(view._actors[body.id][0], actor)
            self.assertEqual(numpy_pose(actor.GetMatrix())[0], (0.3, 0.2, 0.1))
            for camera in ("front", "side", "top", "iso"):
                view.camera(camera)
            screenshot = os.environ.get("VINKULUM_3D_SCREENSHOT")
            if screenshot:
                view.screenshot(screenshot)
        finally:
            view.shutdown()
            view.close()


if __name__ == "__main__":
    unittest.main()
