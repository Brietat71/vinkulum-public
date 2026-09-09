"""Qt/VTK viewport. Rendering and mouse previews never mutate a Project."""

import math

import numpy as np
import vtkmodules.vtkInteractionStyle  # Registers trackball interaction.
import vtkmodules.vtkRenderingOpenGL2
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QVBoxLayout, QWidget
from vtkmodules.vtkCommonCore import vtkPoints, vtkUnsignedCharArray
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkFiltersCore import vtkTubeFilter, vtkPolyDataNormals
from vtkmodules.vtkFiltersGeneral import vtkTransformFilter
from vtkmodules.vtkFiltersSources import (
    vtkArrowSource,
    vtkCubeSource,
    vtkCylinderSource,
    vtkLineSource,
    vtkRegularPolygonSource,
    vtkSphereSource,
)
from vtkmodules.vtkInteractionWidgets import (
    vtkBoxRepresentation,
    vtkBoxWidget2,
    vtkCameraOrientationWidget,
)
from vtkmodules.vtkIOImage import vtkPNGWriter
from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkCellPicker,
    vtkLightKit,
    vtkPolyDataMapper,
    vtkRenderer,
    vtkWindowToImageFilter,
)

from .render_interactor import RenderInteractor


def vtk_matrix(position, orientation):
    matrix = vtkMatrix4x4()
    R = np.array(orientation).reshape(3, 3)
    for i in range(3):
        for j in range(3):
            matrix.SetElement(i, j, R[i, j])
        matrix.SetElement(i, 3, position[i])
    return matrix


def numpy_pose(matrix):
    return (
        tuple(matrix.GetElement(i, 3) for i in range(3)),
        tuple(matrix.GetElement(i, j) for i in range(3) for j in range(3)),
    )


class Viewport(QWidget):
    selected = Signal(object)
    pose_committed = Signal(str, object, object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(400, 300)
        self.setAccessibleName("Scène mécanique 3D")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.view = RenderInteractor(self)
        layout.addWidget(self.view)
        self.renderer = vtkRenderer()
        self.renderer.SetBackground(0.075, 0.084, 0.10)
        self.renderer.SetBackground2(0.16, 0.18, 0.21)
        self.renderer.GradientBackgroundOn()
        self.view.GetRenderWindow().AddRenderer(self.renderer)
        self.view.GetRenderWindow().SetMultiSamples(4)
        self.view.SetInteractorStyle(
            vtkmodules.vtkInteractionStyle.vtkInteractorStyleTrackballCamera()
        )
        self.view.Initialize()
        self.lights = vtkLightKit()
        self.lights.SetKeyLightIntensity(0.7)
        self.lights.SetKeyToFillRatio(2.0)
        self.lights.SetKeyToBackRatio(3.0)
        self.lights.AddLightsToRenderer(self.renderer)
        self.orientation_widget = vtkCameraOrientationWidget()
        self.orientation_widget.SetParentRenderer(self.renderer)
        self.orientation_widget.AnimateOff()
        compass = self.orientation_widget.GetRepresentation()
        compass.SetSize(96, 96)
        compass.SetPadding(16, 16)
        compass.SetXAxisColor(0.82, 0.39, 0.39)
        compass.SetYAxisColor(0.46, 0.70, 0.52)
        compass.SetZAxisColor(0.43, 0.62, 0.88)
        compass.SetNormalizedHandleDia(0.30)
        for axis in ("X", "Y", "Z"):
            for side in ("Plus", "Minus"):
                label = getattr(compass, f"Get{axis}{side}LabelProperty")()
                label.SetFontFamilyToArial()
                label.ItalicOff()
                label.BoldOff()
                label.SetColor(0.94, 0.96, 0.99)
        self.orientation_widget.On()
        self._actors = {}
        self._colors = {}
        self._bounds = {}
        self._body_ids = set()
        self._joint_parts = []
        self._construction = set()
        self._grid_actor = None
        self.grid_visible = True
        self._load_glyphs = []
        self._time = 0.0
        self._local_axes = None
        self._selected = None
        self._editable = True
        self.transform_mode = 3
        self.hidden_objects = set()
        self.dragging = False
        self._closed = False
        self._project = None
        self._poses = {}
        self._sources = []
        self.box = vtkBoxWidget2()
        self.box.SetInteractor(self.view.GetRenderWindow().GetInteractor())
        self.box.ScalingEnabledOff()
        self.box.MoveFacesEnabledOff()
        self.box.TranslationEnabledOn()
        self.box.RotationEnabledOn()
        self.box_rep = vtkBoxRepresentation()
        self.box_rep.SetPlaceFactor(1.15)
        self.box.SetRepresentation(self.box_rep)
        self.box.AddObserver("InteractionEvent", self._preview_pose)
        self.box.AddObserver(
            "StartInteractionEvent", lambda *_: setattr(self, "dragging", True)
        )
        self.box.AddObserver("EndInteractionEvent", self._commit_pose)
        self.view.AddObserver("LeftButtonPressEvent", self._pick, 0.1)

    def _add(self, source, color, identifier=None, width=2.0):
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*color)
        actor.GetProperty().SetLineWidth(width)
        actor.GetProperty().SetSpecular(0.2)
        self.renderer.AddActor(actor)
        if identifier is not None:
            self._actors.setdefault(identifier, []).append(actor)
            self._colors.setdefault(identifier, []).append(color)
        else:
            actor.PickableOff()
        self._sources.append(source)
        return actor

    def set_project(self, project, *, fit=False):
        if project == self._project:
            if fit:
                self.camera("iso")
            return
        selection = self._selected
        self.box.Off()
        self.renderer.RemoveAllViewProps()
        self._actors.clear()
        self._colors.clear()
        self._bounds.clear()
        self._sources.clear()
        self._joint_parts.clear()
        self._construction.clear()
        self._load_glyphs.clear()
        self._time = 0.0
        self._project = project
        self._body_ids = {b.id for b in project.bodies}
        self._poses = {b.id: (b.position, b.orientation) for b in project.bodies}
        for body in project.bodies:
            if body.cad is not None:
                points, cells = vtkPoints(), vtkCellArray()
                for point in body.cad.vertices_m:
                    points.InsertNextPoint(*point)
                for triangle in body.cad.triangles:
                    cells.InsertNextCell(3)
                    for index in triangle:
                        cells.InsertCellPoint(index)
                mesh = vtkPolyData()
                mesh.SetPoints(points)
                mesh.SetPolys(cells)
                normals = vtkPolyDataNormals()
                normals.SetInputData(mesh)
                normals.SetFeatureAngle(35)
                normals.SplittingOn()
                normals.ConsistencyOn()
                normals.Update()
                source = normals
            elif body.shape == "box":
                source = vtkCubeSource()
                source.SetXLength(body.dimensions[0])
                source.SetYLength(body.dimensions[1])
                source.SetZLength(body.dimensions[2])
            elif body.shape == "sphere":
                source = vtkSphereSource()
                source.SetRadius(body.dimensions[0])
                source.SetThetaResolution(64)
                source.SetPhiResolution(48)
            else:
                cylinder = vtkCylinderSource()
                cylinder.SetRadius(body.dimensions[0])
                cylinder.SetHeight(body.dimensions[1])
                cylinder.SetResolution(64)
                transform = vtkTransform()
                transform.RotateX(90)
                source = vtkTransformFilter()
                source.SetInputConnection(cylinder.GetOutputPort())
                source.SetTransform(transform)
            source.Update()
            self._bounds[body.id] = source.GetOutput().GetBounds()
            color = ((0.56, 0.67, 0.81), (0.69, 0.75, 0.80), (0.40, 0.61, 0.64))[
                (len(self._bounds) - 1) % 3
            ]
            actor = self._add(source, color, body.id)
            actor.GetProperty().SetAmbient(0.15)
            actor.GetProperty().SetDiffuse(0.75)
            actor.GetProperty().SetSpecular(0.18)
            actor.GetProperty().SetSpecularPower(55)
            actor.GetProperty().SetEdgeColor(0.22, 0.29, 0.37)
            actor.GetProperty().SetEdgeVisibility(body.shape == "box")
            actor.SetUserMatrix(vtk_matrix(body.position, body.orientation))
        extent = max(
            (np.linalg.norm(b.position) + max(b.dimensions) for b in project.bodies),
            default=1.0,
        )
        self._extent = max(0.1, float(extent))
        self._draw_grid()
        self._local_axes = vtkAxesActor()
        self._local_axes.AxisLabelsOff()
        self._local_axes.PickableOff()
        self._local_axes.VisibilityOff()
        self.renderer.AddActor(self._local_axes)
        invalid = {d.object_id for d in project.diagnostics()}
        for joint in project.joints:
            if any(
                ref is not None and ref not in self._body_ids
                for ref in (joint.a, joint.b)
            ):
                continue
            color = (0.95, 0.35, 0.3) if joint.id in invalid else (0.84, 0.64, 0.32)
            # Each attachment is drawn separately, so mismatched anchors remain visible.
            for side, (ref, point, frame) in enumerate(
                (
                    (joint.a, joint.pa, joint.ra),
                    (joint.b, joint.pb, joint.rb),
                )
            ):
                if joint.kind == "pivot":
                    # Two concentric symbols at the actual attachment positions.
                    # Their radii differ to avoid coincident surfaces; neither is
                    # physical geometry or a modification of the joint anchors.
                    ring = vtkRegularPolygonSource()
                    ring.SetNumberOfSides(48)
                    ring.SetRadius(self._extent * (0.022 if side == 0 else 0.014))
                    ring.GeneratePolygonOff()
                    tube = vtkTubeFilter()
                    tube.SetInputConnection(ring.GetOutputPort())
                    tube.SetRadius(self._extent * 0.0025)
                    tube.SetNumberOfSides(10)
                    actor = self._add(tube, color, joint.id)
                    self._joint_parts.append((actor, ref, point, frame, True))
                else:
                    sphere = vtkSphereSource()
                    sphere.SetRadius(self._extent * (0.014 if side == 0 else 0.010))
                    sphere.SetThetaResolution(24)
                    sphere.SetPhiResolution(16)
                    actor = self._add(sphere, color, joint.id)
                    self._joint_parts.append((actor, ref, point, frame, False))
                if joint.kind in {"pivot", "glissiere"}:
                    axis = vtkLineSource()
                    axis.SetPoint1(0.0, 0.0, 0.0)
                    axis.SetPoint2(0.0, 0.0, (1 if side else -1) * 0.05 * self._extent)
                    tube = vtkTubeFilter()
                    tube.SetInputConnection(axis.GetOutputPort())
                    tube.SetRadius(self._extent * 0.0025)
                    tube.SetNumberOfSides(10)
                    tube.CappingOn()
                    actor = self._add(tube, color, joint.id)
                    self._joint_parts.append((actor, ref, point, frame, True))
            # Link body centre to anchor as a visual attachment, not a mass-bearing rod.
            for ref, point in ((joint.a, joint.pa), (joint.b, joint.pb)):
                if ref is not None:
                    line = vtkLineSource()
                    line.SetPoint1(0.0, 0.0, 0.0)
                    line.SetPoint2(*point)
                    actor = self._add(line, color, joint.id)
                    self._construction.add(actor)
                    self._joint_parts.append(
                        (
                            actor,
                            ref,
                            (0.0, 0.0, 0.0),
                            (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
                            True,
                        )
                    )
        for load in project.loads:
            if load.body not in self._body_ids:
                continue
            sphere = vtkSphereSource()
            sphere.SetRadius(self._extent * 0.018)
            actor = self._add(sphere, (0.86, 0.46, 0.9), load.id)
            self._joint_parts.append(
                (
                    actor,
                    load.body,
                    load.point,
                    (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
                    False,
                )
            )
            for field, color in (
                ("force", (0.98, 0.50, 0.25)),
                ("moment", (0.80, 0.40, 0.95)),
            ):
                arrow = vtkArrowSource()
                arrow.SetTipResolution(20)
                arrow.SetShaftResolution(16)
                glyph = self._add(arrow, color, load.id)
                self._load_glyphs.append((glyph, load, field))
        self.hidden_objects.intersection_update(self._actors)
        self._update_attachments()
        self._apply_visibility()
        self.select(selection if selection in self._actors else None)
        if fit:
            self.camera("iso")
        self.render()

    def _draw_grid(self):
        # One segmented actor: a metric XY grid fading towards its boundary.
        # It stays on z=0 and is never presented as a contact/support surface.
        base = 10 ** math.floor(math.log10(self._extent / 5))
        self.grid_step = (
            min((1, 2, 5, 10), key=lambda n: abs(n * base - self._extent / 5)) * base
        )
        radius = 8 * self.grid_step
        points, lines = vtkPoints(), vtkCellArray()
        colors = vtkUnsignedCharArray()
        colors.SetNumberOfComponents(4)
        for i in range(-8, 9):
            for transpose in (False, True):
                for segment in range(32):
                    a = [i * self.grid_step, -radius + segment * radius / 16, 0.0]
                    b = [a[0], a[1] + radius / 16, 0.0]
                    distance = math.hypot(a[0], (a[1] + b[1]) / 2) / radius
                    alpha = round((65 if i == 0 else 36) * max(0, 1 - distance) ** 1.5)
                    if not alpha:
                        continue
                    if transpose:
                        a[:2], b[:2] = a[1::-1], b[1::-1]
                    lines.InsertNextCell(2)
                    lines.InsertCellPoint(points.InsertNextPoint(*a))
                    lines.InsertCellPoint(points.InsertNextPoint(*b))
                    colors.InsertNextTuple4(125, 143, 163, alpha)
        data = vtkPolyData()
        data.SetPoints(points)
        data.SetLines(lines)
        data.GetCellData().SetScalars(colors)
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(data)
        mapper.SetColorModeToDirectScalars()
        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().LightingOff()
        actor.PickableOff()
        actor.SetVisibility(self.grid_visible)
        self.renderer.AddActor(actor)
        self._grid_actor = actor

    def set_grid_visible(self, visible):
        self.grid_visible = bool(visible)
        if self._grid_actor is not None:
            self._grid_actor.SetVisibility(visible)
        self.render()

    def set_parallel_projection(self, parallel):
        camera = self.renderer.GetActiveCamera()
        if bool(camera.GetParallelProjection()) == bool(parallel):
            return
        # Preserve the apparent scale when switching projection at the focal plane.
        tangent = math.tan(math.radians(camera.GetViewAngle()) / 2)
        if parallel:
            camera.SetParallelScale(camera.GetDistance() * tangent)
        else:
            direction = np.array(camera.GetDirectionOfProjection())
            position = (
                np.array(camera.GetFocalPoint())
                - direction * camera.GetParallelScale() / tangent
            )
            camera.SetPosition(*position)
        camera.SetParallelProjection(parallel)
        self.renderer.ResetCameraClippingRange()
        self.render()

    def _update_attachments(self):
        for actor, reference, point, frame, oriented in self._joint_parts:
            if reference is None:
                position, R = np.zeros(3), np.eye(3)
            else:
                position, rotation = self._poses[reference]
                position, R = np.array(position), np.array(rotation).reshape(3, 3)
            orientation = R @ np.array(frame).reshape(3, 3) if oriented else np.eye(3)
            actor.SetUserMatrix(vtk_matrix(position + R @ point, orientation.flat))
        for actor, load, field in self._load_glyphs:
            values = np.array([law.value(self._time) for law in getattr(load, field)])
            norm = np.linalg.norm(values)
            if not np.isfinite(norm) or norm < 1e-14:
                actor.VisibilityOff()
                continue
            actor.VisibilityOn()
            x = values / norm
            seed = np.eye(3)[int(np.argmin(abs(x)))]
            z = np.cross(x, seed)
            z /= np.linalg.norm(z)
            y = np.cross(z, x)
            position, R = self._poses[load.body]
            R = np.array(R).reshape(3, 3)
            world = np.array(position) + R @ load.point
            frame = np.column_stack((x, y, z)) * (self._extent * 0.22)
            actor.SetUserMatrix(vtk_matrix(world, frame.flat))
        if self._local_axes is not None and self._selected in self._body_ids:
            self._local_axes.SetUserMatrix(vtk_matrix(*self._poses[self._selected]))

    def set_poses(self, poses, time=0.0):
        """Update result transforms without reconstructing actors or touching input data."""
        self._poses = dict(poses)
        self._time = float(time)
        for identifier, (position, orientation) in poses.items():
            if identifier in self._body_ids:
                self._actors[identifier][0].SetUserMatrix(
                    vtk_matrix(position, orientation)
                )
        self._update_attachments()
        self._apply_visibility()
        self.render()

    def set_editable(self, editable):
        self._editable = editable
        self.select(self._selected)

    def select(self, identifier):
        self.box.Off()
        self._selected = identifier
        if self._local_axes is not None:
            self._local_axes.SetVisibility(
                identifier in self._body_ids
                and identifier not in self.hidden_objects
                and self._editable
                and self.transform_mode
            )
            if identifier in self._body_ids:
                body = next(b for b in self._project.bodies if b.id == identifier)
                self._local_axes.SetTotalLength(*([max(body.dimensions) * 0.25] * 3))
                self._local_axes.SetUserMatrix(vtk_matrix(*self._poses[identifier]))
        for key, actors in self._actors.items():
            if key in self._body_ids:
                color = np.array(self._colors[key][0])
                if key == identifier:
                    color = 0.7 * color + 0.3 * np.array((0.6, 0.73, 1.0))
                actors[0].GetProperty().SetColor(*color)
                actors[0].GetProperty().SetEdgeColor(
                    *((0.6, 0.73, 1.0) if key == identifier else (0.22, 0.29, 0.37))
                )
                actors[0].GetProperty().SetLineWidth(2.0 if key == identifier else 1.0)
            else:
                for actor, color in zip(actors, self._colors[key]):
                    actor.GetProperty().SetColor(
                        *((1.0, 0.85, 0.3) if key == identifier else color)
                    )
        if (
            identifier in self._body_ids
            and self._editable
            and self.transform_mode
            and identifier not in self.hidden_objects
        ):
            self.box_rep.PlaceWidget(self._bounds[identifier])
            transform = vtkTransform()
            transform.SetMatrix(self._actors[identifier][0].GetMatrix())
            self.box_rep.SetTransform(transform)
            self.box.On()
        self._apply_visibility()
        self.render()

    def set_transform_mode(self, mode):
        self.transform_mode = mode
        self.box.SetTranslationEnabled(mode in (1, 3))
        self.box.SetRotationEnabled(mode in (2, 3))
        self.select(self._selected)

    def _apply_visibility(self):
        glyphs = {actor for actor, _, _ in self._load_glyphs}
        for identifier, actors in self._actors.items():
            for actor in actors:
                if identifier in self.hidden_objects:
                    actor.VisibilityOff()
                elif actor in self._construction:
                    actor.SetVisibility(identifier == self._selected)
                elif actor not in glyphs:
                    actor.VisibilityOn()

    def set_object_visible(self, identifier, visible):
        if visible:
            self.hidden_objects.discard(identifier)
        else:
            self.hidden_objects.add(identifier)
        self._update_attachments()
        self._apply_visibility()
        self.select(self._selected)

    def isolate(self, identifier):
        self.hidden_objects = set(self._actors) - {identifier}
        self._update_attachments()
        self._apply_visibility()
        self.select(self._selected)

    def show_all(self):
        self.hidden_objects.clear()
        self._update_attachments()
        self._apply_visibility()
        self.select(self._selected)

    def _pick(self, interactor, event):
        picker = vtkCellPicker()
        picker.SetTolerance(0.005)
        x, y = self.view.GetEventPosition()
        picker.Pick(x, y, 0, self.renderer)
        actor = picker.GetActor()
        for identifier, actors in self._actors.items():
            if actor in actors:
                if identifier != self._selected:
                    self.select(identifier)
                    self.selected.emit(identifier)
                return

    def _preview_pose(self, widget, event):
        if self._selected not in self._body_ids or not self._editable:
            return
        transform = vtkTransform()
        self.box_rep.GetTransform(transform)
        pose = numpy_pose(transform.GetMatrix())
        self._poses[self._selected] = pose
        self._actors[self._selected][0].SetUserMatrix(transform.GetMatrix())
        self._update_attachments()
        self._apply_visibility()
        self.render()

    def _commit_pose(self, widget, event):
        self.dragging = False
        if self._selected in self._body_ids and self._editable:
            transform = vtkTransform()
            self.box_rep.GetTransform(transform)
            position, orientation = numpy_pose(transform.GetMatrix())
            self.pose_committed.emit(self._selected, position, orientation)

    def camera(self, direction="iso", selection=False):
        camera = self.renderer.GetActiveCamera()
        vectors = {
            "iso": (1.0, -1.0, 0.8),
            "front": (0.0, -1.0, 0.0),
            "side": (1.0, 0.0, 0.0),
            "top": (0.0, 0.0, 1.0),
        }
        vector = vectors[direction]
        camera.SetFocalPoint(0.0, 0.0, 0.0)
        camera.SetPosition(*vector)
        camera.SetViewUp(*((0.0, 1.0, 0.0) if direction == "top" else (0.0, 0.0, 1.0)))
        camera.SetParallelProjection(direction != "iso")
        self.fit_scene(selection)

    def fit_scene(self, selection=False):
        if self._closed or self.dragging:
            return
        camera = self.renderer.GetActiveCamera()
        if selection and self._selected in self._actors:
            self.renderer.ResetCamera(self._actors[self._selected][0].GetBounds())
        elif self._body_ids:
            bounds = [
                self._actors[identifier][0].GetBounds() for identifier in self._body_ids
            ]
            scene_bounds = tuple(
                (min if index % 2 == 0 else max)(b[index] for b in bounds)
                for index in range(6)
            )
            self.renderer.ResetCamera(scene_bounds)
        else:
            self.renderer.ResetCamera()
        camera.Zoom(0.82)
        self.renderer.ResetCameraClippingRange()
        self.render()

    def render(self):
        if not self._closed and self.isVisible():
            self.view.GetRenderWindow().Render()

    def screenshot(self, path):
        self.view.GetRenderWindow().Render()
        capture = vtkWindowToImageFilter()
        capture.SetInput(self.view.GetRenderWindow())
        capture.ReadFrontBufferOff()
        capture.Update()
        writer = vtkPNGWriter()
        writer.SetFileName(str(path))
        writer.SetInputConnection(capture.GetOutputPort())
        writer.Write()

    def shutdown(self):
        if not self._closed:
            self.box.Off()
            self.orientation_widget.Off()
            self.view.Finalize()
            self._closed = True

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, "orientation_widget"):
            # VTK display coordinates are physical pixels, Qt sizes are logical.
            ratio = self.devicePixelRatioF()
            compass = self.orientation_widget.GetRepresentation()
            compass.SetSize(round(96 * ratio), round(96 * ratio))
            compass.SetPadding(round(16 * ratio), round(16 * ratio))

    def closeEvent(self, event):
        self.shutdown()
        event.accept()
