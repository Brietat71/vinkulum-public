"""Captured CAD/mesh surface display with actual boundary-cell picking."""

import numpy as np
from PySide6.QtCore import QEvent, Qt, Signal
from vtkmodules.vtkCommonCore import vtkPoints, vtkUnsignedCharArray
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkPolyData,
    vtkQuadraticTriangle,
    vtkTriangle,
    vtkUnstructuredGrid,
)
from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkCellPicker,
    vtkDataSetMapper,
    vtkPolyDataMapper,
)

from .boundary_display import boundary_glyphs
from .boundary_view import BoundaryGlyphLayer
from .viewport import Viewport, vtk_matrix

SURFACE_COLORS = {
    "plain": (131, 167, 214),
    "selected": (171, 183, 255),
    "hover": (206, 214, 250),
    "support": (78, 186, 160),
    "pressure": (240, 174, 76),
    "total_force": (221, 113, 132),
    "multiple": (166, 130, 200),
}


class MeshViewport(Viewport):
    surface_picked = Signal(object, bool)
    surface_hovered = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAccessibleName("Captured solid and selectable finite-element faces")
        self._editable = False
        self.box.Off()
        self.selection_enabled = True
        self._press = None
        self._hover = None
        self._selected_faces = set()
        self._face_kinds = {}
        self._cell_faces = []
        self.solid = None
        self.mesh = vtkUnstructuredGrid()
        self.mapper = vtkDataSetMapper()
        self.mapper.SetInputData(self.mesh)
        self.mapper.SetColorModeToDirectScalars()
        self.mapper.SetScalarModeToUseCellData()
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().EdgeVisibilityOff()
        self.actor.GetProperty().SetEdgeColor(0.12, 0.16, 0.22)
        self.actor.GetProperty().SetAmbient(0.16)
        self.actor.GetProperty().SetDiffuse(0.40)
        self.actor.GetProperty().SetSpecular(0.08)
        self.actor.GetProperty().SetSpecularPower(30)
        self.renderer.AddActor(self.actor)
        self.edges_visible = True
        self.edge_mapper = vtkPolyDataMapper()
        self.edge_actor = vtkActor()
        self.edge_actor.SetMapper(self.edge_mapper)
        self.edge_actor.GetProperty().SetColor(0.16, 0.21, 0.28)
        self.edge_actor.GetProperty().SetLineWidth(1.0)
        self.edge_actor.PickableOff()
        self.renderer.AddActor(self.edge_actor)
        self.source_mapper = vtkPolyDataMapper()
        self.source_actor = vtkActor()
        self.source_actor.SetMapper(self.source_mapper)
        self.source_actor.GetProperty().SetColor(0.52, 0.67, 0.85)
        self.source_actor.GetProperty().SetSpecular(0.2)
        self.source_actor.PickableOff()
        self.renderer.AddActor(self.source_actor)
        self.picker = vtkCellPicker()
        self.picker.SetTolerance(0.002)
        self.picker.PickFromListOn()
        self.picker.AddPickList(self.actor)
        self.view.setMouseTracking(True)
        self.view.installEventFilter(self)
        self._frames = {}
        self.glyph_layer = BoundaryGlyphLayer(self, SURFACE_COLORS)

    def _pick(self, interactor, event):
        # Base Viewport picks bodies on press; faces are selected on a click's
        # release so camera orbiting does not change the boundary selection.
        pass

    def set_source(self, body):
        self.glyph_layer.set_glyphs(())
        self._frames = {}
        self._face_kinds = {}
        points, cells = vtkPoints(), vtkCellArray()
        points.SetDataTypeToDouble()
        for point in body.cad.vertices_m:
            points.InsertNextPoint(*point)
        for triangle in body.cad.triangles:
            cells.InsertNextCell(3)
            for node in triangle:
                cells.InsertCellPoint(node)
        data = vtkPolyData()
        data.SetPoints(points)
        data.SetPolys(cells)
        self.source_normals = vtkPolyDataNormals()
        self.source_normals.SetInputData(data)
        self.source_normals.SetFeatureAngle(35)
        self.source_mapper.SetInputConnection(self.source_normals.GetOutputPort())
        self.source_actor.SetUserMatrix(vtk_matrix(body.position, body.orientation))
        self.source_actor.VisibilityOn()
        self.actor.VisibilityOff()
        self.edge_actor.VisibilityOff()
        self.solid = None
        self._cell_faces = []
        self._hover = None
        self._selected_faces.clear()
        self.camera("iso")

    def set_solid(self, solid, *, frames=None, fit=True):
        self.glyph_layer.set_glyphs(())
        self._frames = frames if frames is not None else {}
        self._face_kinds = {}
        self.solid = solid
        self._hover = None
        self._cell_faces = []
        points = vtkPoints()
        points.SetDataTypeToDouble()
        for point in solid.mesh.nodes:
            points.InsertNextPoint(*point)
        self.mesh = vtkUnstructuredGrid()
        self.mesh.SetPoints(points)
        for surface in solid.surfaces:
            for nodes in surface.triangles:
                cell = vtkQuadraticTriangle() if len(nodes) == 6 else vtkTriangle()
                for j, node in enumerate(nodes):
                    cell.GetPointIds().SetId(j, node - 1)
                self.mesh.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
                self._cell_faces.append(surface.id)
        self.mapper.SetInputData(self.mesh)
        # Draw actual finite-element edges, not the internal triangles VTK
        # introduces when tessellating a quadratic boundary face for display.
        edge_points, lines, seen = vtkPoints(), vtkCellArray(), set()
        edge_points.SetDataTypeToDouble()
        coordinates = np.array(solid.mesh.nodes)
        for surface in solid.surfaces:
            for nodes in surface.triangles:
                for j, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
                    key = tuple(sorted((nodes[a], nodes[b])))
                    if key in seen:
                        continue
                    seen.add(key)
                    start, end = coordinates[np.array((nodes[a], nodes[b])) - 1]
                    if len(nodes) == 6:
                        mid = coordinates[nodes[j + 3] - 1]
                        curve = [
                            (1 - t) * (1 - 2 * t) * start
                            + t * (2 * t - 1) * end
                            + 4 * t * (1 - t) * mid
                            for t in np.linspace(0, 1, 5)
                        ]
                    else:
                        curve = (start, end)
                    lines.InsertNextCell(len(curve))
                    for point in curve:
                        lines.InsertCellPoint(edge_points.InsertNextPoint(*point))
        self.edge_mesh = vtkPolyData()
        self.edge_mesh.SetPoints(edge_points)
        self.edge_mesh.SetLines(lines)
        self.edge_mapper.SetInputData(self.edge_mesh)
        self.edge_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        self.edge_mapper.SetRelativeCoincidentTopologyLineOffsetParameters(-1, -1)
        self.edge_actor.SetVisibility(self.edges_visible)
        self.source_actor.VisibilityOff()
        self.actor.VisibilityOn()
        self.set_selection(())
        if fit:
            self.camera("iso")

    def set_edges(self, visible):
        self.edges_visible = bool(visible)
        self.edge_actor.SetVisibility(visible and self.solid is not None)
        self.render()

    def set_selection(self, identifiers):
        self._selected_faces = set(identifiers)
        self._recolor()

    def set_conditions(self, conditions, load_factor=1.0):
        kinds = {}
        for condition in conditions:
            for identifier in condition.surfaces:
                kinds.setdefault(identifier, set()).add(condition.kind)
        self._face_kinds = {
            i: next(iter(v)) if len(v) == 1 else "multiple" for i, v in kinds.items()
        }
        glyphs, invalid = boundary_glyphs(self._frames, conditions, load_factor)
        self.glyph_layer.set_glyphs(glyphs)
        self._recolor()
        return invalid

    def shutdown(self):
        if hasattr(self, "glyph_layer"):
            self.glyph_layer.close()
        super().shutdown()

    def _recolor(self):
        if not self._cell_faces:
            return
        colors = vtkUnsignedCharArray()
        colors.SetName("Boundary condition / selection")
        colors.SetNumberOfComponents(3)
        for identifier in self._cell_faces:
            kind = (
                "selected"
                if identifier in self._selected_faces
                else "hover"
                if identifier == self._hover
                else self._face_kinds.get(identifier, "plain")
            )
            colors.InsertNextTuple3(*SURFACE_COLORS[kind])
        self.mesh.GetCellData().SetScalars(colors)
        self.render()

    def pick_surface(self, x, y):
        """Qt logical coordinates, converted to the actual render-window size."""
        if self.solid is None or self._closed or not self.selection_enabled:
            return None
        # Keep compass interaction independent of a face behind the compass.
        if x >= self.view.width() - 132 and y < 132:
            return None
        width, height = self.view.GetRenderWindow().GetSize()
        if self.view.width() <= 0 or self.view.height() <= 0:
            return None
        self.picker.Pick(
            x * width / self.view.width(),
            height - 1 - y * height / self.view.height(),
            0,
            self.renderer,
        )
        cell = self.picker.GetCellId()
        return (
            self._cell_faces[cell]
            if self.picker.GetActor() == self.actor
            and 0 <= cell < len(self._cell_faces)
            else None
        )

    def eventFilter(self, obj, event):
        if obj is self.view and not self._closed:
            kind = event.type()
            if kind == QEvent.Type.MouseButtonPress:
                self._press = (
                    event.position()
                    if event.button() == Qt.MouseButton.LeftButton
                    else None
                )
            elif (
                kind == QEvent.Type.MouseButtonRelease
                and event.button() == Qt.MouseButton.LeftButton
            ):
                position = event.position()
                if (
                    self._press is not None
                    and (position - self._press).manhattanLength() <= 4
                    and self.selection_enabled
                ):
                    face = self.pick_surface(position.x(), position.y())
                    if not (
                        position.x() >= self.view.width() - 132 and position.y() < 132
                    ):
                        additive = bool(
                            event.modifiers()
                            & (
                                Qt.KeyboardModifier.ControlModifier
                                | Qt.KeyboardModifier.ShiftModifier
                            )
                        )
                        self.surface_picked.emit(face, additive)
                self._press = None
            elif (
                kind == QEvent.Type.MouseMove
                and event.buttons() == Qt.MouseButton.NoButton
            ):
                face = self.pick_surface(event.position().x(), event.position().y())
                if face != self._hover:
                    self._hover = face
                    self._recolor()
                    self.surface_hovered.emit(face)
            elif kind == QEvent.Type.Leave:
                self._hover = None
                self._recolor()
                self.surface_hovered.emit(None)
        return super().eventFilter(obj, event)

    def fit_scene(self, selection=False):
        if self._closed or self.dragging:
            return
        actor = self.actor if self.solid is not None else self.source_actor
        bounds = actor.GetBounds()
        if (
            bounds
            and np.isfinite(bounds).all()
            and all(bounds[i] <= bounds[i + 1] for i in (0, 2, 4))
        ):
            self.renderer.ResetCamera(bounds)
            self.renderer.GetActiveCamera().Zoom(0.9)
            self.renderer.ResetCameraClippingRange()
            self.render()
