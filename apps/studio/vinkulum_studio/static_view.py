"""Static mesh display and read-only tables; geometry scale never changes data."""

import numpy as np
from PySide6.QtCore import QAbstractTableModel, Qt
from vtkmodules.vtkCommonCore import vtkDoubleArray, vtkLookupTable, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkHexahedron, vtkUnstructuredGrid
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
from vtkmodules.vtkRenderingCore import vtkActor, vtkDataSetMapper

from .viewport import Viewport


class StaticTable(QAbstractTableModel):
    def __init__(self, headers=(), rows=(), parent=None):
        super().__init__(parent)
        self.headers = tuple(headers)
        self.rows = np.array(rows, dtype=float, copy=True)
        self.rows.setflags(write=False)

    def rowCount(self, parent=None):
        return 0 if parent is not None and parent.isValid() else len(self.rows)

    def columnCount(self, parent=None):
        return 0 if parent is not None and parent.isValid() else len(self.headers)

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        if role not in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.ToolTipRole):
            return None
        value = float(self.rows[index.row(), index.column()])
        if self.headers[index.column()] in ("Node", "Element", "IP"):
            return str(int(value))
        return repr(value) if role == Qt.ItemDataRole.ToolTipRole else f"{value:.7g}"

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole:
            if section < 0 or (
                orientation == Qt.Orientation.Horizontal
                and section >= len(self.headers)
            ):
                return None
            return (
                self.headers[section]
                if orientation == Qt.Orientation.Horizontal
                else str(section + 1)
            )
        return None


class StaticViewport(Viewport):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAccessibleName("Finite-element mesh and displacement field")
        self._editable = False
        self.box.Off()
        self.mesh = vtkUnstructuredGrid()
        self.mapper = vtkDataSetMapper()
        self.mapper.SetInputData(self.mesh)
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(0.57, 0.71, 0.91)
        self.actor.GetProperty().EdgeVisibilityOn()
        self.actor.GetProperty().SetEdgeColor(0.11, 0.16, 0.22)
        self.actor.GetProperty().SetLineWidth(1)
        self.actor.GetProperty().SetSpecular(0.12)
        self.renderer.AddActor(self.actor)
        self.reference_mapper = vtkDataSetMapper()
        self.reference_mesh = vtkUnstructuredGrid()
        self.reference_mapper.SetInputData(self.reference_mesh)
        self.reference_mapper.ScalarVisibilityOff()
        self.reference_actor = vtkActor()
        self.reference_actor.SetMapper(self.reference_mapper)
        self.reference_actor.GetProperty().SetRepresentationToWireframe()
        self.reference_actor.GetProperty().SetColor(0.67, 0.70, 0.76)
        self.reference_actor.GetProperty().SetOpacity(0.3)
        self.renderer.AddActor(self.reference_actor)
        self.legend = vtkScalarBarActor()
        self.legend.SetTitle("|u| [m]")
        self.legend.SetNumberOfLabels(5)
        self.legend.SetLabelFormat("%.3g")
        self.legend.SetWidth(0.12)
        self.legend.SetHeight(0.48)
        self.legend.SetPosition(0.85, 0.16)
        self.legend.UnconstrainedFontSizeOn()
        self.legend.SetVerticalTitleSeparation(12)
        for prop in (
            self.legend.GetTitleTextProperty(),
            self.legend.GetLabelTextProperty(),
        ):
            prop.SetFontFamilyToArial()
            prop.SetColor(0.90, 0.92, 0.96)
            prop.ItalicOff()
            prop.BoldOff()
            prop.ShadowOff()
            prop.SetFontSize(13)
        self.renderer.AddViewProp(self.legend)
        self.legend.VisibilityOff()
        self.display_positions = None
        self.field_values = None

    @staticmethod
    def _grid(study, coordinates):
        points = vtkPoints()
        points.SetDataTypeToDouble()
        for row in coordinates:
            points.InsertNextPoint(*row)
        mesh = vtkUnstructuredGrid()
        mesh.SetPoints(points)
        for nodes in study.elements:
            cell = vtkHexahedron()
            for i, identifier in enumerate(nodes):
                cell.GetPointIds().SetId(i, identifier - 1)
            mesh.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        return mesh

    def set_study(self, study, report=None, *, scale=1.0, fit=False):
        positions = np.array(study.nodes, dtype=float)
        field = None
        if report is not None:
            displacement = np.array(report["displacements"], dtype=float)
            maximum_component = np.abs(displacement).max()
            field = (
                np.linalg.norm(displacement / maximum_component, axis=1)
                * maximum_component
                if maximum_component > 0
                else np.zeros(len(positions))
            )
            displayed = positions + scale * displacement
            if not np.isfinite(displayed).all() or not np.isfinite(field).all():
                raise ValueError(
                    "The displayed deformation exceeds the numeric range. Reduce its multiplier."
                )
        else:
            displayed = positions
        self.reference_mesh = self._grid(study, positions)
        self.reference_mapper.SetInputData(self.reference_mesh)
        self.field_values = field
        self.display_positions = displayed.copy()
        self.display_positions.setflags(write=False)
        self.mesh = self._grid(study, displayed)
        self.mapper.SetInputData(self.mesh)
        self.mapper.SetScalarVisibility(report is not None)
        self.legend.SetVisibility(report is not None)
        self.reference_actor.SetVisibility(report is not None and scale != 0)
        if report is not None:
            scalar = vtkDoubleArray()
            scalar.SetName("Displacement magnitude [m]")
            for value in self.field_values:
                scalar.InsertNextValue(float(value))
            self.mesh.GetPointData().SetScalars(scalar)
            maximum = float(self.field_values.max())
            self.legend.SetVisibility(maximum > 0)
            # The degenerate zero-field mapper needs a nonempty range, but its
            # legend is hidden. Actual nodal values remain exactly zero.
            extent = maximum if maximum > 0 else 1.0
            colors = np.array(
                (
                    (0.267, 0.005, 0.329),
                    (0.230, 0.322, 0.546),
                    (0.128, 0.567, 0.551),
                    (0.369, 0.789, 0.383),
                    (0.993, 0.906, 0.144),
                )
            )
            lookup = vtkLookupTable()
            lookup.SetNumberOfTableValues(256)
            lookup.SetRange(0, extent)
            for i, t in enumerate(np.linspace(0, 4, 256)):
                left = min(int(t), 3)
                color = colors[left] * (left + 1 - t) + colors[left + 1] * (t - left)
                lookup.SetTableValue(i, *color, 1)
            lookup.Build()
            self.mapper.SetLookupTable(lookup)
            self.mapper.SetScalarRange(0, extent)
            self.legend.SetLookupTable(lookup)
        if fit:
            self.camera("iso")
        else:
            self.renderer.ResetCameraClippingRange()
            self.render()
