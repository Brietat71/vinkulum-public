"""Batched, non-pickable boundary symbols with a stable logical-pixel scale."""

import math

import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkFiltersCore import vtkTubeFilter
from vtkmodules.vtkFiltersSources import vtkArrowSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkGlyph3DMapper


def line_source(segments):
    points, lines = vtkPoints(), vtkCellArray()
    for start, end in segments:
        lines.InsertNextCell(2)
        lines.InsertCellPoint(points.InsertNextPoint(*start))
        lines.InsertCellPoint(points.InsertNextPoint(*end))
    data = vtkPolyData()
    data.SetPoints(points)
    data.SetLines(lines)
    # Physical tube geometry keeps strokes visible on drivers that clamp GL
    # line widths. The glyph scale makes its diameter 1.5 logical pixels.
    tube = vtkTubeFilter()
    tube.SetInputData(data)
    tube.SetRadius(0.025)
    tube.SetNumberOfSides(8)
    tube.CappingOn()
    tube.Update()
    return tube.GetOutput()


def glyph_mapper(source):
    mapper = vtkGlyph3DMapper()
    mapper.SetSourceData(source)
    mapper.SetOrientationArray("direction")
    mapper.SetOrientationModeToDirection()
    mapper.SetScaleArray("symbol_scale")
    mapper.SetScaleModeToScaleByMagnitude()
    mapper.SetScaleFactor(1)
    mapper.ScalarVisibilityOff()
    return mapper


class BoundaryGlyphLayer:
    def __init__(self, viewport, colors):
        self.viewport = viewport
        self.renderer = viewport.renderer
        self.visible = True
        self.groups = {}
        self.glyphs = ()
        self._signature = None
        self.arrow = vtkArrowSource()
        self.arrow.SetTipLength(0.28)
        self.arrow.SetTipRadius(0.10)
        self.arrow.SetShaftRadius(0.025)
        self.arrow.SetTipResolution(12)
        self.arrow.SetShaftResolution(8)
        self.arrow.Update()
        pressure = [((0, 0, 0), (1, 0, 0))]
        for tip in (
            (0.72, 0.13, 0),
            (0.72, -0.13, 0),
            (0.72, 0, 0.13),
            (0.72, 0, -0.13),
        ):
            pressure.append(((1, 0, 0), tip))
        support = [((-0.38, 0, 0), (0.38, 0, 0))]
        for x in (-0.38, 0.38):
            support.extend(
                (((x, -0.13, 0), (x, 0.13, 0)), ((x, 0, -0.13), (x, 0, 0.13)))
            )
        for kind, source in (
            ("pressure", line_source(pressure)),
            ("total_force", self.arrow.GetOutput()),
            ("support", line_source(support)),
        ):
            mapper = glyph_mapper(source)
            actor = vtkActor()
            actor.SetMapper(mapper)
            actor.PickableOff()
            # Symbols must not change the camera's model framing.
            actor.UseBoundsOff()
            actor.GetProperty().SetColor(*(0.55 + 0.45 * v / 255 for v in colors[kind]))
            # These are semantic annotations: lighting must not wash out their
            # colour as the camera or surface normal changes.
            actor.GetProperty().LightingOff()
            self.renderer.AddActor(actor)
            self.groups[kind] = {
                "actor": actor,
                "mapper": mapper,
                "source": source,
            }
        self.set_glyphs(())
        self._observer = self.renderer.AddObserver("StartEvent", self._rescale)

    def set_glyphs(self, glyphs):
        self.glyphs = tuple(glyphs)
        for kind, group in self.groups.items():
            rows = [g for g in self.glyphs if g.kind == kind]
            points = np.array([g.point for g in rows], dtype=float).reshape(-1, 3)
            directions = np.array([g.direction for g in rows], dtype=float).reshape(
                -1, 3
            )
            group.update(
                rows=tuple(rows),
                anchors=points.copy(),
                positions=points,
                normals=np.array([g.normal for g in rows], dtype=float).reshape(-1, 3),
                directions=directions,
                scales=np.ones(len(rows)),
            )
            vtk_points = vtkPoints()
            vtk_points.SetData(numpy_to_vtk(group["positions"], deep=False))
            data = vtkPolyData()
            data.SetPoints(vtk_points)
            for key, name in (("directions", "direction"), ("scales", "symbol_scale")):
                array = numpy_to_vtk(group[key], deep=False)
                array.SetName(name)
                data.GetPointData().AddArray(array)
            group["data"] = data
            group["mapper"].SetInputData(data)
            group["actor"].SetVisibility(self.visible and bool(rows))
        self._signature = None

    def set_visible(self, visible):
        self.visible = bool(visible)
        for group in self.groups.values():
            group["actor"].SetVisibility(self.visible and bool(group["rows"]))
        self._signature = None
        self.viewport.render()

    def _rescale(self, *_):
        if not self.visible or not self.glyphs or self.viewport._closed:
            return
        camera = self.renderer.GetActiveCamera()
        view = self.viewport.view
        horizontal = camera.GetUseHorizontalViewAngle()
        pixels = view.width() if horizontal else view.height()
        if pixels <= 0:
            return
        signature = (
            camera.GetParallelProjection(),
            camera.GetParallelScale(),
            camera.GetViewAngle(),
            camera.GetPosition(),
            camera.GetDirectionOfProjection(),
            horizontal,
            view.width(),
            view.height(),
            view.devicePixelRatioF(),
        )
        if signature == self._signature:
            return
        self._signature = signature
        position = np.array(camera.GetPosition())
        direction = np.array(camera.GetDirectionOfProjection())
        boxes = []
        for kind, group in self.groups.items():
            if not group["rows"]:
                continue
            logical_length = 30 if kind != "support" else 24
            if camera.GetParallelProjection():
                sizes = np.full(
                    len(group["rows"]),
                    logical_length
                    * 2
                    * camera.GetParallelScale()
                    / max(1, view.height()),
                )
            else:
                depths = (group["anchors"] - position) @ direction
                sizes = (
                    logical_length
                    * 2
                    * np.maximum(depths, 0)
                    * math.tan(math.radians(camera.GetViewAngle()) / 2)
                    / pixels
                )
            group["scales"][:] = sizes
            # Keep the symbol outside the body. Inward arrows end at the face;
            # outward arrows start there. A small normal gap prevents z-fighting.
            offsets = (0.45 if kind == "support" else 0.06) * group["normals"]
            if kind != "support":
                inward = np.sum(group["directions"] * group["normals"], axis=1) < 0
                offsets = offsets - inward[:, None] * group["directions"]
            group["positions"][:] = group["anchors"] + sizes[:, None] * offsets
            group["data"].GetPoints().GetData().Modified()
            group["data"].GetPoints().Modified()
            group["data"].GetPointData().GetArray("symbol_scale").Modified()
            # Include symbol extents in clipping, while excluding them from fit.
            boxes.append(
                (
                    np.min(group["anchors"] - 2 * sizes[:, None], axis=0),
                    np.max(group["anchors"] + 2 * sizes[:, None], axis=0),
                )
            )
        bounds = self.renderer.ComputeVisiblePropBounds()
        if boxes and all(bounds[i] <= bounds[i + 1] for i in (0, 2, 4)):
            lo = np.minimum(
                np.array(bounds)[::2], np.min([b[0] for b in boxes], axis=0)
            )
            hi = np.maximum(
                np.array(bounds)[1::2], np.max([b[1] for b in boxes], axis=0)
            )
            self.renderer.ResetCameraClippingRange(
                tuple(np.column_stack((lo, hi)).flat)
            )

    def close(self):
        if self._observer is not None:
            self.renderer.RemoveObserver(self._observer)
            self._observer = None
