from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from shapely.affinity import scale, translate
from shapely.geometry import Polygon

from . import constants as c


@dataclass
class Roi:
    xlims: Optional[Tuple[float, float]] = None
    ylims: Optional[Tuple[float, float]] = None
    center_xy: Optional[Tuple[float, float]] = None
    edge_size: Optional[float] = None
    polygon: Optional[Polygon] = None
    name: str = 'unnamed_roi'
    px_size: float = c.PIXEL_SIZE_MICRONS  # microns per pixel

    def __post_init__(self):

        has_lims = has_center_and_edge = has_polygon = False

        if self.xlims and self.ylims:
            has_lims = True
        elif self.center_xy and self.edge_size:
            has_center_and_edge = True
        elif self.polygon:
            has_polygon = True

        if sum([has_lims, has_center_and_edge, has_polygon]) > 1:
            raise ValueError(
                'Multiple sets of ROI parameters provided. Please provide only one of the following:\n'
                '(1) xlims and ylims, (2) center_xy and edge_size, or (3) polygon.'
            )

        if has_center_and_edge:
            half = self.edge_size / 2
            cx, cy = self.center_xy
            self.xlims = (cx - half, cx + half)
            self.ylims = (cy - half, cy + half)

        elif has_polygon:
            minx, miny, maxx, maxy = self.polygon.bounds
            self.xlims = (minx, maxx)
            self.ylims = (miny, maxy)

    @property
    def width(self):
        return self.xlims[1] - self.xlims[0]

    @property
    def height(self):
        return self.ylims[1] - self.ylims[0]

    @property
    def width_um(self):
        return self.width * self.px_size

    @property
    def height_um(self):
        return self.height * self.px_size

    @property
    def center(self):
        x, y = self.poly.centroid.coords.xy
        return (x[0], y[0])

    @property
    def extent(self):
        return (self.xlims[0], self.xlims[1], self.ylims[0], self.ylims[1])

    @property
    def extent_array(self):
        return np.array(self.extent, dtype=np.int32)

    @property
    def poly(self) -> Polygon:
        x0, x1 = self.xlims
        y0, y1 = self.ylims
        return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])

    def affine(self, scale=1, xoff=0, yoff=0):
        aff_roi = Roi(polygon=self.poly, name=self.name)
        aff_roi = aff_roi.scale(factor=scale)
        aff_roi = aff_roi.translate(xoff=xoff, yoff=yoff)
        return aff_roi

    def scale(self, factor):
        scaled_roi = scale(self.poly, xfact=factor, yfact=factor, origin='center')
        sc_roi = Roi(polygon=scaled_roi, name=None)
        return sc_roi

    def translate(self, xoff=0, yoff=0, relative=True):
        if relative:
            xoff_frac = self.width * xoff
            yoff_frac = self.height * yoff
        else:
            xoff_frac = xoff
            yoff_frac = yoff

        trans_roi = translate(self.poly, xoff=xoff_frac, yoff=yoff_frac)
        tr_roi = Roi(polygon=trans_roi, name=None)
        return tr_roi

    def subtile_roi(self, n=3, label_prefix: str = None, labels='123'):
        scale_factor = 1 / n  # Each sub-ROI should have width and height 1/n of the original

        scaled_roi = self.scale(scale_factor)

        sub_rois = []
        label = 1  # Initialize counter so the top-left sub-ROI becomes 1

        # Reverse the vertical loop: This ensures we start with the top row
        for j in reversed(range(n)):
            for i in range(n):
                # Compute the translation offsets relative to the overall ROI center.
                offset_x = (i - (n - 1) / 2) * (self.width / n)
                offset_y = (j - (n - 1) / 2) * (self.height / n)

                if labels == 'abc':
                    roi_label = chr(65 + (label - 1))
                elif labels == '123':
                    roi_label = str(label)

                sub_roi = scaled_roi.translate(xoff=offset_x, yoff=offset_y, relative=False)
                if label_prefix:
                    roi_label = f'{label_prefix}{roi_label}'
                sub_roi.name = roi_label
                sub_rois.append(sub_roi)

                label += 1

        return sub_rois

    def crop_array(self, array: np.ndarray):

        xlim = (self.extent_array[0:2]).astype(int)
        ylim = (self.extent_array[2:4]).astype(int)

        return array[ylim[0] : ylim[1], xlim[0] : xlim[1]]

    def __repr__(self) -> str:
        return f'Roi {self.name!r}: ({self.width_um:.2f}x{self.height_um:.2f}) microns'
