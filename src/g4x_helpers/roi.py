from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
from shapely.affinity import scale, translate
from shapely.geometry import Polygon

from . import c


@dataclass
class Roi:
    xlims: Optional[Tuple[float, float]] = None
    ylims: Optional[Tuple[float, float]] = None
    extent: Optional[List[float]] = None
    center: Optional[Tuple[float, float]] = None
    edge_size: Optional[float] = None
    polygon: Optional[Polygon] = None
    name: str = 'unnamed_roi'
    px_size: float = c.PIXEL_SIZE_MICRONS  # microns per pixel

    def __post_init__(self):
        # Case 1: build from center + edge_size
        if (
            self.center is not None
            and self.edge_size is not None
            and self.polygon is None
            and not (self.xlims and self.ylims)
            and not self.extent
        ):
            half = self.edge_size / 2
            cx, cy = self.center
            self.xlims = (cx - half, cx + half)
            self.ylims = (cy - half, cy + half)
            self.polygon = self._make_polygon()

        # Case 2: build from xlims/ylims or extent
        if (self.xlims and self.ylims) or self.extent:
            xl = self.xlims if self.xlims else tuple(self.extent[0:2])
            yl = self.ylims if self.ylims else tuple(self.extent[2:4])
            self.xlims, self.ylims = xl, yl
            if self.polygon is None:
                self.polygon = self._make_polygon()

        # Case 3: build from provided polygon
        elif self.polygon is not None and not (self.xlims and self.ylims):
            minx, miny, maxx, maxy = self.polygon.bounds
            self.xlims = (minx, maxx)
            self.ylims = (miny, maxy)

        # Finalize computed properties
        self.width = self.xlims[1] - self.xlims[0]
        self.height = self.ylims[1] - self.ylims[0]
        self.width_um = self.width * self.px_size
        self.height_um = self.height * self.px_size
        self.center = ((self.xlims[0] + self.xlims[1]) / 2, (self.ylims[0] + self.ylims[1]) / 2)
        self.extent_tuple = (self.xlims[0], self.xlims[1], self.ylims[0], self.ylims[1])
        self.extent_array = np.array(self.extent_tuple, dtype=np.int32)
        self.lims = (self.xlims, self.ylims)
        self.order = 'xy'

    def _make_polygon(self) -> Polygon:
        x0, x1 = self.xlims
        y0, y1 = self.ylims
        return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])

    def affine(self, scale=1, xoff=0, yoff=0):
        aff_roi = Roi(polygon=self.polygon, name=self.name)
        aff_roi = aff_roi.scale(factor=scale)
        aff_roi = aff_roi.translate(xoff=xoff, yoff=yoff)
        return aff_roi

    def scale(self, factor):
        scaled_roi = scale(self.polygon, xfact=factor, yfact=factor, origin='center')
        sc_roi = Roi(polygon=scaled_roi, name=None)
        return sc_roi

    def translate(self, xoff=0, yoff=0, relative=True):
        if relative:
            xoff_frac = self.width * xoff
            yoff_frac = self.height * yoff
        else:
            xoff_frac = xoff
            yoff_frac = yoff

        trans_roi = translate(self.polygon, xoff=xoff_frac, yoff=yoff_frac)
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
