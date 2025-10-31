from typing import ClassVar as _ClassVar
from typing import Iterable as _Iterable
from typing import Mapping as _Mapping
from typing import Optional as _Optional
from typing import Union as _Union

from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf.internal import containers as _containers

DESCRIPTOR: _descriptor.FileDescriptor

class UmapEntry(_message.Message):
    __slots__ = ("umapY", "umapX")
    UMAPY_FIELD_NUMBER: _ClassVar[int]
    UMAPX_FIELD_NUMBER: _ClassVar[int]
    umapY: float
    umapX: float
    def __init__(self, umapY: _Optional[float] = ..., umapX: _Optional[float] = ...) -> None: ...

class SingleMask(_message.Message):
    __slots__ = (
        "vertices",
        "color",
        "area",
        "totalCounts",
        "totalGenes",
        "cellId",
        "clusterId",
        "proteinValues",
        "nonzeroGeneIndices",
        "nonzeroGeneValues",
        "umapValues",
    )
    VERTICES_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    AREA_FIELD_NUMBER: _ClassVar[int]
    TOTALCOUNTS_FIELD_NUMBER: _ClassVar[int]
    TOTALGENES_FIELD_NUMBER: _ClassVar[int]
    CELLID_FIELD_NUMBER: _ClassVar[int]
    CLUSTERID_FIELD_NUMBER: _ClassVar[int]
    PROTEINVALUES_FIELD_NUMBER: _ClassVar[int]
    NONZEROGENEINDICES_FIELD_NUMBER: _ClassVar[int]
    NONZEROGENEVALUES_FIELD_NUMBER: _ClassVar[int]
    UMAPVALUES_FIELD_NUMBER: _ClassVar[int]
    vertices: _containers.RepeatedScalarFieldContainer[int]
    color: _containers.RepeatedScalarFieldContainer[int]
    area: int
    totalCounts: int
    totalGenes: int
    cellId: str
    clusterId: str
    proteinValues: _containers.RepeatedScalarFieldContainer[int]
    nonzeroGeneIndices: _containers.RepeatedScalarFieldContainer[int]
    nonzeroGeneValues: _containers.RepeatedScalarFieldContainer[int]
    umapValues: UmapEntry
    def __init__(
        self,
        vertices: _Optional[_Iterable[int]] = ...,
        color: _Optional[_Iterable[int]] = ...,
        area: _Optional[int] = ...,
        totalCounts: _Optional[int] = ...,
        totalGenes: _Optional[int] = ...,
        cellId: _Optional[str] = ...,
        clusterId: _Optional[str] = ...,
        proteinValues: _Optional[_Iterable[int]] = ...,
        nonzeroGeneIndices: _Optional[_Iterable[int]] = ...,
        nonzeroGeneValues: _Optional[_Iterable[int]] = ...,
        umapValues: _Optional[_Union[UmapEntry, _Mapping]] = ...,
    ) -> None: ...

class ColormapEntry(_message.Message):
    __slots__ = ("clusterId", "color")
    CLUSTERID_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    clusterId: str
    color: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, clusterId: _Optional[str] = ..., color: _Optional[_Iterable[int]] = ...) -> None: ...

class Metadata(_message.Message):
    __slots__ = ("proteinNames", "geneNames")
    PROTEINNAMES_FIELD_NUMBER: _ClassVar[int]
    GENENAMES_FIELD_NUMBER: _ClassVar[int]
    proteinNames: _containers.RepeatedScalarFieldContainer[str]
    geneNames: _containers.RepeatedScalarFieldContainer[str]
    def __init__(
        self, proteinNames: _Optional[_Iterable[str]] = ..., geneNames: _Optional[_Iterable[str]] = ...
    ) -> None: ...

class CellMasks(_message.Message):
    __slots__ = ("cellMasks", "colormap", "metadata", "numberOfCells")
    CELLMASKS_FIELD_NUMBER: _ClassVar[int]
    COLORMAP_FIELD_NUMBER: _ClassVar[int]
    METADATA_FIELD_NUMBER: _ClassVar[int]
    NUMBEROFCELLS_FIELD_NUMBER: _ClassVar[int]
    cellMasks: _containers.RepeatedCompositeFieldContainer[SingleMask]
    colormap: _containers.RepeatedCompositeFieldContainer[ColormapEntry]
    metadata: Metadata
    numberOfCells: int
    def __init__(
        self,
        cellMasks: _Optional[_Iterable[_Union[SingleMask, _Mapping]]] = ...,
        colormap: _Optional[_Iterable[_Union[ColormapEntry, _Mapping]]] = ...,
        metadata: _Optional[_Union[Metadata, _Mapping]] = ...,
        numberOfCells: _Optional[int] = ...,
    ) -> None: ...
