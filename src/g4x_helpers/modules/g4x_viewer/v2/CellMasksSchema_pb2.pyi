from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import (
    ClassVar as _ClassVar,
    Iterable as _Iterable,
    Mapping as _Mapping,
    Optional as _Optional,
    Union as _Union,
)

DESCRIPTOR: _descriptor.FileDescriptor

class UmapEntry(_message.Message):
    __slots__ = ('umapY', 'umapX')
    UMAPY_FIELD_NUMBER: _ClassVar[int]
    UMAPX_FIELD_NUMBER: _ClassVar[int]
    umapY: float
    umapX: float
    def __init__(self, umapY: _Optional[float] = ..., umapX: _Optional[float] = ...) -> None: ...

class SingleMask(_message.Message):
    __slots__ = (
        'vertices',
        'color',
        'area',
        'totalCounts',
        'totalGenes',
        'cellId',
        'clusterId',
        'proteins',
        'umapValues',
    )
    class ProteinsEntry(_message.Message):
        __slots__ = ('key', 'value')
        KEY_FIELD_NUMBER: _ClassVar[int]
        VALUE_FIELD_NUMBER: _ClassVar[int]
        key: str
        value: float
        def __init__(self, key: _Optional[str] = ..., value: _Optional[float] = ...) -> None: ...

    VERTICES_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    AREA_FIELD_NUMBER: _ClassVar[int]
    TOTALCOUNTS_FIELD_NUMBER: _ClassVar[int]
    TOTALGENES_FIELD_NUMBER: _ClassVar[int]
    CELLID_FIELD_NUMBER: _ClassVar[int]
    CLUSTERID_FIELD_NUMBER: _ClassVar[int]
    PROTEINS_FIELD_NUMBER: _ClassVar[int]
    UMAPVALUES_FIELD_NUMBER: _ClassVar[int]
    vertices: _containers.RepeatedScalarFieldContainer[float]
    color: _containers.RepeatedScalarFieldContainer[int]
    area: str
    totalCounts: str
    totalGenes: str
    cellId: str
    clusterId: str
    proteins: _containers.ScalarMap[str, float]
    umapValues: UmapEntry
    def __init__(
        self,
        vertices: _Optional[_Iterable[float]] = ...,
        color: _Optional[_Iterable[int]] = ...,
        area: _Optional[str] = ...,
        totalCounts: _Optional[str] = ...,
        totalGenes: _Optional[str] = ...,
        cellId: _Optional[str] = ...,
        clusterId: _Optional[str] = ...,
        proteins: _Optional[_Mapping[str, float]] = ...,
        umapValues: _Optional[_Union[UmapEntry, _Mapping]] = ...,
    ) -> None: ...

class ColormapEntry(_message.Message):
    __slots__ = ('clusterId', 'color')
    CLUSTERID_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    clusterId: str
    color: _containers.RepeatedScalarFieldContainer[int]
    def __init__(self, clusterId: _Optional[str] = ..., color: _Optional[_Iterable[int]] = ...) -> None: ...

class CellMasks(_message.Message):
    __slots__ = ('cellMasks', 'colormap', 'numberOfCells')
    CELLMASKS_FIELD_NUMBER: _ClassVar[int]
    COLORMAP_FIELD_NUMBER: _ClassVar[int]
    NUMBEROFCELLS_FIELD_NUMBER: _ClassVar[int]
    cellMasks: _containers.RepeatedCompositeFieldContainer[SingleMask]
    colormap: _containers.RepeatedCompositeFieldContainer[ColormapEntry]
    numberOfCells: int
    def __init__(
        self,
        cellMasks: _Optional[_Iterable[_Union[SingleMask, _Mapping]]] = ...,
        colormap: _Optional[_Iterable[_Union[ColormapEntry, _Mapping]]] = ...,
        numberOfCells: _Optional[int] = ...,
    ) -> None: ...
