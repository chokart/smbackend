from pydantic import BaseModel
from typing import List, Optional
from enum import Enum

class NodeType(str, Enum):
    ALIMENTO = "initial"
    PROCESO = "process"
    CONCENTRADO = "concentrate"
    RELAVE = "tailing"
    PRODUCTO = "product"

class ElementData(BaseModel):
    name: str
    ley: float
    ley_error: float
    ley_fixed: bool = False
    not_measured: bool = False

class FlowData(BaseModel):
    tonelaje: float
    tonelaje_error: float
    tonelaje_fixed: bool = False
    not_measured: bool = False
    elementos: List[ElementData]

class Flow(BaseModel):
    id: str
    source: str
    target: str
    data: FlowData

class Node(BaseModel):
    id: str
    label: str
    node_type: NodeType

class ReconciliationRequest(BaseModel):
    nodes: List[Node]
    flows: List[Flow]

class CorrectedElementData(BaseModel):
    name: str
    ley_corregida: float
    ley_corregida_error: float

class CorrectedFlowData(BaseModel):
    tonelaje_corregido: float
    tonelaje_corregido_error: float
    elementos: List[CorrectedElementData]

class ReconciledFlow(BaseModel):
    id: str
    corrected_data: CorrectedFlowData

class ReconciliationResponse(BaseModel):
    flows: List[ReconciledFlow]
    report: Optional[dict] = None