from pydantic import BaseModel
from typing import List, Optional
from enum import Enum

class NodeType(str, Enum):
    ALIMENTO = "initial"
    PROCESO = "process"
    CONCENTRADO = "concentrate"
    RELAVE = "tailing"
    PRODUCTO = "product"
    AGUA = "water_source" # Opcional: para fuentes de agua pura

class WaterFlowData(BaseModel):
    # En balance de agua, ahora ingresamos el Flujo Másico de Sólido Seco (tph)
    # y el Porcentaje de Sólidos (%S). La pulpa y el agua se calculan por balance.
    tonelaje_solido: Optional[float] = None # Flujo másico de sólido seco (tph)
    tonelaje_solido_error: Optional[float] = None # Desviación estándar
    tonelaje_solido_fixed: bool = False
    
    # Campos legacy opcionales para compatibilidad con diagramas guardados anteriormente
    tonelaje_pulpa: Optional[float] = None
    tonelaje_pulpa_error: Optional[float] = None
    tonelaje_pulpa_fixed: Optional[bool] = None

    not_measured: bool = False
    
    porcentaje_solidos: float # % de sólidos en peso (Cp)
    porcentaje_solidos_error: float
    porcentaje_solidos_fixed: bool = False
    porcentaje_solidos_not_measured: bool = False

    def get_solid_tonnage(self) -> float:
        if self.tonelaje_solido is not None:
            return float(self.tonelaje_solido)
        if self.tonelaje_pulpa is not None:
            cp = self.porcentaje_solidos or 30.0
            return float(self.tonelaje_pulpa * (cp / 100.0))
        return 0.0

    def get_solid_error(self) -> float:
        if self.tonelaje_solido_error is not None:
            return float(self.tonelaje_solido_error)
        if self.tonelaje_pulpa_error is not None:
            cp = self.porcentaje_solidos or 30.0
            return float(self.tonelaje_pulpa_error * (cp / 100.0))
        return 0.1

    def is_solid_fixed(self) -> bool:
        if self.tonelaje_solido_fixed:
            return True
        if self.tonelaje_pulpa_fixed:
            return True
        return False

class WaterFlow(BaseModel):
    id: str
    source: str
    target: str
    data: WaterFlowData

class WaterNode(BaseModel):
    id: str
    label: str
    node_type: NodeType

class WaterReconciliationRequest(BaseModel):
    nodes: List[WaterNode]
    flows: List[WaterFlow]

class CorrectedWaterFlowData(BaseModel):
    tonelaje_solido_corregido: float
    tonelaje_solido_error: float
    
    porcentaje_solidos_corregido: float
    porcentaje_solidos_error: float
    
    # Datos calculados a partir de sólido seco y % sólidos
    tonelaje_pulpa_calculado: float
    tonelaje_agua_calculado: float

    # Propiedades legacy para compatibilidad con interfaces o consumidores anteriores
    tonelaje_pulpa_corregido: Optional[float] = None
    tonelaje_pulpa_error: Optional[float] = None
    tonelaje_solido_calculado: Optional[float] = None

class ReconciledWaterFlow(BaseModel):
    id: str
    corrected_data: CorrectedWaterFlowData

class WaterReconciliationResponse(BaseModel):
    flows: List[ReconciledWaterFlow]
    report: Optional[dict] = None
