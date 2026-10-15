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

    # Campos para flujos de solo agua (Agua de dilución / adición)
    is_water_only: bool = False
    tonelaje_agua: Optional[float] = None
    tonelaje_agua_error: Optional[float] = None
    tonelaje_agua_fixed: bool = False

    not_measured: bool = False
    
    porcentaje_solidos: float = 0.0 # % de sólidos en peso (Cp)
    porcentaje_solidos_error: float = 0.5
    porcentaje_solidos_fixed: bool = False
    porcentaje_solidos_not_measured: bool = False

    def get_solid_tonnage(self) -> float:
        if self.is_water_only:
            return 0.0
        if self.tonelaje_solido is not None:
            return float(self.tonelaje_solido)
        if self.tonelaje_pulpa is not None:
            cp = self.porcentaje_solidos or 30.0
            return float(self.tonelaje_pulpa * (cp / 100.0))
        return 0.0

    def get_solid_error(self) -> float:
        if self.is_water_only:
            return 0.0
        if self.tonelaje_solido_error is not None and self.tonelaje_solido_error > 0:
            return float(self.tonelaje_solido_error)
        if self.tonelaje_pulpa_error is not None and self.tonelaje_pulpa_error > 0:
            cp = self.porcentaje_solidos or 30.0
            return float(self.tonelaje_pulpa_error * (cp / 100.0))
        s_val = self.get_solid_tonnage()
        return max(s_val * 0.02 + 0.1, 0.1)

    def is_solid_fixed(self) -> bool:
        if self.is_water_only:
            return True
        if self.tonelaje_solido_fixed:
            return True
        if self.tonelaje_pulpa_fixed:
            return True
        return False

    def get_water_tonnage(self) -> float:
        if self.is_water_only:
            return float(self.tonelaje_agua or 0.0)
        s = self.get_solid_tonnage()
        cp = self.porcentaje_solidos or 30.0
        if cp > 0:
            return float(s * (100.0 / cp - 1.0))
        return 0.0

    def get_water_error(self) -> float:
        if self.is_water_only:
            if self.tonelaje_agua_error is not None and self.tonelaje_agua_error > 0:
                return float(self.tonelaje_agua_error)
            w = self.get_water_tonnage()
            return max(w * 0.02 + 0.1, 0.5)
        return 0.5

    def is_water_fixed(self) -> bool:
        if self.is_water_only:
            return bool(self.tonelaje_agua_fixed)
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
    tonelaje_agua_error: Optional[float] = None

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
