import logging
from typing import Optional, Dict, Any, List
from pydantic import BaseModel
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from models import ReconciliationRequest, ReconciliationResponse
from reconciliation import reconcile_data_qp, _get_balance_system, format_equations_to_strings
from water_models import WaterReconciliationRequest, WaterReconciliationResponse
from water_reconciliation import reconcile_water_data, _format_physical_equations
from pulp_models import PulpCalculationRequest, PulpCalculationResponse
from pulp_logic import calculate_pulp_properties
import database as db

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("suiteminerals")

app = FastAPI(
    title="Suite Minerals API",
    version="1.0.0",
    description="Motor de cálculo, balance metalúrgico y reconciliación de datos con persistencia SQLite",
)

# Configuración de CORS estandarizada para consumo web
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- ESQUEMAS PARA DIAGRAMAS ---
class DiagramCreateRequest(BaseModel):
    name: str
    description: Optional[str] = ""
    module: Optional[str] = "metallurgical"
    data: Dict[str, Any]

class DiagramUpdateRequest(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    data: Optional[Dict[str, Any]] = None

@app.get("/health")
def health_check():
    """Endpoint de salud para balanceadores de carga y Dokploy."""
    return {"status": "ok"}

@app.get("/")
def read_root():
    return {"message": "Servidor de Balance Metalúrgico Suite Minerals funcionando", "version": "1.0.0"}

# --- ENDPOINTS DE DIAGRAMAS (SQLITE) ---
@app.get("/diagrams")
def list_saved_diagrams(module: Optional[str] = None):
    try:
        return db.list_diagrams(module=module)
    except Exception as e:
        logger.exception("Error al listar diagramas:")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/diagrams/{diagram_id}")
def get_saved_diagram(diagram_id: int):
    diagram = db.get_diagram_by_id(diagram_id)
    if not diagram:
        raise HTTPException(status_code=404, detail="Diagrama no encontrado")
    return diagram

@app.post("/diagrams")
def create_saved_diagram(payload: DiagramCreateRequest):
    try:
        return db.create_diagram(
            name=payload.name,
            description=payload.description or "",
            module=payload.module or "metallurgical",
            data=payload.data,
        )
    except Exception as e:
        logger.exception("Error al guardar diagrama:")
        raise HTTPException(status_code=500, detail=str(e))

@app.put("/diagrams/{diagram_id}")
def update_saved_diagram(diagram_id: int, payload: DiagramUpdateRequest):
    try:
        updated = db.update_diagram(
            diagram_id=diagram_id,
            name=payload.name,
            description=payload.description,
            data=payload.data,
        )
        if not updated:
            raise HTTPException(status_code=404, detail="Diagrama no encontrado")
        return updated
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("Error al actualizar diagrama:")
        raise HTTPException(status_code=500, detail=str(e))

@app.delete("/diagrams/{diagram_id}")
def delete_saved_diagram(diagram_id: int):
    try:
        success = db.delete_diagram(diagram_id)
        if not success:
            raise HTTPException(status_code=404, detail="Diagrama no encontrado")
        return {"success": True, "message": "Diagrama eliminado"}
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("Error al eliminar diagrama:")
        raise HTTPException(status_code=500, detail=str(e))

# --- ENDPOINTS DE RECONCILIACIÓN Y CÁLCULOS ---
@app.post("/reconcile", response_model=ReconciliationResponse)
def reconcile(request: ReconciliationRequest):
    try:
        reconciled_flows, report = reconcile_data_qp(request)
        return ReconciliationResponse(flows=reconciled_flows, report=report)
    except ValueError as e:
        logger.warning(f"Error de validación en reconciliación: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Error inesperado durante la optimización:")
        raise HTTPException(status_code=500, detail=f"Error interno en el servidor: {e}")

@app.post("/reconcile/water", response_model=WaterReconciliationResponse)
def reconcile_water(request: WaterReconciliationRequest):
    try:
        reconciled_flows, report = reconcile_water_data(request)
        return WaterReconciliationResponse(flows=reconciled_flows, report=report)
    except ValueError as e:
        logger.warning(f"Error de validación en balance de agua: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Error inesperado en balance de agua:")
        raise HTTPException(status_code=500, detail=f"Error interno en el servidor: {e}")

@app.post("/preview-equations")
def preview_equations(request: ReconciliationRequest):
    try:
        A, b, var_map, _, _, _ = _get_balance_system(request)
        equations = format_equations_to_strings(A, b, var_map)
        return {"equations": equations}
    except Exception as e:
        logger.exception("Error al generar ecuaciones:")
        raise HTTPException(status_code=500, detail=f"Error al generar ecuaciones: {e}")

@app.post("/preview-water-equations")
def preview_water_equations(request: WaterReconciliationRequest):
    try:
        equations = _format_physical_equations(request)
        return {"equations": equations}
    except Exception as e:
        logger.exception("Error al generar ecuaciones de agua:")
        raise HTTPException(status_code=500, detail=f"Error al generar ecuaciones de agua: {e}")

@app.post("/pulp-calculate", response_model=PulpCalculationResponse)
def pulp_calculate(request: PulpCalculationRequest):
    try:
        return calculate_pulp_properties(request)
    except ValueError as e:
        logger.warning(f"Error en datos de pulpa: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Error en el cálculo de pulpa:")
        raise HTTPException(status_code=500, detail=f"Error en el cálculo de pulpa: {e}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
