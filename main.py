import logging
from typing import Optional, Dict, Any, List
from pydantic import BaseModel
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from models import ReconciliationRequest, ReconciliationResponse
from reconciliation import reconcile_data_qp, _get_balance_system, format_equations_to_strings
from water_models import WaterReconciliationRequest, WaterReconciliationResponse
from water_reconciliation import reconcile_water_data, _format_physical_equations
from datetime import datetime
from pulp_models import PulpCalculationRequest, PulpCalculationResponse
from pulp_logic import calculate_pulp_properties
from level_models import LevelReadingCreate, LevelReadingResponse
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

# --- ENDPOINTS TELEMETRÍA SENSOR DE NIVEL (ESP32) ---
@app.post("/level/reading", response_model=LevelReadingResponse)
def receive_level_reading(reading: LevelReadingCreate):
    """
    Recibe la medición transmitida por el ESP32 vía WiFi HTTP POST.
    """
    try:
        saved = db.save_level_reading(
            device_id=reading.device_id or "esp32_tanque_1",
            distance_cm=round(reading.distance_cm, 2),
            level_cm=round(reading.level_cm, 2),
            percentage=round(reading.percentage, 2)
        )
        return LevelReadingResponse(
            id=saved["id"],
            device_id=saved["device_id"],
            distance_cm=saved["distance_cm"],
            level_cm=saved["level_cm"],
            percentage=saved["percentage"],
            created_at=saved["created_at"],
            is_online=True,
            seconds_ago=0.0
        )
    except Exception as e:
        logger.exception("Error al registrar lectura de nivel:")
        raise HTTPException(status_code=500, detail=f"Error guardando telemetría: {e}")

@app.get("/level/latest", response_model=Optional[LevelReadingResponse])
def get_latest_level(device_id: Optional[str] = None):
    """
    Obtiene la última lectura del sensor de nivel y determina si el ESP32 está online.
    """
    latest = db.get_latest_level_reading(device_id=device_id)
    if not latest:
        return None
    
    # Evaluar estado de conexión según la antigüedad de la última lectura (umbral: 15 seg)
    seconds_ago = None
    is_online = False
    try:
        created_dt = datetime.fromisoformat(latest["created_at"])
        delta = datetime.now() - created_dt
        seconds_ago = max(0.0, delta.total_seconds())
        is_online = seconds_ago <= 15.0
    except Exception:
        pass

    return LevelReadingResponse(
        id=latest["id"],
        device_id=latest["device_id"],
        distance_cm=latest["distance_cm"],
        level_cm=latest["level_cm"],
        percentage=latest["percentage"],
        created_at=latest["created_at"],
        is_online=is_online,
        seconds_ago=round(seconds_ago, 1) if seconds_ago is not None else None
    )

@app.get("/level/history", response_model=List[LevelReadingResponse])
def get_level_history(device_id: Optional[str] = None, limit: int = 50):
    """
    Historial cronológico de lecturas para visualización gráfica en el frontend.
    """
    history = db.get_level_history(device_id=device_id, limit=limit)
    res = []
    now = datetime.now()
    for row in history:
        sec = None
        try:
            dt = datetime.fromisoformat(row["created_at"])
            sec = max(0.0, (now - dt).total_seconds())
        except Exception:
            pass
        res.append(LevelReadingResponse(
            id=row["id"],
            device_id=row["device_id"],
            distance_cm=row["distance_cm"],
            level_cm=row["level_cm"],
            percentage=row["percentage"],
            created_at=row["created_at"],
            is_online=(sec <= 15.0 if sec is not None else False),
            seconds_ago=round(sec, 1) if sec is not None else None
        ))
    return res

@app.post("/level/simulate", response_model=LevelReadingResponse)
def simulate_level_reading(reading: LevelReadingCreate):
    """
    Permite enviar una lectura simulada desde el frontend para probar la UI sin el ESP32 físico.
    """
    return receive_level_reading(reading)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
