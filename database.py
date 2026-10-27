import sqlite3
import json
import os
from datetime import datetime
from typing import List, Optional, Dict, Any

DB_DIR = os.getenv("DB_DIR", os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(DB_DIR, "diagrams.db")

def get_connection():
    """Obtiene una conexión a la base de datos SQLite con row_factory para diccionarios y timeout extendido."""
    conn = sqlite3.connect(DB_PATH, timeout=15.0)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    """Inicializa las tablas de diagramas y lecturas de sensores si no existen."""
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS diagrams (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                description TEXT DEFAULT '',
                module TEXT DEFAULT 'metallurgical',
                data TEXT NOT NULL,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS level_readings (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                device_id TEXT NOT NULL,
                distance_cm REAL NOT NULL,
                level_cm REAL NOT NULL,
                percentage REAL NOT NULL,
                created_at TEXT NOT NULL
            )
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_level_readings_dev_date 
            ON level_readings(device_id, created_at DESC)
        """)
        conn.commit()

def list_diagrams(module: Optional[str] = None) -> List[Dict[str, Any]]:
    """Lista todos los diagramas (resumen sin la data pesada)."""
    with get_connection() as conn:
        cursor = conn.cursor()
        if module:
            cursor.execute("""
                SELECT id, name, description, module, created_at, updated_at
                FROM diagrams
                WHERE module = ?
                ORDER BY updated_at DESC
            """, (module,))
        else:
            cursor.execute("""
                SELECT id, name, description, module, created_at, updated_at
                FROM diagrams
                ORDER BY updated_at DESC
            """)
        rows = cursor.fetchall()
        return [dict(row) for row in rows]

def get_diagram_by_id(diagram_id: int) -> Optional[Dict[str, Any]]:
    """Obtiene un diagrama completo por su ID incluyendo los nodos y flujos."""
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT id, name, description, module, data, created_at, updated_at
            FROM diagrams
            WHERE id = ?
        """, (diagram_id,))
        row = cursor.fetchone()
        if not row:
            return None
        res = dict(row)
        res["data"] = json.loads(res["data"])
        return res

def create_diagram(name: str, description: str, module: str, data: Dict[str, Any]) -> Dict[str, Any]:
    """Crea un nuevo diagrama en la base de datos."""
    now = datetime.now().isoformat()
    json_data = json.dumps(data)
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            INSERT INTO diagrams (name, description, module, data, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (name, description, module, json_data, now, now))
        conn.commit()
        diagram_id = cursor.lastrowid
        return {
            "id": diagram_id,
            "name": name,
            "description": description,
            "module": module,
            "data": data,
            "created_at": now,
            "updated_at": now
        }

def update_diagram(diagram_id: int, name: Optional[str] = None, description: Optional[str] = None, data: Optional[Dict[str, Any]] = None) -> Optional[Dict[str, Any]]:
    """Actualiza un diagrama existente."""
    existing = get_diagram_by_id(diagram_id)
    if not existing:
        return None

    now = datetime.now().isoformat()
    new_name = name if name is not None else existing["name"]
    new_desc = description if description is not None else existing["description"]
    new_data = json.dumps(data) if data is not None else json.dumps(existing["data"])

    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            UPDATE diagrams
            SET name = ?, description = ?, data = ?, updated_at = ?
            WHERE id = ?
        """, (new_name, new_desc, new_data, now, diagram_id))
        conn.commit()

    return get_diagram_by_id(diagram_id)

def delete_diagram(diagram_id: int) -> bool:
    """Elimina un diagrama por ID."""
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM diagrams WHERE id = ?", (diagram_id,))
        conn.commit()
        return cursor.rowcount > 0

# --- LECTURAS DE TELEMETRÍA (SENSOR DE NIVEL ESP32) ---
def save_level_reading(device_id: str, distance_cm: float, level_cm: float, percentage: float) -> Dict[str, Any]:
    """Guarda una nueva lectura de nivel enviada por el ESP32."""
    now = datetime.now().isoformat()
    with get_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            INSERT INTO level_readings (device_id, distance_cm, level_cm, percentage, created_at)
            VALUES (?, ?, ?, ?, ?)
        """, (device_id, distance_cm, level_cm, percentage, now))
        conn.commit()
        reading_id = cursor.lastrowid
        return {
            "id": reading_id,
            "device_id": device_id,
            "distance_cm": distance_cm,
            "level_cm": level_cm,
            "percentage": percentage,
            "created_at": now
        }

def get_latest_level_reading(device_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Obtiene la última lectura registrada por el ESP32."""
    with get_connection() as conn:
        cursor = conn.cursor()
        if device_id:
            cursor.execute("""
                SELECT id, device_id, distance_cm, level_cm, percentage, created_at
                FROM level_readings
                WHERE device_id = ?
                ORDER BY id DESC LIMIT 1
            """, (device_id,))
        else:
            cursor.execute("""
                SELECT id, device_id, distance_cm, level_cm, percentage, created_at
                FROM level_readings
                ORDER BY id DESC LIMIT 1
            """)
        row = cursor.fetchone()
        return dict(row) if row else None

def get_level_history(device_id: Optional[str] = None, limit: int = 50) -> List[Dict[str, Any]]:
    """Obtiene el historial reciente de lecturas para visualización en gráfico."""
    with get_connection() as conn:
        cursor = conn.cursor()
        if device_id:
            cursor.execute("""
                SELECT id, device_id, distance_cm, level_cm, percentage, created_at
                FROM level_readings
                WHERE device_id = ?
                ORDER BY id DESC LIMIT ?
            """, (device_id, limit))
        else:
            cursor.execute("""
                SELECT id, device_id, distance_cm, level_cm, percentage, created_at
                FROM level_readings
                ORDER BY id DESC LIMIT ?
            """, (limit,))
        rows = cursor.fetchall()
        # Retornamos en orden cronológico ascendente para el gráfico
        res = [dict(r) for r in rows]
        res.reverse()
        return res

# Inicializar tabla al importar
init_db()
