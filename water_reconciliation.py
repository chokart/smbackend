import numpy as np
from scipy.optimize import minimize
from water_models import WaterReconciliationRequest, NodeType
from typing import List, Dict, Any, Tuple

def _get_variable_mapping(request: WaterReconciliationRequest):
    """
    Mapea cada flujo a índices en el vector de estado y genera valores iniciales inteligentes.
    Evita inicializaciones en cero para que las derivadas parciales del Jacobiano no colapsen.
    """
    var_map = {}
    initial_guess = []
    bounds = []
    var_idx = 0

    # Heurística para flujos no medidos
    measured_pulps = [f.data.tonelaje_pulpa for f in request.flows if not f.data.not_measured and f.data.tonelaje_pulpa > 0]
    avg_pulp = float(np.mean(measured_pulps)) if measured_pulps else 100.0

    for flow in request.flows:
        # Variable: Tonelaje Pulpa (Q)
        var_map[f"q_{flow.id}"] = var_idx
        q_val = flow.data.tonelaje_pulpa
        if flow.data.not_measured and q_val <= 0:
            q_val = avg_pulp
        initial_guess.append(q_val)
        bounds.append((flow.data.tonelaje_pulpa, flow.data.tonelaje_pulpa) if flow.data.tonelaje_pulpa_fixed else (0.0, None))
        var_idx += 1

        # Variable: Porcentaje Sólidos (Cp)
        var_map[f"cp_{flow.id}"] = var_idx
        cp_val = flow.data.porcentaje_solidos
        # Si es no medido y está en cero, usar semilla de 30% para garantizar gradiente activo
        if getattr(flow.data, 'porcentaje_solidos_not_measured', False) and cp_val <= 0 and not flow.data.porcentaje_solidos_fixed:
            cp_val = 30.0
        initial_guess.append(cp_val)
        
        # Límites para % Sólidos: entre 0 y 100%
        if flow.data.porcentaje_solidos_fixed:
            bounds.append((flow.data.porcentaje_solidos, flow.data.porcentaje_solidos))
        else:
            bounds.append((0.0, 100.0))
        var_idx += 1

    return var_map, np.array(initial_guess, dtype=float), bounds

def _get_balance_nodes(request: WaterReconciliationRequest):
    """Identifica nodos que deben cumplir balance (tienen entradas y salidas, o son de proceso)."""
    balance_nodes = []
    for node in request.nodes:
        has_in = any(f.target == node.id for f in request.flows)
        has_out = any(f.source == node.id for f in request.flows)
        if (has_in and has_out) or node.node_type == NodeType.PROCESO:
            balance_nodes.append(node)
    return balance_nodes

def _get_constraints(y: np.ndarray, request: WaterReconciliationRequest, var_map: dict):
    """Calcula los residuos de las restricciones de balance de pulpa y de sólidos secos."""
    balance_nodes = _get_balance_nodes(request)
    residuals = []

    for node in balance_nodes:
        total_mass_balance = 0.0
        solids_mass_balance = 0.0

        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            Q = y[q_idx]
            Cp = y[cp_idx]

            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                total_mass_balance += factor * Q
                solids_mass_balance += factor * (Q * Cp / 100.0)

        residuals.append(total_mass_balance)
        residuals.append(solids_mass_balance)

    return np.array(residuals, dtype=float)

def _get_jacobian(y: np.ndarray, request: WaterReconciliationRequest, var_map: dict):
    """Calcula la matriz Jacobiana exacta de las restricciones no lineales."""
    balance_nodes = _get_balance_nodes(request)
    num_vars = len(y)
    num_constraints = len(balance_nodes) * 2
    J = np.zeros((num_constraints, num_vars))

    row = 0
    for node in balance_nodes:
        # Fila 1: Balance de masa total de pulpa (Lineal)
        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                J[row, q_idx] = factor
        row += 1

        # Fila 2: Balance de sólidos secos (No lineal: d/dQ = factor * Cp/100, d/dCp = factor * Q/100)
        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            Q = y[q_idx]
            Cp = y[cp_idx]
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                J[row, q_idx] = factor * Cp / 100.0
                J[row, cp_idx] = factor * Q / 100.0
        row += 1

    return J

def reconcile_water_data(request: WaterReconciliationRequest):
    """
    Ejecuta la reconciliación del balance hídrico y de pulpa usando SLSQP
    con Jacobiano analítico y matriz de covarianza posterior Gauss-Markov.
    """
    var_map, y_init, bounds = _get_variable_mapping(request)

    measured_indices, measured_values, variances = [], [], []
    for flow in request.flows:
        # Pulpa (Q)
        if not flow.data.not_measured and not flow.data.tonelaje_pulpa_fixed:
            measured_indices.append(var_map[f"q_{flow.id}"])
            measured_values.append(flow.data.tonelaje_pulpa)
            # Protección contra varianza cero
            err = flow.data.tonelaje_pulpa_error if flow.data.tonelaje_pulpa_error > 0 else max(flow.data.tonelaje_pulpa * 0.02 + 0.1, 0.1)
            variances.append(err ** 2)

        # % Sólidos (Cp)
        cp_not_measured = getattr(flow.data, 'porcentaje_solidos_not_measured', False)
        if not cp_not_measured and not flow.data.porcentaje_solidos_fixed:
            measured_indices.append(var_map[f"cp_{flow.id}"])
            measured_values.append(flow.data.porcentaje_solidos)
            # Protección contra varianza cero
            err = flow.data.porcentaje_solidos_error if flow.data.porcentaje_solidos_error > 0 else max(flow.data.porcentaje_solidos * 0.01 + 0.2, 0.2)
            variances.append(err ** 2)

    if not measured_indices:
        return _build_response(request, y_init, np.zeros_like(y_init), var_map, True)

    W = np.diag(1.0 / np.array(variances, dtype=float))
    m_vals = np.array(measured_values, dtype=float)

    def objective(y):
        diff = y[measured_indices] - m_vals
        return 0.5 * diff.T @ W @ diff

    def objective_grad(y):
        grad = np.zeros_like(y)
        diff = y[measured_indices] - m_vals
        grad[measured_indices] = W @ diff
        return grad

    constraints = {
        'type': 'eq',
        'fun': _get_constraints,
        'jac': _get_jacobian,
        'args': (request, var_map)
    }

    res = minimize(
        objective, y_init, jac=objective_grad, method='SLSQP',
        bounds=bounds, constraints=constraints, tol=1e-8,
        options={'maxiter': 250, 'ftol': 1e-9}
    )

    y_hat = res.x

    # --- INCERTIDUMBRE POSTERIOR (GAUSS-MARKOV) ---
    V_in = np.zeros((len(y_hat), len(y_hat)))
    for i, idx in enumerate(measured_indices):
        V_in[idx, idx] = variances[i]

    J = _get_jacobian(y_hat, request, var_map)
    try:
        if J.size > 0:
            S = J @ V_in @ J.T
            V_out = V_in - V_in @ J.T @ np.linalg.pinv(S) @ J @ V_in
            std_devs = np.sqrt(np.maximum(np.diag(V_out), 0.0))
        else:
            std_devs = np.sqrt(np.diag(V_in))
    except Exception:
        std_devs = np.zeros_like(y_hat)

    # Forzar error cero en variables fijas
    for flow in request.flows:
        if flow.data.tonelaje_pulpa_fixed:
            std_devs[var_map[f"q_{flow.id}"]] = 0.0
        if flow.data.porcentaje_solidos_fixed:
            std_devs[var_map[f"cp_{flow.id}"]] = 0.0

    return _build_response(request, y_hat, std_devs, var_map, success=res.success)

def _build_response(request: WaterReconciliationRequest, y_hat: np.ndarray, std_devs: np.ndarray, var_map: dict, success=True):
    reconciled_flows = []
    for flow in request.flows:
        q_idx = var_map[f"q_{flow.id}"]
        cp_idx = var_map[f"cp_{flow.id}"]
        q_val = float(y_hat[q_idx])
        cp_val = float(y_hat[cp_idx])

        reconciled_flows.append({
            "id": flow.id,
            "corrected_data": {
                "tonelaje_pulpa_corregido": q_val,
                "tonelaje_pulpa_error": float(std_devs[q_idx]),
                "porcentaje_solidos_corregido": cp_val,
                "porcentaje_solidos_error": float(std_devs[cp_idx]),
                "tonelaje_solido_calculado": float(q_val * cp_val / 100.0),
                "tonelaje_agua_calculado": float(q_val * (1.0 - cp_val / 100.0))
            }
        })

    # Residuos nodales finales
    balance_nodes = _get_balance_nodes(request)
    residuals = []
    for node in balance_nodes:
        in_pulp = sum(y_hat[var_map[f"q_{f.id}"]] for f in request.flows if f.target == node.id)
        out_pulp = sum(y_hat[var_map[f"q_{f.id}"]] for f in request.flows if f.source == node.id)
        in_solids = sum(y_hat[var_map[f"q_{f.id}"]] * y_hat[var_map[f"cp_{f.id}"]] / 100.0 for f in request.flows if f.target == node.id)
        out_solids = sum(y_hat[var_map[f"q_{f.id}"]] * y_hat[var_map[f"cp_{f.id}"]] / 100.0 for f in request.flows if f.source == node.id)
        residuals.append({
            "node": node.label,
            "residual_pulpa": float(in_pulp - out_pulp),
            "residual_solidos": float(in_solids - out_solids),
            "residual_agua": float((in_pulp - in_solids) - (out_pulp - out_solids))
        })

    report = {
        "success": bool(success),
        "methodology": [
            "Optimización no lineal SLSQP con gradiente y Jacobiano analítico de alta precisión.",
            "Detección topológica automática de nodos de conservación de masa de agua y pulpa.",
            "Estimación de covarianza e incertidumbre posterior mediante proyección de Gauss-Markov.",
            "Cálculo estequiométrico exacto de tonelaje de agua y sólidos por flujo reconciliado."
        ],
        "equations": _format_physical_equations(request),
        "node_residuals": residuals
    }
    return reconciled_flows, report

def _format_physical_equations(request: WaterReconciliationRequest) -> List[str]:
    """Genera ecuaciones algebraicas legibles para previsualización."""
    nodes = _get_balance_nodes(request)
    eqs = []
    for node in nodes:
        in_f = [f"Q_{f.id}" for f in request.flows if f.target == node.id]
        out_f = [f"Q_{f.id}" for f in request.flows if f.source == node.id]
        in_str = " + ".join(in_f) if in_f else "0"
        out_str = " - ".join(out_f) if out_f else "0"
        eqs.append(f"[{node.label}] Balance de Pulpa: {in_str} - {out_str} = 0")

        in_s = [f"(Cp_{f.id}/100)*Q_{f.id}" for f in request.flows if f.target == node.id]
        out_s = [f"(Cp_{f.id}/100)*Q_{f.id}" for f in request.flows if f.source == node.id]
        in_s_str = " + ".join(in_s) if in_s else "0"
        out_s_str = " - ".join(out_s) if out_s else "0"
        eqs.append(f"[{node.label}] Balance de Sólidos Secos: {in_s_str} - {out_s_str} = 0")
    return eqs
