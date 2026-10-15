import numpy as np
from scipy.optimize import minimize
from water_models import WaterReconciliationRequest, NodeType
from typing import List, Dict, Any, Tuple

def _get_variable_mapping(request: WaterReconciliationRequest):
    """
    Mapea cada flujo a índices en el vector de estado (S: tonelaje de sólido seco, Cp: % sólidos)
    y genera valores iniciales inteligentes.
    """
    var_map = {}
    initial_guess = []
    bounds = []
    var_idx = 0

    # Heurística para flujos de sólido no medidos
    measured_solids = [f.data.get_solid_tonnage() for f in request.flows if not f.data.not_measured and f.data.get_solid_tonnage() > 0]
    avg_solid = float(np.mean(measured_solids)) if measured_solids else 30.0

    for flow in request.flows:
        s_val = flow.data.get_solid_tonnage()
        is_fixed = flow.data.is_solid_fixed()

        # Variable: Tonelaje Sólido Seco (S)
        var_map[f"s_{flow.id}"] = var_idx
        if flow.data.not_measured and s_val <= 0:
            s_val = avg_solid
        initial_guess.append(s_val)
        bounds.append((s_val, s_val) if is_fixed else (0.0, None))
        var_idx += 1

        # Variable: Porcentaje Sólidos (Cp)
        var_map[f"cp_{flow.id}"] = var_idx
        cp_val = flow.data.porcentaje_solidos
        # Si es no medido y está en cero o muy bajo, usar semilla de 30%
        if getattr(flow.data, 'porcentaje_solidos_not_measured', False) and cp_val <= 0 and not flow.data.porcentaje_solidos_fixed:
            cp_val = 30.0
        elif cp_val <= 0:
            cp_val = 30.0
        initial_guess.append(cp_val)
        
        # Límites para % Sólidos: entre 0.1 y 100% para evitar división por cero en pulpa
        if flow.data.porcentaje_solidos_fixed:
            bounds.append((flow.data.porcentaje_solidos, flow.data.porcentaje_solidos))
        else:
            bounds.append((0.1, 100.0))
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
    """Calcula los residuos de las restricciones de balance de sólidos secos y de pulpa total."""
    balance_nodes = _get_balance_nodes(request)
    residuals = []

    for node in balance_nodes:
        solids_mass_balance = 0.0
        pulp_mass_balance = 0.0

        for flow in request.flows:
            s_idx = var_map[f"s_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            S = y[s_idx]
            Cp = max(y[cp_idx], 0.01) # Protección contra cero

            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                # Balance de sólidos secos (Lineal)
                solids_mass_balance += factor * S
                # Balance de pulpa total: Q = S / (Cp / 100) = S * 100 / Cp
                pulp_mass_balance += factor * (S * 100.0 / Cp)

        residuals.append(solids_mass_balance)
        residuals.append(pulp_mass_balance)

    return np.array(residuals, dtype=float)

def _get_jacobian(y: np.ndarray, request: WaterReconciliationRequest, var_map: dict):
    """Calcula la matriz Jacobiana exacta de las restricciones de sólidos y pulpa."""
    balance_nodes = _get_balance_nodes(request)
    num_vars = len(y)
    num_constraints = len(balance_nodes) * 2
    J = np.zeros((num_constraints, num_vars))

    row = 0
    for node in balance_nodes:
        # Fila 1: Balance de sólidos secos (Lineal: d/dS = factor, d/dCp = 0)
        for flow in request.flows:
            s_idx = var_map[f"s_{flow.id}"]
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                J[row, s_idx] = factor
        row += 1

        # Fila 2: Balance de pulpa total: Q = S * 100 / Cp
        # d/dS = factor * 100 / Cp
        # d/dCp = - factor * 100 * S / Cp^2
        for flow in request.flows:
            s_idx = var_map[f"s_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            S = y[s_idx]
            Cp = max(y[cp_idx], 0.01)
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                J[row, s_idx] = factor * 100.0 / Cp
                J[row, cp_idx] = - factor * 100.0 * S / (Cp ** 2)
        row += 1

    return J

def reconcile_water_data(request: WaterReconciliationRequest):
    """
    Ejecuta la reconciliación del balance hídrico basada en tonelaje de sólido seco y % de sólidos
    usando SLSQP con Jacobiano analítico y matriz de covarianza posterior Gauss-Markov.
    """
    var_map, y_init, bounds = _get_variable_mapping(request)

    measured_indices, measured_values, variances = [], [], []

    for flow in request.flows:
        s_val = flow.data.get_solid_tonnage()
        is_fixed = flow.data.is_solid_fixed()

        # Sólidos Secos (S)
        if not flow.data.not_measured and not is_fixed:
            measured_indices.append(var_map[f"s_{flow.id}"])
            measured_values.append(s_val)
            # Protección contra varianza cero
            err = flow.data.get_solid_error()
            if err <= 0:
                err = max(s_val * 0.02 + 0.1, 0.1)
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
            S_mat = J @ V_in @ J.T
            V_out = V_in - V_in @ J.T @ np.linalg.pinv(S_mat) @ J @ V_in
            std_devs = np.sqrt(np.maximum(np.diag(V_out), 0.0))
        else:
            std_devs = np.sqrt(np.diag(V_in))
    except Exception:
        std_devs = np.zeros_like(y_hat)

    # Forzar error cero en variables fijas
    for flow in request.flows:
        if flow.data.is_solid_fixed():
            std_devs[var_map[f"s_{flow.id}"]] = 0.0
        if flow.data.porcentaje_solidos_fixed:
            std_devs[var_map[f"cp_{flow.id}"]] = 0.0

    return _build_response(request, y_hat, std_devs, var_map, success=res.success)

def _build_response(request: WaterReconciliationRequest, y_hat: np.ndarray, std_devs: np.ndarray, var_map: dict, success=True):
    reconciled_flows = []
    for flow in request.flows:
        s_idx = var_map[f"s_{flow.id}"]
        cp_idx = var_map[f"cp_{flow.id}"]
        s_val = float(y_hat[s_idx])
        cp_val = float(y_hat[cp_idx])
        s_err = float(std_devs[s_idx])
        cp_err = float(std_devs[cp_idx])

        # Pulpa calculada: Q = S / (Cp / 100)
        q_calc = float(s_val * 100.0 / cp_val) if cp_val > 0 else 0.0
        # Agua calculada: W = Q - S
        w_calc = float(q_calc - s_val)

        # Propagación de error para pulpa calculada: dQ = sqrt( (100/Cp * dS)^2 + (-100*S/Cp^2 * dCp)^2 )
        q_err = float(np.sqrt((100.0 / cp_val * s_err) ** 2 + (100.0 * s_val / (cp_val ** 2) * cp_err) ** 2)) if cp_val > 0 else 0.0

        reconciled_flows.append({
            "id": flow.id,
            "corrected_data": {
                "tonelaje_solido_corregido": s_val,
                "tonelaje_solido_error": s_err,
                "porcentaje_solidos_corregido": cp_val,
                "porcentaje_solidos_error": cp_err,
                "tonelaje_pulpa_calculado": q_calc,
                "tonelaje_agua_calculado": w_calc,
                # Compatibilidad hacia atrás:
                "tonelaje_solido_calculado": s_val,
                "tonelaje_pulpa_corregido": q_calc,
                "tonelaje_pulpa_error": q_err
            }
        })

    # Residuos nodales finales
    balance_nodes = _get_balance_nodes(request)
    residuals = []
    for node in balance_nodes:
        in_solids = sum(y_hat[var_map[f"s_{f.id}"]] for f in request.flows if f.target == node.id)
        out_solids = sum(y_hat[var_map[f"s_{f.id}"]] for f in request.flows if f.source == node.id)
        
        in_pulp = sum(y_hat[var_map[f"s_{f.id}"]] * 100.0 / max(y_hat[var_map[f"cp_{f.id}"]], 0.01) for f in request.flows if f.target == node.id)
        out_pulp = sum(y_hat[var_map[f"s_{f.id}"]] * 100.0 / max(y_hat[var_map[f"cp_{f.id}"]], 0.01) for f in request.flows if f.source == node.id)

        in_water = in_pulp - in_solids
        out_water = out_pulp - out_solids

        residuals.append({
            "node": node.label,
            "residual_solidos": float(in_solids - out_solids),
            "residual_pulpa": float(in_pulp - out_pulp),
            "residual_agua": float(in_water - out_water)
        })

    report = {
        "success": bool(success),
        "methodology": [
            "Optimización no lineal SLSQP basada en conservación de masa de sólidos secos y balance hídrico.",
            "Variables independientes de medición y ajuste: Tonelaje de Sólidos Secos (S) y Porcentaje de Sólidos (%S / Cp).",
            "Cálculo estequiométrico exacto de caudal de pulpa y agua a partir de sólidos y concentración.",
            "Estimación de covarianza e incertidumbre posterior mediante proyección analítica de Gauss-Markov."
        ],
        "equations": _format_physical_equations(request),
        "node_residuals": residuals
    }
    return reconciled_flows, report

def _format_physical_equations(request: WaterReconciliationRequest) -> List[str]:
    """Genera ecuaciones algebraicas legibles para previsualización con sólido seco y % sólidos."""
    nodes = _get_balance_nodes(request)
    eqs = []
    for node in nodes:
        in_s = [f"S_{f.id}" for f in request.flows if f.target == node.id]
        out_s = [f"S_{f.id}" for f in request.flows if f.source == node.id]
        in_s_str = " + ".join(in_s) if in_s else "0"
        out_s_str = " - ".join(out_s) if out_s else "0"
        eqs.append(f"[{node.label}] Balance de Sólidos Secos: {in_s_str} - {out_s_str} = 0")

        in_p = [f"(100/Cp_{f.id})*S_{f.id}" for f in request.flows if f.target == node.id]
        out_p = [f"(100/Cp_{f.id})*S_{f.id}" for f in request.flows if f.source == node.id]
        in_p_str = " + ".join(in_p) if in_p else "0"
        out_p_str = " - ".join(out_p) if out_p else "0"
        eqs.append(f"[{node.label}] Balance de Pulpa Total: {in_p_str} - {out_p_str} = 0")
    return eqs
