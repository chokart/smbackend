import numpy as np
from scipy.optimize import minimize
from water_models import WaterReconciliationRequest, NodeType
from typing import List, Dict, Any, Tuple

def _get_variable_mapping(request: WaterReconciliationRequest):
    """
    Mapea cada flujo a índices en el vector de estado y genera valores iniciales inteligentes:
    - Flujo de Pulpa: S (sólido seco) y Cp (% sólidos)
    - Flujo de Solo Agua: W (caudal de agua pura)
    """
    var_map = {}
    initial_guess = []
    bounds = []
    var_idx = 0

    # Heurística para estimar flujos de sólido no medidos
    measured_solids = [
        f.data.get_solid_tonnage() for f in request.flows 
        if not f.data.is_water_only and not f.data.not_measured and f.data.get_solid_tonnage() > 0
    ]
    avg_solid = float(np.mean(measured_solids)) if measured_solids else 50.0

    # Heurística para flujos de agua pura no medidos
    measured_waters = [
        f.data.get_water_tonnage() for f in request.flows 
        if f.data.is_water_only and not f.data.not_measured and f.data.get_water_tonnage() > 0
    ]
    avg_water = float(np.mean(measured_waters)) if measured_waters else 100.0

    for flow in request.flows:
        if flow.data.is_water_only:
            # Variable única para flujos de solo agua: W (Caudal de Agua en tph)
            var_map[f"w_{flow.id}"] = var_idx
            w_val = flow.data.get_water_tonnage()
            if flow.data.not_measured and w_val <= 0:
                w_val = avg_water
            initial_guess.append(w_val)
            bounds.append((w_val, w_val) if flow.data.is_water_fixed() else (0.0, None))
            var_idx += 1
        else:
            # Flujo estándar de pulpa mineral
            s_val = flow.data.get_solid_tonnage()
            is_fixed = flow.data.is_solid_fixed()

            # Variable 1: Tonelaje Sólido Seco (S)
            var_map[f"s_{flow.id}"] = var_idx
            if flow.data.not_measured and s_val <= 0:
                s_val = avg_solid
            initial_guess.append(s_val)
            bounds.append((s_val, s_val) if is_fixed else (0.0, None))
            var_idx += 1

            # Variable 2: Porcentaje Sólidos (Cp)
            var_map[f"cp_{flow.id}"] = var_idx
            cp_val = flow.data.porcentaje_solidos
            if getattr(flow.data, 'porcentaje_solidos_not_measured', False) and cp_val <= 0 and not flow.data.porcentaje_solidos_fixed:
                cp_val = 30.0
            elif cp_val <= 0:
                cp_val = 30.0
            initial_guess.append(cp_val)
            
            # Límites para % Sólidos: entre 0.1 y 100%
            if flow.data.porcentaje_solidos_fixed:
                bounds.append((flow.data.porcentaje_solidos, flow.data.porcentaje_solidos))
            else:
                bounds.append((0.1, 100.0))
            var_idx += 1

    return var_map, np.array(initial_guess, dtype=float), bounds

def _get_balance_nodes(request: WaterReconciliationRequest):
    """
    Identifica nodos que deben cumplir balance de conservación de masa.
    Un nodo solo puede y debe balancear si tiene tanto flujos de entrada como de salida.
    Los nodos con solo salidas son fuentes (alimentos/agua) y con solo entradas son sumideros (productos).
    """
    balance_nodes = []
    for node in request.nodes:
        has_in = any(f.target == node.id for f in request.flows)
        has_out = any(f.source == node.id for f in request.flows)
        if has_in and has_out:
            balance_nodes.append(node)
    return balance_nodes

def _get_constraints(y: np.ndarray, request: WaterReconciliationRequest, var_map: dict):
    """
    Calcula los residuos de las dos leyes de conservación fundamentales e independientes:
    1. Balance de Masa de Sólidos Secos (tph)
    2. Balance de Masa de Agua (tph)
    (El balance de pulpa total Q = S + W se satisface automáticamente).
    """
    balance_nodes = _get_balance_nodes(request)
    residuals = []

    for node in balance_nodes:
        solids_mass_balance = 0.0
        water_mass_balance = 0.0

        for flow in request.flows:
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor == 0:
                continue

            if flow.data.is_water_only:
                # Flujo de solo agua: S = 0, W = y[w_idx]
                w_idx = var_map[f"w_{flow.id}"]
                W = y[w_idx]
                water_mass_balance += factor * W
            else:
                # Flujo de pulpa: S = y[s_idx], W = S * (100 / Cp - 1)
                s_idx = var_map[f"s_{flow.id}"]
                cp_idx = var_map[f"cp_{flow.id}"]
                S = y[s_idx]
                Cp = max(y[cp_idx], 0.01)

                solids_mass_balance += factor * S
                water_mass_balance += factor * (S * (100.0 / Cp - 1.0))

        residuals.append(solids_mass_balance)
        residuals.append(water_mass_balance)

    return np.array(residuals, dtype=float)

def _get_jacobian(y: np.ndarray, request: WaterReconciliationRequest, var_map: dict):
    """Calcula la matriz Jacobiana exacta de las restricciones de balance de sólidos secos y de agua."""
    balance_nodes = _get_balance_nodes(request)
    num_vars = len(y)
    num_constraints = len(balance_nodes) * 2
    J = np.zeros((num_constraints, num_vars))

    row = 0
    for node in balance_nodes:
        # Fila 1: Balance de sólidos secos
        for flow in request.flows:
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                if not flow.data.is_water_only:
                    s_idx = var_map[f"s_{flow.id}"]
                    J[row, s_idx] = factor
        row += 1

        # Fila 2: Balance de agua
        for flow in request.flows:
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                if flow.data.is_water_only:
                    w_idx = var_map[f"w_{flow.id}"]
                    J[row, w_idx] = factor * 1.0
                else:
                    s_idx = var_map[f"s_{flow.id}"]
                    cp_idx = var_map[f"cp_{flow.id}"]
                    S = y[s_idx]
                    Cp = max(y[cp_idx], 0.01)
                    # dW/dS = factor * (100 / Cp - 1)
                    J[row, s_idx] = factor * (100.0 / Cp - 1.0)
                    # dW/dCp = factor * (- 100 * S / Cp^2)
                    J[row, cp_idx] = - factor * 100.0 * S / (Cp ** 2)
        row += 1

    return J

def reconcile_water_data(request: WaterReconciliationRequest):
    """
    Ejecuta la reconciliación del balance hídrico basada en tonelaje de sólido seco, % de sólidos y agua pura
    usando SLSQP con Jacobiano analítico y matriz de covarianza posterior Gauss-Markov.
    """
    var_map, y_init, bounds = _get_variable_mapping(request)

    measured_indices, measured_values, variances = [], [], []

    for flow in request.flows:
        if flow.data.is_water_only:
            # Flujo de Agua Pura (W)
            if not flow.data.not_measured and not flow.data.is_water_fixed():
                measured_indices.append(var_map[f"w_{flow.id}"])
                w_val = flow.data.get_water_tonnage()
                measured_values.append(w_val)
                err = flow.data.get_water_error()
                variances.append(err ** 2)
        else:
            # Sólidos Secos (S)
            s_val = flow.data.get_solid_tonnage()
            is_fixed = flow.data.is_solid_fixed()
            if not flow.data.not_measured and not is_fixed:
                measured_indices.append(var_map[f"s_{flow.id}"])
                measured_values.append(s_val)
                err = flow.data.get_solid_error()
                if err <= 0:
                    err = max(s_val * 0.02 + 0.1, 0.1)
                variances.append(err ** 2)

            # % Sólidos (Cp)
            cp_not_measured = getattr(flow.data, 'porcentaje_solidos_not_measured', False)
            if not cp_not_measured and not flow.data.porcentaje_solidos_fixed:
                measured_indices.append(var_map[f"cp_{flow.id}"])
                measured_values.append(flow.data.porcentaje_solidos)
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
        bounds=bounds, constraints=constraints, tol=1e-7,
        options={'maxiter': 250, 'ftol': 1e-8}
    )

    y_hat = res.x

    # Comprobación de cumplimiento físico de restricciones
    constraints_residuals = _get_constraints(y_hat, request, var_map)
    max_residual = float(np.max(np.abs(constraints_residuals))) if len(constraints_residuals) > 0 else 0.0
    is_success = bool(res.success or (max_residual < 0.05))

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
        if flow.data.is_water_only:
            if flow.data.is_water_fixed():
                std_devs[var_map[f"w_{flow.id}"]] = 0.0
        else:
            if flow.data.is_solid_fixed():
                std_devs[var_map[f"s_{flow.id}"]] = 0.0
            if flow.data.porcentaje_solidos_fixed:
                std_devs[var_map[f"cp_{flow.id}"]] = 0.0

    return _build_response(request, y_hat, std_devs, var_map, success=is_success)

def _build_response(request: WaterReconciliationRequest, y_hat: np.ndarray, std_devs: np.ndarray, var_map: dict, success=True):
    reconciled_flows = []
    for flow in request.flows:
        if flow.data.is_water_only:
            w_idx = var_map[f"w_{flow.id}"]
            w_val = float(y_hat[w_idx])
            w_err = float(std_devs[w_idx])
            s_val = 0.0
            s_err = 0.0
            cp_val = 0.0
            cp_err = 0.0
            q_calc = w_val
            q_err = w_err
        else:
            s_idx = var_map[f"s_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            s_val = float(y_hat[s_idx])
            cp_val = float(y_hat[cp_idx])
            s_err = float(std_devs[s_idx])
            cp_err = float(std_devs[cp_idx])

            # Pulpa calculada: Q = S / (Cp / 100)
            q_calc = float(s_val * 100.0 / cp_val) if cp_val > 0 else 0.0
            # Agua calculada: W = Q - S
            w_val = float(q_calc - s_val)
            w_err = float(np.sqrt((s_err * (100.0 / cp_val - 1.0)) ** 2 + (s_val * 100.0 / (cp_val ** 2) * cp_err) ** 2)) if cp_val > 0 else 0.0
            q_err = float(np.sqrt((100.0 / cp_val * s_err) ** 2 + (100.0 * s_val / (cp_val ** 2) * cp_err) ** 2)) if cp_val > 0 else 0.0

        reconciled_flows.append({
            "id": flow.id,
            "corrected_data": {
                "tonelaje_solido_corregido": s_val,
                "tonelaje_solido_error": s_err,
                "porcentaje_solidos_corregido": cp_val,
                "porcentaje_solidos_error": cp_err,
                "tonelaje_pulpa_calculado": q_calc,
                "tonelaje_agua_calculado": w_val,
                "tonelaje_agua_error": w_err,
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
        in_solids, out_solids = 0.0, 0.0
        in_water, out_water = 0.0, 0.0
        in_pulp, out_pulp = 0.0, 0.0

        for f in request.flows:
            factor = 1.0 if f.target == node.id else (-1.0 if f.source == node.id else 0.0)
            if factor == 0:
                continue

            if f.data.is_water_only:
                w_f = float(y_hat[var_map[f"w_{f.id}"]])
                s_f = 0.0
            else:
                s_f = float(y_hat[var_map[f"s_{f.id}"]])
                cp_f = max(float(y_hat[var_map[f"cp_{f.id}"]]), 0.01)
                w_f = float(s_f * (100.0 / cp_f - 1.0))

            q_f = s_f + w_f

            if factor > 0:
                in_solids += s_f
                in_water += w_f
                in_pulp += q_f
            else:
                out_solids += s_f
                out_water += w_f
                out_pulp += q_f

        residuals.append({
            "node": node.label,
            "residual_solidos": float(in_solids - out_solids),
            "residual_pulpa": float(in_pulp - out_pulp),
            "residual_agua": float(in_water - out_water)
        })

    report = {
        "success": bool(success),
        "methodology": [
            "Optimización no lineal SLSQP con restricciones desacopladas de conservación de sólidos secos y balance hídrico.",
            "Soporte integral para corrientes bifásicas de pulpa (Sólido y %S) y corrientes monofásicas de agua pura (Agua de Dilución).",
            "Cálculo estequiométrico exacto de caudal de pulpa y agua sin singularidades en cero sólidos.",
            "Estimación de covarianza e incertidumbre posterior mediante proyección analítica de Gauss-Markov."
        ],
        "equations": _format_physical_equations(request),
        "node_residuals": residuals
    }
    return reconciled_flows, report

def _format_physical_equations(request: WaterReconciliationRequest) -> List[str]:
    """Genera ecuaciones algebraicas legibles para previsualización con sólido seco y agua."""
    nodes = _get_balance_nodes(request)
    eqs = []
    for node in nodes:
        # Balance de Sólidos Secos
        in_s = [f"S_{f.id}" for f in request.flows if f.target == node.id and not f.data.is_water_only]
        out_s = [f"S_{f.id}" for f in request.flows if f.source == node.id and not f.data.is_water_only]
        in_s_str = " + ".join(in_s) if in_s else "0"
        out_s_str = " - ".join(out_s) if out_s else "0"
        eqs.append(f"[{node.label}] Balance de Sólidos Secos: {in_s_str} - {out_s_str} = 0")

        # Balance de Agua
        in_w, out_w = [], []
        for f in request.flows:
            if f.data.is_water_only:
                term = f"W_{f.id}"
            else:
                term = f"W(S_{f.id}, Cp_{f.id})"
            if f.target == node.id:
                in_w.append(term)
            elif f.source == node.id:
                out_w.append(term)

        in_w_str = " + ".join(in_w) if in_w else "0"
        out_w_str = " - ".join(out_w) if out_w else "0"
        eqs.append(f"[{node.label}] Balance de Agua: {in_w_str} - {out_w_str} = 0")
    return eqs
