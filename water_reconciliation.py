import numpy as np
from scipy.optimize import minimize
from water_models import WaterReconciliationRequest, NodeType, WaterFlow
from typing import List

def _get_water_balance_system(request: WaterReconciliationRequest, y_k: np.ndarray = None):
    """
    Construye el sistema de ecuaciones para el balance de agua/sólidos.
    Variables: [FlujoPulpa_1, %Solidos_1, FlujoPulpa_2, %Solidos_2, ...]
    """
    var_map, initial_guess, bounds, fixed_vars_indices = {}, [], [], set()
    var_idx = 0
    
    # 1. Definir variables y límites
    for flow in request.flows:
        # Variable 1: Tonelaje Total de Pulpa (Q)
        var_map[f"q_{flow.id}"] = var_idx
        initial_guess.append(flow.data.tonelaje_pulpa)
        # Límite: No negativo. Si está fijo, se fija el valor.
        bounds.append((flow.data.tonelaje_pulpa, flow.data.tonelaje_pulpa) if flow.data.tonelaje_pulpa_fixed else (0, None))
        if flow.data.tonelaje_pulpa_fixed: fixed_vars_indices.add(var_idx)
        var_idx += 1
        
        # Variable 2: Porcentaje de Sólidos (%S)
        var_map[f"cp_{flow.id}"] = var_idx
        initial_guess.append(flow.data.porcentaje_solidos)
        # Límite: 0 a 100%. Si está fijo, se fija el valor.
        bounds.append((flow.data.porcentaje_solidos, flow.data.porcentaje_solidos) if flow.data.porcentaje_solidos_fixed else (0, 100))
        if flow.data.porcentaje_solidos_fixed: fixed_vars_indices.add(var_idx)
        var_idx += 1

    if y_k is None:
        y_k = np.array(initial_guess)

    # 2. Definir Restricciones (Balance de Masa)
    process_nodes = [n for n in request.nodes if n.node_type == NodeType.PROCESO]
    num_constraints = len(process_nodes) * 2  # 1 Total + 1 Sólidos por nodo
    
    A = np.zeros((num_constraints, len(y_k)))
    b = np.zeros(num_constraints)
    constraint_idx = 0

    for node in process_nodes:
        # A) Balance de Masa Total (Entradas = Salidas) -> Sum(Q_in) - Sum(Q_out) = 0
        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            
            factor = 0
            if flow.target == node.id: factor = 1  # Entrada
            elif flow.source == node.id: factor = -1 # Salida
            
            if factor != 0:
                A[constraint_idx, q_idx] += factor
        
        b[constraint_idx] = 0 # Siempre 0 para estado estacionario
        constraint_idx += 1

        # B) Balance de Sólidos (Entradas = Salidas) -> Sum(Q_in * Cp_in/100) - Sum(Q_out * Cp_out/100) = 0
        # Linealización Taylor: f(x) ~ f(xk) + grad(f) * (x - xk) = 0
        # A * x = grad(f) * xk - f(xk)
        
        current_imbalance_f_xk = 0
        grad_dot_xk = 0

        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            
            Q_k = y_k[q_idx] # Valor iteración anterior
            Cp_k = y_k[cp_idx] # Valor iteración anterior
            
            factor = 0
            if flow.target == node.id: factor = 1
            elif flow.source == node.id: factor = -1
            
            if factor != 0:
                # Derivadas parciales (Gradiente en A)
                dq = factor * Cp_k / 100.0
                dcp = factor * Q_k / 100.0
                
                A[constraint_idx, q_idx] += dq
                A[constraint_idx, cp_idx] += dcp
                
                # Cálculos para b
                term_val = factor * Q_k * Cp_k / 100.0
                current_imbalance_f_xk += term_val
                grad_dot_xk += (dq * Q_k + dcp * Cp_k)
        
        # b = grad(f) * xk - f(xk)
        # Nota matemática: grad_dot_xk para bilineal x*y es x*y + y*x = 2xy.
        # f(xk) es xy. Entonces b = 2xy - xy = xy.
        # Por lo tanto, b es simplemente el valor acumulado actual, pero calculado explícitamente.
        b[constraint_idx] = grad_dot_xk - current_imbalance_f_xk
        constraint_idx += 1
            
    return A, b, var_map, initial_guess, bounds, fixed_vars_indices

def format_water_equations(A: np.ndarray, b: np.ndarray, var_map: dict) -> List[str]:
    """Genera strings legibles de las ecuaciones."""
    equations = []
    idx_to_name = {v: k for k, v in var_map.items()}

    for i, row in enumerate(A):
        equation_str = ""
        for j, coeff in enumerate(row):
            if abs(coeff) > 1e-9:
                var_name = idx_to_name.get(j, f"var_{j}")
                sign = "+" if coeff >= 0 else "-"
                equation_str += f" {sign} {abs(coeff):.4f} * {var_name}"
        
        if equation_str.strip().startswith("+"):
            equation_str = equation_str.strip()[1:].strip()
        
        equation_str += f" = {b[i]:.4f}"
        equations.append(equation_str.strip())
        
    return equations

def reconcile_water_data(request: WaterReconciliationRequest):
    MAX_ITERATIONS = 30
    TOLERANCE = 1e-7

    A, b, var_map, initial_guess, bounds, fixed_vars_indices = _get_water_balance_system(request)
    y_k = np.array(initial_guess)

    # Construir matrices de varianza y pesos
    measured_indices, measured_values, variances = [], [], []
    
    for flow in request.flows:
        # Flujo Pulpa
        if not flow.data.not_measured and not flow.data.tonelaje_pulpa_fixed:
            idx = var_map[f"q_{flow.id}"]
            measured_indices.append(idx)
            measured_values.append(flow.data.tonelaje_pulpa)
            # Evitar varianza 0
            err = flow.data.tonelaje_pulpa_error if flow.data.tonelaje_pulpa_error > 0 else 1e-6
            variances.append(err**2)
            
        # % Sólidos
        if not flow.data.porcentaje_solidos_not_measured and not flow.data.porcentaje_solidos_fixed:
            idx = var_map[f"cp_{flow.id}"]
            measured_indices.append(idx)
            measured_values.append(flow.data.porcentaje_solidos)
            err = flow.data.porcentaje_solidos_error if flow.data.porcentaje_solidos_error > 0 else 1e-6
            variances.append(err**2)

    x_measured = np.array(measured_values)
    V_inv = np.linalg.inv(np.diag(variances)) if variances else np.eye(len(x_measured))

    H = np.zeros((len(x_measured), len(y_k)))
    for i, idx in enumerate(measured_indices):
        H[i, idx] = 1

    A_final, b_final = None, None
    
    # Bucle de Optimización (Iterativo para linealizar restricciones bilineales)
    for iteration in range(MAX_ITERATIONS):
        # Actualizar linealización con los valores actuales (y_k)
        A, b, _, _, _, _ = _get_water_balance_system(request, y_k=y_k)

        def objective(y):
            if not measured_values: return 0.0
            error_vector = H @ y - x_measured
            return error_vector.T @ V_inv @ error_vector

        # Restricción: A*y = b (aproximadamente, en realidad b se ajusta en el bucle)
        # Nota: La lógica de b en tu código original parece requerir revisión si b no es 0.
        # Asumiremos la lógica estándar SQP: A * y = b_linearized
        constraints = [{'type': 'eq', 'fun': lambda y: A @ y - b}]
        
        result = minimize(objective, y_k, bounds=bounds, constraints=constraints, method='SLSQP')

        if not result.success:
            print(f"Warning: Optimizacion no convergió en iter {iteration}: {result.message}")
        
        y_new = result.x
        A_final, b_final = A, b
        
        if np.linalg.norm(y_new - y_k) < TOLERANCE:
            y_k = y_new
            break
        y_k = y_new

    y_hat = y_k
    
    # Calcular errores ajustados (simplificado: asumiendo independientes)
    # En un sistema real, usaríamos la matriz de covarianza posterior (J^T V^-1 J)^-1
    std_devs_y = np.zeros_like(y_hat) 

    # Construir respuesta
    reconciled_flows = []
    for flow in request.flows:
        q_idx = var_map[f"q_{flow.id}"]
        cp_idx = var_map[f"cp_{flow.id}"]
        
        q_val = y_hat[q_idx]
        cp_val = y_hat[cp_idx]
        
        # Calcular componentes derivados
        # Sólido = Q * Cp / 100
        # Agua = Q * (1 - Cp/100)
        t_solido = q_val * (cp_val / 100.0)
        t_agua = q_val * (1.0 - (cp_val / 100.0))
        
        reconciled_flows.append({
            "id": flow.id,
            "corrected_data": {
                "tonelaje_pulpa_corregido": q_val,
                "tonelaje_pulpa_error": 0.0, # Pendiente de cálculo riguroso
                "porcentaje_solidos_corregido": cp_val,
                "porcentaje_solidos_error": 0.0,
                "tonelaje_solido_calculado": t_solido,
                "tonelaje_agua_calculado": t_agua
            }
        })

    report = {
        "methodology": [
            "Balance de Agua y Sólidos realizado mediante optimización cuadrática.",
            "Variables ajustadas: Flujo Másico de Pulpa y Porcentaje de Sólidos (%S).",
            "Restricciones: Conservación de masa total y conservación de masa sólida."
        ],
        "equations": format_water_equations(A_final, b_final, var_map)
    }

    return reconciled_flows, report
