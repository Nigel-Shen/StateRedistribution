import numpy as np
from grid import Grid
from typing import Callable

def initialize(f: Callable, grid: Grid) -> np.ndarray:
    F = np.zeros((len(grid.cells), (grid.P + 1) * (grid.P + 2) // 2))
    for i in range(len(grid.cells)):
        fval = f(grid.quad_nodes[0][i], grid.quad_nodes[1][i])
        F[i, :] = (
            (fval * np.array(grid.quad_weights[i]))
            @ grid.basis_vols[i]
            / np.sum(grid.quad_weights[i])
        )
    return F

def evaluate_at_nodes(F: list, grid: Grid) -> list:
    fval = []
    for i in range(len(grid.cells)):
        fval.append(F[i] @ grid.basis_vols[i].T)
    return fval

def compute_error(grid: Grid, F: np.ndarray, f: Callable) -> float:
    fval = evaluate_at_nodes(F, grid)
    L2error = 0
    Linferror = 0
    for i in range(len(grid.cells)):
        fval[i] = fval[i] - f(grid.quad_nodes[0][i], grid.quad_nodes[1][i])
        L2error += np.sum(fval[i] ** 2 * grid.quad_weights[i])
        linferror = np.max(np.abs(fval[i]))
        if linferror > Linferror:
            Linferror = linferror
    L2error = np.sqrt(L2error)
    return L2error, Linferror

def volume_integral(U: np.ndarray, grid: Grid, cell_idx: int, mode: str) -> np.ndarray:
    if mode == "scalar":
        basis_vol   = grid.basis_vols[cell_idx]
        dbasis_vol  = np.concatenate((grid.dxbasis_vols[cell_idx].T, grid.dybasis_vols[cell_idx].T))
        ws          = grid.quad_weights[cell_idx]
        K           = dbasis_vol @ np.diag(ws) @ basis_vol / np.sum(ws)
        return K @ U[cell_idx, :]
    elif mode == "vector":
        basis_vol       = grid.basis_vols[cell_idx]
        dims            = np.shape(basis_vol)
        basis_vol_db    = np.zeros((dims[0]*2, dims[1]*2))
        basis_vol_db[0:dims[0], 0:dims[1]] = basis_vol
        basis_vol_db[dims[0]:, dims[1]:]   = basis_vol
        dbasis_vol      = np.concatenate((grid.dxbasis_vols[cell_idx].T, grid.dybasis_vols[cell_idx].T), axis=1)
        ws              = grid.quad_weights[cell_idx] * 2 # Replicate weights, NOT multiplying by 2!
        K               = dbasis_vol @ np.diag(ws) @ basis_vol_db / np.sum(grid.quad_weights[cell_idx])
        return K @ np.concatenate((U[0][cell_idx, :], U[1][cell_idx, :]))

def compute_flux_U(U: np.ndarray, grid: Grid, face_idx: int, flux_type: str, g: Callable, h: Callable, beta: np.float64 = 0.5) -> np.ndarray:
    cells_idx = grid.faces_lr[face_idx]
    U_pm = []
    if cells_idx[0] >= 0:
        idx:int          = grid.cells[cells_idx[0]][1:].index(face_idx)
        basis_edge_plus  = grid.basis_edges[cells_idx[0]][idx]
        U_pm.append(basis_edge_plus @ U[cells_idx[0], :])
    if cells_idx[1] >= 0:
        idx:int          = grid.cells[cells_idx[1]][1:].index(face_idx)
        basis_edge_minus = grid.basis_edges[cells_idx[1]][idx]
        U_pm.append(basis_edge_minus @ U[cells_idx[1], :])
    U_pm = np.array(U_pm)
    
    match flux_type:
        case "LDG":
            if ((- 2) in cells_idx) or ((-1) in cells_idx): # Dirichlet
                flux_U = g(grid.quad_nodes_surf[face_idx])
            # elif (- 1) in cells_idx: # Neumann
                # flux_U = U_pm[0]
            else:
                flux_U = np.sum(U_pm, axis=0) / 2 - beta * np.diff(U_pm, axis=0)[0]
                
        case "central":
            if ((- 2) in cells_idx) or ((-1) in cells_idx): # Dirichlet
                flux_U = g(grid.quad_nodes_surf[face_idx])
            # elif (- 1) in cells_idx: # Neumann
                # flux_U = U_pm[0]
            else:
                flux_U = np.sum(U_pm, axis=0) / 2
                
        case "IP":
            if ((- 2) in cells_idx) or ((-1) in cells_idx): # Dirichlet
                flux_U = g(grid.quad_nodes_surf[face_idx])
            # elif (- 1) in cells_idx: # Neumann
                # flux_U = U_pm[0]
            else:
                flux_U = np.sum(U_pm, axis=0) / 2
    
    return flux_U


def compute_flux_sigma(sigma_x: np.ndarray, sigma_y: np.ndarray, U: np.ndarray, 
                       grid: Grid, face_idx: int, flux_type: str, g: Callable, h: Callable, 
                       beta: float = 0.5, tau: float = 1.0, tau_D: float = 1.0) -> np.ndarray:
    
    cells_idx = grid.faces_lr[face_idx]
    v_idx = grid.faces[face_idx]
    h_e = np.linalg.norm(grid.vertices[v_idx[1]] - grid.vertices[v_idx[-1]])
    tau = tau / h_e
    tau_D = tau_D / h_e
    sigma_x_pm = []
    sigma_y_pm = []
    U_pm = []
    
    if cells_idx[0] >= 0:
        idx:int          = grid.cells[cells_idx[0]][1:].index(face_idx)
        basis_edge_plus  = grid.basis_edges[cells_idx[0]][idx]
        sigma_x_pm.append(basis_edge_plus @ sigma_x[cells_idx[0], :])
        sigma_y_pm.append(basis_edge_plus @ sigma_y[cells_idx[0], :])
        U_pm.append(basis_edge_plus @ U[cells_idx[0], :])
    if cells_idx[1] >= 0:
        idx:int          = grid.cells[cells_idx[1]][1:].index(face_idx)
        basis_edge_minus = grid.basis_edges[cells_idx[1]][idx]
        sigma_x_pm.append(basis_edge_minus @ sigma_x[cells_idx[1], :])
        sigma_y_pm.append(basis_edge_minus @ sigma_y[cells_idx[1], :])
        U_pm.append(basis_edge_minus @ U[cells_idx[1], :])
        
    sigma_x_pm = np.array(sigma_x_pm)
    sigma_y_pm = np.array(sigma_y_pm)
    U_pm = np.array(U_pm)
    normal = grid.normals[face_idx] / np.linalg.norm(grid.normals[face_idx], axis=1)[:, None]
    match flux_type:
        case "LDG":
            if ((-2) in cells_idx) or ((-1) in cells_idx): # Dirichlet
                flux_sigma_x = sigma_x_pm[0, :] - tau_D * (U_pm[0, :] - g(grid.quad_nodes_surf[face_idx])) * normal[:, 0]
                flux_sigma_y = sigma_y_pm[0, :] - tau_D * (U_pm[0, :] - g(grid.quad_nodes_surf[face_idx])) * normal[:, 1]
            else:
                jump_sigma = np.diff(sigma_x_pm, axis=0)[0] * normal[:, 0] + np.diff(sigma_y_pm, axis=0)[0] * normal[:, 1]
                jump_U = np.diff(U_pm, axis=0)[0, :, None] * normal
                flux_sigma_x = np.sum(sigma_x_pm, axis=0) / 2 + beta * jump_sigma * normal[:, 0] - tau * jump_U[:, 0]
                flux_sigma_y = np.sum(sigma_y_pm, axis=0) / 2 + beta * jump_sigma * normal[:, 1] - tau * jump_U[:, 1]
                
        case "central":
            if ((-2) in cells_idx) or ((-1) in cells_idx): # Dirichlet
                flux_sigma_x = sigma_x_pm[0, :] - tau_D * (U_pm[0, :] - g(grid.quad_nodes_surf[face_idx])) * normal[:, 0]
                flux_sigma_y = sigma_y_pm[0, :] - tau_D * (U_pm[0, :] - g(grid.quad_nodes_surf[face_idx])) * normal[:, 1]
            else:
                jump_sigma = np.diff(sigma_x_pm, axis=0)[0] * normal[:, 0] + np.diff(sigma_y_pm, axis=0)[0] * normal[:, 1]
                jump_U = np.diff(U_pm, axis=0)[0, :, None] * normal
                flux_sigma_x = np.sum(sigma_x_pm, axis=0) / 2 - tau * jump_U[:, 0]
                flux_sigma_y = np.sum(sigma_y_pm, axis=0) / 2 - tau * jump_U[:, 1]
        
        case "IP":
            if ((-2) in cells_idx) or ((-1) in cells_idx): # Dirichlet
                flux_sigma_x = sigma_x_pm[0, :] - tau_D * (U_pm[0, :] - g(grid.quad_nodes_surf[face_idx])) * normal[:, 0]
                flux_sigma_y = sigma_y_pm[0, :] - tau_D * (U_pm[0, :] - g(grid.quad_nodes_surf[face_idx])) * normal[:, 1]
            else:
                jump_sigma = np.diff(sigma_x_pm, axis=0)[0] * normal[:, 0] + np.diff(sigma_y_pm, axis=0)[0] * normal[:, 1]
                jump_U = np.diff(U_pm, axis=0)[0, :, None] * normal
                flux_sigma_x = np.sum(sigma_x_pm, axis=0) / 2 - tau * jump_U[:, 0]
                flux_sigma_y = np.sum(sigma_y_pm, axis=0) / 2 - tau * jump_U[:, 1]
                
                
    return flux_sigma_x, flux_sigma_y

def surface_integral(flux_U: list, grid: Grid, cell_idx: int, mode: str) -> np.ndarray:
    if mode == "scalar":
        integral    = 0
        basis_edge  = grid.basis_edges[cell_idx]
        face_idx    = grid.cells[cell_idx][1:]
        for i in range(len(face_idx)):
            M1  = basis_edge[i] * grid.normals[face_idx[i]][:, 0, None]
            M2  = basis_edge[i] * grid.normals[face_idx[i]][:, 1, None]
            M   = np.concatenate((M1, M2), axis=1)
            if grid.reverse_flag[cell_idx][i]:
                M = - M
            ws = grid.quad_weights_surf[face_idx[i]]
            area = np.sum(grid.quad_weights[cell_idx])
            integral += M.T @ np.diag(ws) @ flux_U[face_idx[i]] / area
        return integral
    if mode == "vector":
        integral    = 0
        basis_edge  = grid.basis_edges[cell_idx]
        face_idx    = grid.cells[cell_idx][1:]
        for i in range(len(face_idx)):
            normal = grid.normals[face_idx[i]]
            M = basis_edge[i]
            if grid.reverse_flag[cell_idx][i]:
                M = - M
            ws = grid.quad_weights_surf[face_idx[i]]
            area = np.sum(grid.quad_weights[cell_idx])
            integral += M.T @ np.diag(ws) @ (flux_U[0][face_idx[i]] * normal[:, 0] + flux_U[1][face_idx[i]] * normal[:, 1]) / area
        return integral
 
def compute_sigma(U: list, grid: Grid, g: Callable, h: Callable, flux_type: str = "LDG", beta: float = 0.5) -> tuple:
    sigma_x = np.zeros((len(grid.cells), (grid.P + 1) * (grid.P + 2) // 2))
    sigma_y = np.zeros((len(grid.cells), (grid.P + 1) * (grid.P + 2) // 2))
    
    flux_U = []
    for face_idx in range(len(grid.faces)):
        flux_U.append(compute_flux_U(U, grid, face_idx, flux_type, g, h, beta=beta))
    
    for cell_idx in range(len(grid.cells)):
        v_int = volume_integral(U, grid, cell_idx, "scalar")
        s_int = surface_integral(flux_U, grid, cell_idx, "scalar")
        sigma_local = s_int - v_int 
        L = len(sigma_local)
        sigma_x[cell_idx, :] = sigma_local[0:L//2]
        sigma_y[cell_idx, :] = sigma_local[L//2:]
    
    return sigma_x, sigma_y

def compute_F(U: np.ndarray, sigma_x: np.ndarray, sigma_y: np.ndarray, grid: Grid, g: Callable, h: Callable, 
              flux_type: str = "LDG", beta: float = 0.5, tau: float = 1.0, tau_D: float = 1.0, Ugrad_x=None, Ugrad_y=None) -> np.ndarray:
    
    F = np.zeros((len(grid.cells), (grid.P + 1) * (grid.P + 2) // 2))
    flux_sigma = [[], []]
    for face_idx in range(len(grid.faces)):
        if flux_type == "IP":
            fsigma_x, fsigma_y = compute_flux_sigma(Ugrad_x, Ugrad_y, U, grid, face_idx, flux_type, g, h, beta=beta, tau=tau, tau_D=tau_D)
        else:
            fsigma_x, fsigma_y = compute_flux_sigma(sigma_x, sigma_y, U, grid, face_idx, flux_type, g, h, beta=beta, tau=tau, tau_D=tau_D)
        flux_sigma[0].append(fsigma_x)
        flux_sigma[1].append(fsigma_y)
        
    for cell_idx in range(len(grid.cells)):
        v_int = volume_integral([sigma_x, sigma_y], grid, cell_idx, "vector")
        s_int = surface_integral(flux_sigma, grid, cell_idx, "vector")
        F_local = - s_int + v_int
        F[cell_idx, :] = F_local
    return F

def compute_Ugrad(U: np.ndarray, grid: Grid) -> tuple:
    Ugrad_x = np.zeros((len(grid.cells), (grid.P + 1) * (grid.P + 2) // 2))
    Ugrad_y = np.zeros((len(grid.cells), (grid.P + 1) * (grid.P + 2) // 2))
    for cell_idx in range(len(grid.cells)):
        basis_vol   = grid.basis_vols[cell_idx]
        dxbasis_vol = grid.dxbasis_vols[cell_idx]
        dybasis_vol = grid.dybasis_vols[cell_idx]
        ws          = grid.quad_weights[cell_idx]
        Ugrad_x[cell_idx, :] = basis_vol.T @ np.diag(ws) @ dxbasis_vol @ U[cell_idx] / np.sum(ws)
        Ugrad_y[cell_idx, :] = basis_vol.T @ np.diag(ws) @ dybasis_vol @ U[cell_idx] / np.sum(ws)
    return Ugrad_x, Ugrad_y