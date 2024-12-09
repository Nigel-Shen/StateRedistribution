import numpy as np
import scipy.sparse.linalg as sp
from rich import *
from grid import Grid
from typing import Callable
from discontinuousGalerkin import *
from plot import *


def get_mv(grid: Grid, g: Callable, h:Callable, srd: bool, flux_type: str, **kwargs):
    beta = kwargs.get("beta", 0.5)
    tau = kwargs.get("tau", 1)
    tau_D = kwargs.get("tau_D", tau)
    mode = kwargs.get("mode", "elliptic")
    def mv(u):
        U = np.reshape(u, (grid.Nc, (grid.P + 1) * (grid.P + 2) // 2))
        
        if mode == "parabolic" and srd:
            U = grid.Redistribute(grid.Coarsen(U))
            
        sigma_x, sigma_y = compute_sigma(U, grid, g, h, flux_type=flux_type, beta=beta)
        if flux_type == "IP":
            Ugrad_x, Ugrad_y = compute_Ugrad(U, grid)
            
        if srd:
            sigma_x = grid.Redistribute(grid.Coarsen(sigma_x))
            sigma_y = grid.Redistribute(grid.Coarsen(sigma_y))
            if flux_type == "IP":
                Ugrad_x = grid.Redistribute(grid.Coarsen(Ugrad_x))
                Ugrad_y = grid.Redistribute(grid.Coarsen(Ugrad_y))
                
        if flux_type == "IP":
            F = compute_F(U, sigma_x, sigma_y, grid, g, h, flux_type=flux_type, beta=beta, tau=tau, tau_D=tau_D, Ugrad_x=Ugrad_x, Ugrad_y=Ugrad_y)
        else:
            F = compute_F(U, sigma_x, sigma_y, grid, g, h, flux_type=flux_type, beta=beta, tau=tau, tau_D=tau_D)
        
        if mode == "parabolic" and srd:
            F = grid.Redistribute(grid.Coarsen(F))
        
        for i in range(len(F)):
            F[i, :] = F[i, :] * np.sum(grid.quad_weights[i]) / (grid.h[0] * grid.h[1])
        
        return F.flatten()

    return mv

def solve_Poisson(gridname: str, Nx: int, Ny: int, P: int, f: Callable, g: Callable, h: Callable, **kwargs):
    srd = kwargs.get("srd", False)
    flux_type = kwargs.get("flux_type", "LDG")
    plot_matrix = kwargs.get("plot_matrix", False)
    plot_sol = kwargs.get("plot_sol", False)
    threshold = kwargs.get("threshold", 0.5)
    beta = kwargs.get("beta", 0.5)
    tau = kwargs.get("tau", 1)
    tau_D = kwargs.get("tau_D", tau)
    hxy = kwargs.get("hxy", None)
    mode = kwargs.get("mode", "elliptic")

    grid, A, b_hat = form_linear_system(gridname, Nx, Ny, P, g, h, 
                                        srd=srd, flux_type=flux_type, 
                                        plot_matrix=plot_matrix, 
                                        threshold=threshold, 
                                        beta=beta, tau=tau, tau_D=tau_D, hxy=hxy, mode=mode)

    F = initialize(f, grid)
    for i in range(len(F)):
        F[i, :] = F[i, :] * np.sum(grid.quad_weights[i]) / (grid.h[0] * grid.h[1])
    F = F.flatten() - b_hat
    U, _ = sp.cg(A, F, tol=1e-11)
    U = np.reshape(U, (grid.Nc, (P + 1) * (P + 2) // 2))

    if plot_sol:
        Uval = evaluate_at_nodes(U, grid)
        plot(grid, Uval)

    return U, grid

def form_linear_system(gridname: str, Nx: int, Ny: int, P: int, g: Callable, h: Callable, **kwargs) -> tuple:
    srd = kwargs.get("srd", False)
    flux_type = kwargs.get("flux_type", "LDG")
    plot_matrix = kwargs.get("plot_matrix", False)
    threshold = kwargs.get("threshold", 0.5)
    beta = kwargs.get("beta", 0.5)
    tau = kwargs.get("tau", 1)
    tau_D = kwargs.get("tau_D", tau)
    hxy = kwargs.get("hxy", None)
    mode = kwargs.get("mode", "elliptic")
    
    grid = Grid(gridname, P, hxy=hxy)
    grid.find_basis()
    grid.get_mask(Nx=Nx, Ny=Ny)
    grid.preprocess(threshold=threshold, mode=mode)
    grid.find_merged_basis()
    if srd:
        grid.get_coarsen_M()
        grid.get_redistribute_M()
    A = sp.LinearOperator(
        (grid.Nc * (P + 1) * (P + 2) // 2, grid.Nc * (P + 1) * (P + 2) // 2),
        matvec=get_mv(
            grid, g, h, srd, flux_type, beta=beta, tau=tau, tau_D=tau_D, mode=mode
        ),
    )
    b_hat = A * np.zeros(grid.Nc * (P + 1) * (P + 2) // 2)

    def mv_b(u):
        return b_hat
    
    B = sp.LinearOperator((grid.Nc * (P + 1) * (P + 2) // 2, grid.Nc * (P + 1) * (P + 2) // 2),
            matvec=mv_b
        )
    
    A = A - B
    if plot_matrix:
        Adense = A * np.eye(grid.Nc * (P + 1) * (P + 2) // 2)
        print(f"Condition number of A is {np.linalg.cond(Adense)}")
        plot_stiffness_matrix(Adense)
    
    return grid, A, b_hat


def Euler_forward(dt: float, U: np.ndarray, eps: float, A: sp.LinearOperator, b_hat: np.ndarray, grid: Grid) -> np.ndarray:
    F = A * U.flatten() + b_hat
    F = np.reshape(F, (grid.Nc, (grid.P + 1) * (grid.P + 2) // 2))
    for i in range(len(F)):
        F[i, :] = F[i, :] * (grid.h[0] * grid.h[1]) / np.sum(grid.quad_weights[i])
    new_U = U - F * eps * dt
    return new_U

def RK2(dt: float, U: np.ndarray, eps: float, A: sp.LinearOperator, b_hat: np.ndarray, grid: Grid) -> np.ndarray:
    U_1 = Euler_forward(dt, U, eps, A, b_hat, grid)
    new_U = U / 2 + Euler_forward(dt, U_1, eps, A, b_hat, grid) / 2
    return new_U

def RK3(dt: float, U: np.ndarray, eps: float, A: sp.LinearOperator, b_hat: np.ndarray, grid: Grid) -> np.ndarray:
    U_1 = Euler_forward(dt, U, eps, A, b_hat, grid)
    U_2 = 3 * U / 4 + Euler_forward(dt, U_1, eps, A, b_hat, grid) / 4
    new_U = U / 3 + 2 * Euler_forward(dt, U_2, eps, A, b_hat, grid) / 3
    return new_U