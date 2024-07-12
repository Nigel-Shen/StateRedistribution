import numpy as np
import scipy.sparse.linalg as sp
from grid import Grid
from typing import Callable
from discontinuousGalerkin import *
from plot import *


def get_mv(grid: Grid, g: Callable, h:Callable, srd: bool, flux_type: str, **kwargs):
    beta = kwargs.get("beta", 0.5)
    tau = kwargs.get("tau", 1)
    tau_D = kwargs.get("tau_D", tau)
    def mv(u):
        U = np.reshape(u, (grid.Nc, (grid.P + 1) * (grid.P + 2) // 2))
            
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
        for i in range(len(F)):
            F[i, :] = F[i, :] * np.sum(grid.quad_weights[i]) / (grid.h[0] * grid.h[1])
        
        return F.flatten()

    return mv

def solve_Poisson(gridname: str, Nx: int, Ny: int, P: int, f: Callable, g: Callable, h: Callable, **kwargs):
    basis_mode = kwargs.get("basis_mode", "nontensor")
    srd = kwargs.get("srd", False)
    flux_type = kwargs.get("flux_type", "LDG")
    plot_matrix = kwargs.get("plot_matrix", False)
    plot_sol = kwargs.get("plot_sol", False)
    threshold = kwargs.get("threshold", 0.5)
    beta = kwargs.get("beta", 0.5)
    tau = kwargs.get("tau", 1)
    tau_D = kwargs.get("tau_D", tau)

    grid = Grid(gridname, P)
    grid.find_basis(mode=basis_mode)
    grid.get_mask(Nx=Nx, Ny=Ny)
    grid.preprocess(threshold=threshold)
    grid.find_merged_basis(mode=basis_mode)
    if srd:
        grid.get_coarsen_M()
        grid.get_redistribute_M()

    A = sp.LinearOperator(
        (grid.Nc * (P + 1) * (P + 2) // 2, grid.Nc * (P + 1) * (P + 2) // 2),
        matvec=get_mv(
            grid, g, h, srd, flux_type, beta=beta, tau=tau, tau_D=tau_D
        ),
    )
    b_hat = A * np.zeros(grid.Nc * (P + 1) * (P + 2) // 2)
    
    if plot_matrix:
        Adense = A * np.eye(grid.Nc * (P + 1) * (P + 2) // 2) - np.array([b_hat]).T
        print(f"Condition number of A is {np.linalg.cond(Adense)}")
        plot_stiffness_matrix(Adense)

    def mv_b(u):
        return b_hat
    B = sp.LinearOperator((grid.Nc * (P + 1) * (P + 2) // 2, grid.Nc * (P + 1) * (P + 2) // 2),
            matvec=mv_b
        )

    F = initialize(f, grid)
    for i in range(len(F)):
        F[i, :] = F[i, :] * np.sum(grid.quad_weights[i]) / (grid.h[0] * grid.h[1])
    F = F.flatten() - b_hat
    U, _ = sp.cg(A-B, F, tol=1e-11)
    U = np.reshape(U, (grid.Nc, (P + 1) * (P + 2) // 2))

    if plot_sol:
        Uval = evaluate_at_nodes(U, grid)
        plot(grid, Uval)

    return U, grid