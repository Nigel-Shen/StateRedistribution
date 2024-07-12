import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from matplotlib import rcParams
from grid import Grid

sns.set()
sns.set_style("whitegrid")
rcParams["font.family"] = "serif"
blue_cmap = sns.color_palette("Blues", as_cmap=True)
ice_fire = sns.color_palette("icefire", as_cmap=True)
spectral = sns.color_palette("Spectral", as_cmap=True)

def pyplot_flip(field: np.ndarray) -> np.ndarray:
    return np.flip(field.T, axis=0)

def plot(grid: Grid, fval: list, **kwargs) -> None:
    xlim: int = kwargs.get("xlim", [0, grid.h[0]*grid.Nx])
    ylim: int = kwargs.get("ylim", [0, grid.h[1]*grid.Ny])
    res: int  = kwargs.get("res", 24)
    cmap      = kwargs.get("cmap", ice_fire)
    figsize   = kwargs.get("figsize", (5, 4))
    dpi: int  = kwargs.get("dpi", 250)
    plot_grid: bool = kwargs.get("plot_grid", True)
    plot_range = kwargs.get("plot_range", [slice(0, None), slice(1, None)])
    
    field = grid.evaluate_on_mesh(fval, res)
    vmax = np.max(np.abs(field))
    fig, ax = plt.subplots(dpi=dpi, figsize = figsize)
    im = ax.imshow(
                pyplot_flip(field)[plot_range[0], plot_range[1]],
                cmap=cmap, 
                vmin=-vmax, vmax=vmax, 
                extent=[0, grid.h[0]*grid.Nx, 0, grid.h[1]*grid.Ny]
            )
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    plt.colorbar(im)
    ax.axis("equal")
    ax.grid(False)
    if plot_grid:
        for face in grid.faces:
            face = face[1:]
            for i in range(len(face)-1):
                vx = grid.vertices[[face[i], face[i+1]], 0]
                vy = grid.vertices[[face[i], face[i+1]], 1]
                ax.plot(vx, vy, c="w", linewidth=1)
    plt.show()
    
def plot_error(grid: Grid, fval: list, f_true: callable, **kwargs) -> None:
    xlim: int = kwargs.get("xlim", grid.h[0]*grid.Nx)
    ylim: int = kwargs.get("ylim", grid.h[1]*grid.Ny)
    res: int  = kwargs.get("res", 24)
    cmap      = kwargs.get("cmap", spectral)
    figsize   = kwargs.get("figsize", (5, 4))
    dpi: int  = kwargs.get("dpi", 250)
    plot_grid: bool = kwargs.get("plot_grid", False)
    
    field = []
    for i in range(len(fval)):
        field.append(fval[i] - f_true(grid.quad_nodes[0][i], grid.quad_nodes[1][i]))
    field = np.abs(grid.evaluate_on_mesh(field, res=res)) + 1e-16
    fig, ax = plt.subplots(dpi=dpi, figsize = figsize)
    im = ax.imshow(
                pyplot_flip(field),
                cmap=cmap, 
                extent=[0, xlim, 0, ylim],
                norm=LogNorm(),
            )
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    plt.colorbar(im)
    ax.axis("equal")
    ax.grid(False)
    if plot_grid:
        for face in grid.faces:
            face = face[1:]
            for i in range(len(face)-1):
                vx = grid.vertices[[face[i], face[i+1]], 0]
                vy = grid.vertices[[face[i], face[i+1]], 1]
                ax.plot(vx, vy, c="w", linewidth=1)
    plt.show()
    
def plot_stiffness_matrix(A: np.ndarray):
    fig, ax = plt.subplots(1, 2, dpi=250, figsize=(11, 4))
    sns.heatmap(
        np.abs(A),
        norm=LogNorm(vmin=10, vmax=np.max(np.abs(A))),
        ax=ax[0],
        cmap=blue_cmap,
    )
    sns.heatmap(np.abs(A - A.T), norm=LogNorm(vmin=1e-16), ax=ax[1], cmap=blue_cmap)
    ax[0].axis("equal")
    ax[1].axis("equal")
    plt.show()