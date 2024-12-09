import numpy as np
from numba import jit, njit
from numba.typed import List
from typing import Callable
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from copy import copy
from scipy.sparse import csr_matrix
from vandermonde_t import get_lagrange_legendre, get_legendre_vandermonde


def read_file(filename: str, skip_header: bool = True) -> list:
    file1 = open(filename, "r")
    Lines = file1.readlines()

    if skip_header:
        Lines = Lines[1:] # skip the first line if it is a header

    lines = []
    for line in Lines:
        lines.append([eval(i) for i in line.split(" ")])
    return lines

class Grid:
    def __init__(self, filename: str, P: int, **kwargs) -> None:
        self.cells: list = read_file(filename + "_c.txt") # face indices of each cell
        self.Nc = len(self.cells) 
        self.vertices: np.ndarray = np.array(read_file(filename + "_v.txt")) # coordinates of vertices
        self.faces: list = read_file(filename + "_f.txt") # vertex indices of each face
        self.faces_lr: list = read_file(filename + "_flr.txt") # cells indices of each face, -1 -> Dirichlet BC, -2 -> Neumann BC
        self.cut_flag: list = [i[0] for i in self.cells]

        quad = read_file(filename + "_q.txt", skip_header=False)
        quad_nodes_x, quad_nodes_y, quad_weights = quad[0::3], quad[1::3], quad[2::3]

        self.quad_nodes: list = [quad_nodes_x, quad_nodes_y]
        self.quad_weights: list = quad_weights
        self.P: int = P

        h = kwargs.get("hxy", None)
        if h == None:
            sample_cell_idx: int = self.cut_flag.index(0)
            fidx: list = self.cells[sample_cell_idx][1:]
            vidx: list = []
            for idx in fidx:
                vidx += self.faces[idx][1:]
            vidx = list(set(vidx))
            sample_xc = self.vertices[vidx, :]

            x_max, x_min, y_max, y_min = (
                np.max(sample_xc[:, 0]),
                np.min(sample_xc[:, 0]),
                np.max(sample_xc[:, 1]),
                np.min(sample_xc[:, 1]),
            )
            self.h = np.array([x_max - x_min, y_max - y_min])
        else:
            self.h = h

    def get_vertices(self, cell_idx: int) -> tuple:
        face_idx = self.cells[cell_idx][1:]
        vertex_idx = []
        reverse_flag = []
        for idx in face_idx:
            vidx = self.faces[idx][1:]
            if len(vertex_idx) != 0:
                if vertex_idx[-1][0] in vidx:
                    vertex_idx[-1].reverse()
                    reverse_flag.append(True)
                else:
                    reverse_flag.append(False)
            vertex_idx.append(vidx)
        if vertex_idx[-1][0] in vertex_idx[0]:
            vertex_idx[-1].reverse()
            reverse_flag.append(True)
        else:
            reverse_flag.append(False)
        xc = []
        for idx in vertex_idx:
            for i in range(len(idx) - 1):
                xc.append(self.vertices[idx[i], :])
        return np.array(xc), reverse_flag

    def find_basis(self, **kwargs) -> None:
        self.basis_vols: list = []
        self.dxbasis_vols: list = []
        self.dybasis_vols: list = []
        self.reverse_flag: list = []
        self.quad_nodes_surf: list = []
        self.quad_weights_surf: list = []
        self.basis_edges: list = []
        self.normals: list = []
        
        for face in self.faces:
            vertices = self.vertices[face[1:], :]
            Q = len(face) - 2
            N = int(np.ceil(((self.P + 1) * (Q + 1) + 1) / 2))
            Phi, dPhi, ws = get_lagrange_legendre(Q, N)
            x_gl = Phi @ vertices
            normal = dPhi @ np.array([vertices[:, 1], - vertices[:, 0]]).T
            self.quad_nodes_surf.append(x_gl)
            self.quad_weights_surf.append(ws)
            self.normals.append(normal)
            
        for i in range(len(self.cells)):
            xq, yq, wq = (
                self.quad_nodes[0][i],
                self.quad_nodes[1][i],
                self.quad_weights[i],
            )
            xc, reverse_flag = self.get_vertices(i)
            self.reverse_flag.append(reverse_flag)

            x_max, x_min, y_max, y_min = (
                np.max(xc[:, 0]),
                np.min(xc[:, 0]),
                np.max(xc[:, 1]),
                np.min(xc[:, 1]),
            )
            xs = (2 * np.array(xq) - (x_min + x_max)) / (x_max - x_min)
            ys = (2 * np.array(yq) - (y_min + y_max)) / (y_max - y_min)
            
            V, dxV, dyV = get_legendre_vandermonde(xs, ys, self.P, get_gradient=True)
                
            Q, R = np.linalg.qr(np.diag(np.sqrt(wq)) @ V)
            assert np.size(R, 0) == np.size(R, 1), "Not enough quadrature points"
            Rinv = np.linalg.inv(R)
            self.basis_vols.append(np.diag(np.sqrt(np.sum(wq) / wq)) @ Q)
            
            self.dxbasis_vols.append(
                np.sqrt(np.sum(wq)) * 2 * dxV @ Rinv / (x_max - x_min)
            )
            self.dybasis_vols.append(
                np.sqrt(np.sum(wq)) * 2 * dyV @ Rinv / (y_max - y_min)
            )
                
            basis_edge = []
            for j in range(len(self.cells[i][1:])):
                face_idx = self.cells[i][j+1]
                x_edge = copy(self.quad_nodes_surf[face_idx])
                x_edge[:, 0] = (2 * x_edge[:, 0] - (x_min + x_max)) / (x_max - x_min)
                x_edge[:, 1] = (2 * x_edge[:, 1] - (y_min + y_max)) / (y_max - y_min)
                V_edge = get_legendre_vandermonde(x_edge[:, 0], x_edge[:, 1], self.P)
                basis_edge.append(np.sqrt(np.sum(wq)) * V_edge @ Rinv)
            self.basis_edges.append(basis_edge)
        
    def get_mask(self, **kwargs) -> None:
        self.Nx, self.Ny = (np.max(self.vertices, axis=0) // self.h).astype(int) + 1
        self.Nx = kwargs.get("Nx", self.Nx)
        self.Ny = kwargs.get("Ny", self.Ny)
        self.mask = np.ones((self.Nx, self.Ny), dtype=int) * (-1)
        self.inv_mask = np.zeros((self.Nc, 2))
        for i in range(len(self.cells)):
            xc, _ = self.get_vertices(i)
            xymean = np.mean(xc, axis=0)
            nx, ny = (xymean // self.h).astype(int)
            self.mask[nx, ny] = i
            self.inv_mask[i, :] = np.array([nx, ny], dtype=int)
            
    def preprocess(self, **kwargs) -> None:
        threshold = kwargs.get("threshold", 1)
        mode = kwargs.get("mode", "elliptic")
        self.threshold = threshold
    
        merged_ind: list = []
        overlap: list    = [[] for _ in range(self.Nc)]
        for i in range(self.Nc):
            overlap[i].append(i)
            merged_ind.append([i])
        
        self.alphas = []
        for i in range(len(self.cells)):
            if self.cut_flag[i] == 0:
                alpha = 1
            else:
                perimeter = 0
                volume = np.sum(self.quad_weights[i])
                for face_idx in self.cells[i][1:]:
                    perimeter += np.sum(np.sqrt(np.sum(self.normals[face_idx] ** 2, axis=1)) * self.quad_weights_surf[face_idx])
                if mode == "elliptic":
                    alpha = 16 * volume / perimeter ** 2
                elif mode == "parabolic":
                    alpha = volume / (self.h[0] * self.h[1])
                
            self.alphas.append(alpha)
        
        cut_idx = 0    
        for i in range(len(self.cells)):
            if self.alphas[i] < threshold: # needs stabilization
                overlap[i].append(self.Nc + cut_idx)
                merged_ind.append([i])
                
                # First try normal merging
                mean_normal = np.array([0.0, 0])
                for face_idx in self.cells[i][1:]:
                    face = self.faces[face_idx][1:]
                    if -1 in self.faces_lr[face_idx] or -2 in self.faces_lr[face_idx]:
                        edge_len = np.linalg.norm(self.vertices[face[0], :] - self.vertices[face[-1], :])
                        mean_normal -= edge_len * np.sum(self.normals[face_idx] * self.quad_weights_surf[face_idx][:, None], axis=0)
                if np.abs(mean_normal[0]) >= np.abs(mean_normal[1]):
                    mean_normal[1] = 0
                else:
                    mean_normal[0] = 0
                direction = mean_normal / np.linalg.norm(mean_normal)
                neighbor_loc = (self.inv_mask[i, :] + direction).astype(int)
                neighbor_idx = self.mask[neighbor_loc[0], neighbor_loc[1]]
                
                
                if np.sum(self.quad_weights[i]) + np.sum(self.quad_weights[neighbor_idx]) >= self.h[0] * self.h[1] * threshold: # needs modification
                    overlap[neighbor_idx].append(self.Nc + cut_idx)
                    merged_ind[-1].append(neighbor_idx)
                else: # If normal merging does not work, try central merging
                    for m in range(-1, 2):
                        for n in range(-1, 2):
                            if m == 0 and n == 0:
                                continue
                            else:
                                neighbor_loc = (self.inv_mask[i, :] + np.array([m, n])).astype(int)
                                if neighbor_loc[0] < 0 or neighbor_loc[0] >= self.Nx or neighbor_loc[1] < 0 or neighbor_loc[1] >= self.Ny:
                                    continue
                                neighbor_idx = self.mask[neighbor_loc[0], neighbor_loc[1]]
                                if neighbor_idx == -1:
                                    continue
                                else:
                                    overlap[neighbor_idx].append(self.Nc + cut_idx)
                                    merged_ind[-1].append(neighbor_idx)
                                    
                cut_idx += 1
                
        self.merged_ind = merged_ind
        self.overlap = overlap
        self.weight = csr_matrix((self.Nc, len(self.merged_ind)))

        for i in range(self.Nc):
            if self.alphas[i] <= threshold:
                self.weight[i, i] = np.min([1, self.alphas[i]]) # np.sum(self.quad_weights[i]) / (self.h[0] * self.h[1])
            else:
                self.weight[i, i] = 1 / len(self.overlap[i])
            for j in self.overlap[i]:
                if i != j:
                    self.weight[i, j] = (1 - self.weight[i, i]) / (
                        len(self.overlap[i]) - 1
                    )
        
        # for i in range(self.Nc):
        #     for j in self.overlap[i]:
        #         if i != j:
        #             self.weight[i, j] = (1 - self.alphas[merged_ind[j][0]] / (threshold)) / (
        #                 len(self.overlap[i]) - 1
        #             )
        #     self.weight[i, i] = 1 - np.sum(self.weight[i, :])
                    
    def find_merged_basis(self, **kwargs) -> None:
        
        self.mbasis_vols = []
        for i in range(len(self.merged_ind)):
            if self.alphas[self.merged_ind[i][0]] >= self.threshold:
                self.mbasis_vols.append(self.basis_vols[self.merged_ind[i][0]])
            elif len(self.merged_ind[i]) == 1:
                self.mbasis_vols.append(self.basis_vols[self.merged_ind[i][0]])
            else:
                xq = np.concatenate([self.quad_nodes[0][j] for j in self.merged_ind[i]])
                yq = np.concatenate([self.quad_nodes[1][j] for j in self.merged_ind[i]])
                wq = np.concatenate([np.array(self.quad_weights[j]) * self.weight[j, i] for j in self.merged_ind[i]])
                
                xc = np.concatenate([self.get_vertices(j)[0] for j in self.merged_ind[i]])

                x_max, x_min, y_max, y_min = (
                    np.max(xc[:, 0]),
                    np.min(xc[:, 0]),
                    np.max(xc[:, 1]),
                    np.min(xc[:, 1]),
                )
                xs = (2 * np.array(xq) - (x_min + x_max)) / (x_max - x_min)
                ys = (2 * np.array(yq) - (y_min + y_max)) / (y_max - y_min)
                
                V = get_legendre_vandermonde(xs, ys, self.P)
                    
                Q, _ = np.linalg.qr(np.diag(np.sqrt(wq)) @ V)
                self.mbasis_vols.append(np.diag(np.sqrt(np.sum(wq) / wq)) @ Q)
    
    def get_coarsen_M(self) -> None:
        coarsen_M: list = []
        p = (self.P + 1) * (self.P + 2) // 2
        for i in range(len(self.merged_ind)):  # loop over merged neighborhoods
            M = np.zeros((p, p * len(self.merged_ind[i])))
            H = 0
            row_idx = [0, 0]
            for K in range(len(self.merged_ind[i])):  # loop over unmerged elements inside the neighborhood
                k       = self.merged_ind[i][K]
                row_idx = [row_idx[1], row_idx[1] + len(self.quad_weights[k])]
                H       += np.sum(self.quad_weights[k]) * self.weight[k, i]
                M[:, K * p : (K + 1) * p] = (
                    self.weight[k, i] * self.mbasis_vols[i][row_idx[0]:row_idx[1], :].T @ np.diag(self.quad_weights[k]) @ self.basis_vols[k]
                )
            
            if H != 0:
                coarsen_M.append(M / H)
            else:
                coarsen_M.append(np.zeros((p, p * len(self.merged_ind[i]))))
        self.coarsen_M = coarsen_M

    def Coarsen(self, U: np.ndarray) -> np.ndarray:
        p: int          = (self.P + 1) * (self.P + 2) // 2
        Q: np.ndarray   = np.zeros((len(self.merged_ind), p))
        
        for i in range(len(self.merged_ind)):
            c = np.zeros(p * len(self.merged_ind[i]))
            for K in range(len(self.merged_ind[i])):
                k = self.merged_ind[i][K]
                c[K * p : (K + 1) * p] = U[k]
            Q[i, :] = self.coarsen_M[i] @ c
            
        return Q 
        
    def get_redistribute_M(self) -> None:
        redistribute_M: list = []
        p = (self.P + 1) * (self.P + 2) // 2
        for i in range(self.Nc):  # loop over elements
            M = np.zeros((p, p * len(self.overlap[i])))
            H = np.sum(self.quad_weights[i])
            for K in range(len(self.overlap[i])):
                k = self.overlap[i][K]  # the index of current merged neighborhood
                I = self.merged_ind[k].index(i)
                row_idx = [0, 0]
                for n in range(I + 1):
                    row_idx = [row_idx[1], row_idx[1] + len(self.quad_weights[self.merged_ind[k][n]])]
                M[:, K * p : (K + 1) * p] = (
                    self.weight[i, k] * self.basis_vols[i].T @ np.diag(self.quad_weights[i]) @ self.mbasis_vols[k][row_idx[0]:row_idx[1], :]
                )
            redistribute_M.append(M / H)
        self.redistribute_M = redistribute_M


    def Redistribute(self, Q: np.ndarray) -> np.ndarray:
        p: int            = (self.P + 1) * (self.P + 2) // 2
        U_SRD: np.ndarray = np.zeros((self.Nc, p))
        for i in range(self.Nc):
            c = np.zeros(p * len(self.overlap[i]))
            for k in range(len(self.overlap[i])):
                c[k * p: (k + 1) * p] = Q[self.overlap[i][k]]
            U_SRD[i, :] = self.redistribute_M[i] @ c
        return U_SRD
    
    def evaluate_on_mesh(self, f: list|Callable, res: float) -> np.ndarray:
        fvals = np.zeros((self.Nx * res, self.Ny * res))
        xs = np.linspace(-1, 1, res, endpoint=False) + 1 / res
        xs, ys = np.meshgrid(xs, xs)
        xs = xs.flatten()
        ys = ys.flatten()

        for i in range(self.Nx):
            for j in range(self.Ny):
                if self.mask[i, j] == -1:
                    continue
                else:
                    idx = self.mask[i, j]
                    if isinstance(f, list):
                        f_local = np.array(f[idx])
                        xnodes, ynodes = self.quad_nodes[0][idx], self.quad_nodes[1][idx]

                        x_max, x_min, y_max, y_min = (
                            (i + 1) * self.h[0],
                            i * self.h[0],
                            (j + 1) * self.h[1],
                            j * self.h[1],
                        )
                        xnodes = (2 * np.array(xnodes) - (x_min + x_max)) / (x_max - x_min)
                        ynodes = (2 * np.array(ynodes) - (y_min + y_max)) / (y_max - y_min)
                        Vnodes = get_legendre_vandermonde(xnodes, ynodes, self.P)
                        Vvals = get_legendre_vandermonde(xs, ys, self.P)
                        
                        fval = Vvals @ (np.linalg.pinv(Vnodes) @ f_local)
                        if self.cut_flag[idx] == 1:
                            vc, _ = self.get_vertices(idx)
                            vc[:, 0] = (2 * vc[:, 0] - (x_min + x_max)) / (x_max - x_min)
                            vc[:, 1] = (2 * vc[:, 1] - (y_min + y_max)) / (y_max - y_min)
                            polygon = Polygon(vc)
                            for k in range(len(fval)):
                                point = Point(xs[k], ys[k])
                                if not polygon.contains(point):
                                    fval[k] = 0
                        fval = fval.reshape(res, res).T
                        fvals[i * res : (i + 1) * res, j * res : (j + 1) * res] = fval
                    else:
                        x_max, x_min, y_max, y_min = (
                            (i + 1) * self.h[0],
                            i * self.h[0],
                            (j + 1) * self.h[1],
                            j * self.h[1],
                        )
                        xnodes = (xs * (x_max - x_min) + (x_min + x_max)) / 2
                        ynodes = (ys * (y_max - y_min) + (y_min + y_max)) / 2
                        fval = f(xnodes, ynodes)
                        if self.cut_flag[idx] == 1:
                            vc, _ = self.get_vertices(idx)
                            vc[:, 0] = (2 * vc[:, 0] - (x_min + x_max)) / (x_max - x_min)
                            vc[:, 1] = (2 * vc[:, 1] - (y_min + y_max)) / (y_max - y_min)
                            polygon = Polygon(vc)
                            for k in range(len(fval)):
                                point = Point(xs[k], ys[k])
                                if not polygon.contains(point):
                                    fval[k] = 0
                        fval = fval.reshape(res, res).T
                        fvals[i * res : (i + 1) * res, j * res : (j + 1) * res] = fval
        return fvals