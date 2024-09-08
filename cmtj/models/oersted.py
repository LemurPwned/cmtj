import math
from dataclasses import dataclass
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from numba import njit, prange
from tqdm import tqdm

two_pi = math.pi * 2


@njit
def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


@njit
def ampere_law(I, p1, p2):
    r = distance(p1, p2)
    return I / (two_pi * r)


@dataclass
class Block:
    ix: int
    iy: int
    iz: int
    dx: float
    dy: float
    dz: float
    Hlocal: float = 0
    I: float = 0

    def __post_init__(self):
        self.x = (self.ix + 1) * self.dx / 2.0  # compute center point
        self.y = (self.iy + 1) * self.dy / 2.0  # compute center point
        self.z = (self.iz + 1) * self.dz / 2.0  # compute center point
        self.area = self.dx * self.dy
        self.dl = self.dz
        return self.dz

    def distance_sqr_from(self, other_block: "Block"):
        return (self.x - other_block.x) ** 2 + (self.y - other_block.y) ** 2 + (self.z - other_block.z) ** 2

    def set_I(self, I):
        self.I = I

    def set_j(self, j):
        self.I = j * self.area

    def add_H(self, H):
        self.Hlocal += H

    def biot_savart(self, other_block: "Block"):
        r = distance((self.x, self.y), (other_block.x, other_block.y))
        if r < 1e-15:
            return 0.0
        H = other_block.I * other_block.area * self.dl / r**2
        return H / (4 * math.pi)

    def ampere_law(self, other_block: "Block"):
        return ampere_law(other_block.I, (self.x, self.y), (other_block.x, other_block.y))

    def __eq__(self, __o: "Block") -> bool:
        return self.ix == __o.ix and self.iy == __o.iy and self.iz == __o.iz


class Structure:
    def __init__(self, maxX, maxY, maxZ, dx, dy, dz, method: Literal["ampere", "biot-savart"] = "ampere") -> None:
        self.maxX = maxX
        self.maxY = maxY
        self.maxZ = maxZ
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.Xsize = math.ceil(maxX / dx)
        self.Ysize = math.ceil(maxY / dy)
        self.Zsize = max(math.ceil(maxZ / dz), 1)
        print(f"Creating {self.Xsize}x{self.Ysize}x{self.Zsize} blocks")
        if (method == "ampere") and (self.Zsize > 1):
            raise ValueError("Wasting compute with z dim non-zero!")

        self.blocks = self.init_blocks()
        self.borders = []
        self.labels = []

    def set_region_I_idx(self, I, min_y_indx, max_y_indx=-1):
        if max_y_indx == -1:
            max_y_indx = self.Ysize
        I_mag = I / (self.Xsize * (min(max_y_indx + 1, self.Ysize) - min_y_indx))
        # print(f"Setting I={I*1e3:.2f}mA in region {min_y_indx}:{max_y_indx}")
        # print(f"Unit I={I_mag*1e6:.2f}uA")
        for yindx in prange(min_y_indx, min(max_y_indx, self.Ysize)):
            for x in range(self.Xsize):
                for zindx in range(self.Zsize):
                    self.blocks[x, yindx, zindx].set_I(I_mag)

    def set_region_I(self, I, min_y, max_y=-1, label=None):
        # convert to index
        self.borders.append(min_y)
        self.labels.append(label)
        min_y_indx = math.ceil(min_y / self.dy)
        min_y_indx, max_y_indx = self.__ycoords2indx(min_y, max_y)
        self.set_region_I_idx(I, min_y_indx, max_y_indx)

    def init_blocks(self):
        null_blocks = np.empty((self.Xsize, self.Ysize, self.Zsize), dtype=Block)
        for ix in prange(self.Xsize):
            for iy in range(self.Ysize):
                for iz in range(self.Zsize):
                    null_blocks[ix, iy, iz] = Block(ix, iy, iz, self.dx, self.dy, self.dz)
        return null_blocks

    def reset(self):
        self.blocks = self.init_blocks()

    def compute_blocks(self):
        all_blocks = self.blocks.flatten()
        for i in tqdm(prange(len(all_blocks)), mininterval=0.5):
            for j in range(i + 1, len(all_blocks)):
                block = all_blocks[i]
                other_block = all_blocks[j]
                # use symmetry to reduce complexity
                Hval1 = block.ampere_law(other_block)
                Hval2 = other_block.ampere_law(block)
                block.Hlocal += Hval1
                other_block.Hlocal += Hval2

    def __ycoords2indx(self, min_coord_y, max_coord_y):
        min_y_indx = math.ceil(min_coord_y / self.dy)
        max_y_indx = self.maxY if max_coord_y == -1 else math.ceil(max_coord_y / self.dy)
        return min_y_indx, max_y_indx

    def get_region_contributions_idx(self, min_y_indx, max_y_indx=-1):
        if max_y_indx == -1:
            max_y_indx = self.Ysize
        H = 0
        for yindx in prange(min_y_indx, min(max_y_indx, self.Ysize)):
            for x in range(self.Xsize):
                for zindx in range(self.Zsize):
                    H += self.blocks[x, yindx, zindx].Hlocal
        return H

    def compute_region_contribution(self, source_min_y, source_max_y, target_min_y, target_max_y):
        source_min_y_indx, source_max_y_indx = self.__ycoords2indx(source_min_y, source_max_y)
        target_min_y_indx, target_max_y_indx = self.__ycoords2indx(target_min_y, target_max_y)

        return self.compute_region_contribution_idx(
            source_min_y_indx, source_max_y_indx, target_min_y_indx, target_max_y_indx
        )

    def compute_region_contribution_idx(
        self, source_min_y_indx, source_max_y_indx, target_min_y_indx, target_max_y_indx
    ):
        print(f"Computing H from {source_min_y_indx}:{source_max_y_indx} to {target_min_y_indx}:{target_max_y_indx}")
        total_H = 0
        block_src = self.blocks[:, source_min_y_indx:source_max_y_indx, :].flatten()
        block_targ = self.blocks[:, target_min_y_indx:target_max_y_indx, :].flatten()
        print(f"Source blocks: {block_src.shape}")

        print(f"Target blocks: {block_targ.shape}")
        for source_block in block_src:
            for target_block in block_targ:
                if source_block == target_block:
                    continue
                H = target_block.ampere_law(source_block)
                total_H += H
        return total_H

    def show_field(self, log=False):
        field = np.zeros((self.Xsize, self.Ysize))
        for block in self.blocks.flatten():
            field[block.ix, block.iy] = block.Hlocal
        with plt.style.context(["nature", "science"]):
            fig, ax = plt.subplots(dpi=300)

            img = ax.pcolormesh(np.log(field).T, cmap="viridis") if log else ax.pcolormesh(field.T, cmap="viridis")
            # add colorbar
            ax.set_xticklabels([f"{x*1e9:.2f}" for x in np.linspace(0, self.Xsize * self.dx, self.Xsize)])
            ax.set_yticklabels([f"{y*1e9:.2f}" for y in np.linspace(0, self.Ysize * self.dy, self.Ysize)])
            ax.set_xlabel("x (nm)")
            ax.set_ylabel("y (nm)")
            # add colorbar
            for unq_border, _label in zip(self.borders, self.labels):
                ax.axhline(unq_border / self.dy, color="crimson")
            fig.colorbar(img, ax=ax)
        return field
