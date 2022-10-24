import math
from dataclasses import dataclass

import numpy as np
from numba import njit, prange
from tqdm import tqdm

two_pi = math.pi * 2


@njit
def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


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
        self.x = (self.ix + 1) * self.dx / 2.  # compute center point
        self.y = (self.iy + 1) * self.dy / 2.  # compute center point
        self.z = (self.iz + 1) * self.dz / 2.  # compute center point
        self.area = self.dx * self.dy
        self.dl = self.dz
        return self.dz

    def distance_sqr_from(self, other_block: 'Block'):
        return (self.x - other_block.x)**2 + (self.y - other_block.y)**2 + (
            self.z - other_block.z)**2

    def set_I(self, I):
        self.I = I

    def set_j(self, j):
        self.I = j * self.area

    def biot_savart(self, other_block: 'Block'):
        r = self.distance_sqr_from(other_block)
        if r < 1e-15:
            return 0.
        H = other_block.I * other_block.area * self.dl / r**2
        return H / (4 * np.pi)

    def ampere_law(self, other_block: 'Block'):
        return ampere_law(other_block.I, (self.x, self.y),
                          (other_block.x, other_block.y))

    def __eq__(self, __o: 'Block') -> bool:
        if (self.ix == __o.ix) and (self.iy == __o.iy) and (self.iz == __o.iz):
            return True
        return False


class Structure:

    def __init__(self,
                 maxX,
                 maxY,
                 maxZ,
                 dx,
                 dy,
                 dz,
                 method: Literal['ampere', 'biot-savart'] = 'ampere') -> None:
        self.maxX = maxX
        self.maxY = maxY
        self.maxZ = maxZ
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.Xsize = int(maxX / dx)
        self.Ysize = int(maxY / dy)
        self.Zsize = max(int(maxZ / dz), 1)
        print(f"Creating {self.Xsize}x{self.Ysize}x{self.Zsize} blocks")
        if (method == 'ampere') and (self.Zsize > 1):
            raise ValueError("Wasting compute with z dim non-zero!")

        self.blocks = self.init_blocks()

    def set_region_I(self, min_y, max_y, I):
        if max_y == -1:
            max_y = self.Ysize
        for yindx in prange(min_y, min(max_y + 1, self.Ysize)):
            for x in range(self.Xsize):
                for zindx in range(self.Zsize):
                    self.blocks[x, yindx, zindx].set_I(I)

    def init_blocks(self):
        null_blocks = np.empty((self.Xsize, self.Ysize, self.Zsize),
                               dtype=Block)
        for ix in prange(self.Xsize):
            for iy in range(self.Ysize):
                for iz in range(self.Zsize):
                    null_blocks[ix, iy, iz] = Block(ix, iy, iz, self.dx,
                                                    self.dy, self.dz)
        return null_blocks

    def iterate_blocks(self):
        for block in self.blocks.flatten():
            yield block

    def compute_blocks(self):
        all_blocks = self.blocks.flatten()
        for i in tqdm(range(len(all_blocks)), mininterval=0.5):
            for j in range(i + 1, len(all_blocks)):
                block = all_blocks[i]
                other_block = all_blocks[j]
                if block == other_block:
                    raise ValueError("Trying to compute the same block")
                # use symmetry to reduce complexity
                block.Hlocal += block.ampere_law(other_block)
                other_block.Hlocal += other_block.ampere_law(block)
