from typing import final
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.fft import fft
import pandas as pd
import time as tm

final_PIMM = []
start = tm.time()
with open('PIMM_res.csv', "r") as f:
    lines = f.read().split("\n")
    for line in lines:
        if not line:
            continue
        els = np.asarray(line.split(";")[1:], dtype=np.float32)
        yf = np.abs(fft(els))
        final_PIMM.append(yf[1:len(yf) // 2])
final_PIMM = np.asarray(final_PIMM)
stop = tm.time()
print(f"Time: {stop-start:.2f}")
fig, ax = plt.subplots(figsize=(10, 10))
im = ax.imshow(final_PIMM.T,
               origin='lower',
               vmin=np.min(final_PIMM),
               vmax=np.mean(final_PIMM) + 3*np.std(final_PIMM))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
print(np.asarray(final_PIMM).shape)
print(final_PIMM.dtype)
fig.savefig("PIMM.png")
