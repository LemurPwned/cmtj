import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

df = pd.read_csv('VSD_res.csv', sep=';')

fig, ax = plt.subplots(figsize=(12, 12))

res_img = []
for f in sorted(df['f'].unique(), reverse=True):
    tmp = df.loc[df['f'] == f].sort_values(by='H')
    # print(tmp['H'])
    res_img.append(tmp['Vmix'].tolist())

ax.imshow(res_img)
ax.set_xlabel("H [A/m]")
ax.set_ylabel("Vmix [V]")
ax.set_title("VSD map frequency dispersion")
ax.legend()
fig.savefig("VSD.png")