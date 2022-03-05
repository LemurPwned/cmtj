import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.signal.windows import hann


def find_max_f_frequency(freqs: np.ndarray, values: np.ndarray,
                         frequency: float):
    # take between 0 and max
    freqs_pos_freq = np.argwhere((freqs <= 1.001 * frequency)
                                 & (freqs >= 0.999 * frequency)
                                 & (freqs > 0)).ravel()
    max_freq_indx = np.argmax(values[freqs_pos_freq])
    max_value = np.max(values[freqs_pos_freq])
    max_freq = freqs[freqs_pos_freq][max_freq_indx]
    return max_value, max_freq


df = pd.read_csv("Torque_res.csv", sep=";")

if "phase" not in df.columns:
    f = 0.8e9
    tstep = 1e-11
    Hscan = []
    amp_diagram = []
    phase_diagram = []
    for H, Hgroup in df.groupby("H"):
        Hscan.append(H / 1e3)
        t = Hgroup.sort_values(by="indx")["Vmix"].values
        # t *= hann(len(t))
        # print(len(t), H)
        y = fft(t) * (2 / len(t))
        amp = np.abs(y)
        phase = np.angle(y)
        freqs = fftfreq(len(y), tstep)
        y = y[:len(y) // 2]
        freqs = freqs[:len(freqs) // 2]

        max_phase2f, max_freq = find_max_f_frequency(freqs, phase, 2 * f)
        print("Diff", abs(max_freq - 2 * f))
        max_amp1f, _ = find_max_f_frequency(freqs, amp, f)
        max_amp2f, _ = find_max_f_frequency(freqs, amp, 2 * f)

        amp_diagram.append(max_amp1f)
        phase_diagram.append(max_amp2f * np.cos(max_phase2f))
    amp_diagram = np.asarray(amp_diagram)
    phase_diagram = np.asarray(phase_diagram)
elif "frequency" in df.columns:

    cf = 0.01e9
    Hscan = df["H"].unique()
    amp_diagram, phase_diagram = [], []
    for H in Hscan:
        sub_h = df.loc[df["H"] == H]
        max_value, max_freq = find_max_f_frequency(sub_h["f"].values,
                                                   sub_h["Vmix"].values, cf)

        max_value_phase, max_freq_phrase = find_max_f_frequency(
            sub_h["f"].values, sub_h["phase"].values, 2 * cf)
        if abs(max_freq - cf):
            print("AMP DIFF", max_freq)
        if abs(max_freq_phrase - 2 * cf):
            print("PHASE DIFF", max_freq_phrase)

        amp_diagram.append(max_value)
        phase_diagram.append(max_value_phase)
    amp_diagram = np.asarray(amp_diagram)
    phase_diagram = np.asarray(phase_diagram)
else:
    amp_diagram = df["Vmix"]
    phase_diagram = df["phase"]
    Hscan = df["H"] / 1e3

amp_diagram -= amp_diagram[0]
print(amp_diagram[0])
normalise = False
v = "Hy"

OeToAm = 79.57747
data = pd.read_csv("./Ta_CoFeB_harmonics-new/5360_1f_hx_hy.csv", sep=",")
exp_data1f_H = data[v] * OeToAm / 1e3
exp_data1f_amp = data[f"{v}_1f_voltage"]
exp_data1f_amp -= exp_data1f_amp[0]

data = pd.read_csv("./Ta_CoFeB_harmonics-new/5360_2f_hx_hy.csv", sep=",")
exp_data2f_H = data[v] * OeToAm / 1e3
exp_data2f_amp = data[f"{v}_2f_voltage"]
exp_data2f_amp -= exp_data2f_amp[0]
# data = pd.read_csv(f'./torque-data/5360_66_1f_{v}.txt', sep='\t')
# data = data.stack().str.replace(',', '.').unstack()
# data['Field'] = data['Field'].astype(np.float32) / 1e3
# data['AC Voltage'] = data['AC Voltage'].astype(np.float32)

# exp_data1f_H = data['Field']
# exp_data1f_amp = data['AC Voltage'] - data['AC Voltage'][0]
# # 2f
# data = pd.read_csv(f'./torque-data/5360_66_2f_{v}.txt', sep='\t')
# data = data.stack().str.replace(',', '.').unstack()
# data['Field'] = data['Field'].astype(np.float32) / 1e3
# data['Phase'] = data['Phase'].astype(np.float32)
# exp_data2f_H = data['Field']
# exp_data2f_amp = data['Phase'] - data['Phase'][0]

if normalise:
    print(exp_data1f_amp.max() / amp_diagram.max())
    exp_data1f_amp /= exp_data1f_amp.max()
    exp_data2f_amp /= exp_data2f_amp.max()
    amp_diagram /= abs(amp_diagram.min())
    phase_diagram /= phase_diagram.max()

#
with plt.style.context(["science", "no-latex"]):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
    ax1.plot(Hscan, amp_diagram, "*-", color="r", label="Simulation")
    ax1.plot(exp_data1f_H, exp_data1f_amp, "-.", label="Experiment")
    ax1.set_ylabel("Amplitude [mV]")
    ax1.axhline(y=0, color="g", linestyle="--")
    ax1.legend()
    ax1.set_title("1st Harmonics")

    ax2.plot(Hscan, phase_diagram, "*-", color="r", label="Simulation")
    ax2.plot(exp_data2f_H, exp_data2f_amp, "-.", label="Experiment")
    ax2.legend()
    ax2.set_xlabel("H [kA/m]")
    ax2.set_ylabel("Phase signal [mV]")
    ax2.set_title("2nd Harmonics")
    fig.savefig("Spectra_1f_2f.png")
print("Done plotting")
