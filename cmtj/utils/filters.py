from scipy.signal import butter, lfilter


class Filters:
    @staticmethod
    def butter_bandpass_filter(data, pass_freq, fs, order=5):
        nyq = 0.5 * fs
        if isinstance(pass_freq, float):
            if pass_freq == 0:
                pass_freq = 0.1
                try:
                    b, a = butter(order, [0.9 * pass_freq / nyq, pass_freq / nyq],
                                btype='bandpass',
                                analog=False)
                except ValueError:
                    print(fs, pass_freq, nyq, 0.9 * pass_freq / nyq, pass_freq / nyq)
                    raise ValueError("Error in filtering")
        elif isinstance(pass_freq, tuple):
             b, a = butter(order, [pass_freq[0], pass_freq[1]],
                                btype='bandpass',
                                analog=False)
        y = lfilter(b, a, data, zi=None)
        return y

    @staticmethod
    def butter_lowpass_filter(data, cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = lfilter(b, a, data, zi=None)
        return y