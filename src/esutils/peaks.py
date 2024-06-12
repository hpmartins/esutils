import numpy as np
import scipy.signal as signal

def highest_peak(x, y):
    i_peaks, _ = signal.find_peaks(y)
    i_max_peak = i_peaks[np.argmax(y[i_peaks])]
    return x[i_max_peak], y[i_max_peak]
