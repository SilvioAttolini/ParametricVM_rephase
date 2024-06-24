import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, csd, coherence


def spat_cohers_2(complete_vms_time, params, n_Vms):
    # Parameters for the cross-PSD calculation
    window = 'hamming'  # Window function (e.g., Hamming window)
    nperseg = 1024  # Length of each segment
    noverlap = 512  # Number of overlapping samples
    nfft = 1024  # Number of FFT points
    fs = params['Fs']

    micA = 0
    micB = micA + 1
    i = 0
    while micB < n_Vms:
        print(f"pair {micA + 1}-{micB + 1}")

        # create current signals
        x = complete_vms_time[micA, :]
        y = complete_vms_time[micB, :]

        # Calculate the PSD using the Welch method
        # f, pxx = welch(x, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)
        # _, pyy = welch(y, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)

        # Plot the PSD
        # plt.figure()
        # plt.semilogy(f, pxx)  # Use semilogarithmic scale for better visualization
        # plt.xlabel('Frequency (Hz)')
        # plt.ylabel('Power/Frequency (dB/Hz)')
        # plt.title('Power Spectral Density')
        # plt.grid(True)
        # plt.show()

        # Calculate the cross-PSD using the Welch method
        f, Pxy = csd(x, y, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)

        # Plot the cross-PSD
        # plt.figure()
        # plt.semilogy(f, np.abs(Pxy))  # Use semilogarithmic scale for better visualization
        # plt.xlabel('Frequency (Hz)')
        # plt.ylabel('Cross Power/Frequency (dB/Hz)')
        # plt.title('Cross Power Spectral Density')
        # plt.grid(True)
        # plt.show()

        # Calculate the coherence
        Cxy = coherence(x, y, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)[1]

        # Plot the coherence
        plt.figure()
        plt.plot(f, Cxy)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Coherence')
        plt.title('Coherence between x and y')
        plt.grid(True)
        plt.show()

        # loop
        micA = micB + 1
        micB = micA + 1
        i = i + 1

    return
