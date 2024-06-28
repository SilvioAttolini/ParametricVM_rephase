import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sinc, jv
from scipy.signal import stft
from mix_signals2 import mix_signals2


def version2():
    """ """

    # Initialization
    Fs = 8000  # Sample frequency (Hz)
    c = 340  # Sound velocity (m/s)
    K = 256  # FFT length
    M = 2  # Number of sensors
    d = 0.2  # Inter sensor distance (m)
    type_nf = 'spherical'  # Type of noise field: 'spherical' or 'cylindrical'
    L = 10 * Fs  # Data length

    # Generate M mutually independent input signals of length L
    np.random.seed(1)
    n = np.random.randn(L, M)

    # Generate matrix with desired spatial coherence
    ww = 2 * np.pi * Fs * np.arange(K // 2 + 1) / K
    DC = np.zeros((M, M, K // 2 + 1))
    for p in range(M):
        for q in range(M):
            if p == q:
                DC[p, q, :] = np.ones(K // 2 + 1)
            else:
                if type_nf.lower() == 'spherical':
                    DC[p, q, :] = sinc(ww * abs(p - q) * d / (c * np.pi))
                elif type_nf.lower() == 'cylindrical':
                    DC[p, q, :] = jv(0, ww * abs(p - q) * d / c)
                else:
                    raise ValueError('Unknown noise field.')

    # Generate sensor signals with desired spatial coherence
    x = mix_signals2(n, DC, 'eigen')

    # Compare desired and generated coherence
    K_eval = 256
    ww = 2 * np.pi * Fs * np.arange(K_eval // 2 + 1) / K_eval
    sc_theory = np.zeros((M - 1, K // 2 + 1))
    sc_generated = np.zeros((M - 1, K // 2 + 1))

    # Calculate STFT and PSD of all output signals
    f, t, X = stft(x, fs=Fs, window='hann', nperseg=K_eval, noverlap=int(0.75 * K_eval), nfft=K_eval)
    X = X[:K_eval // 2 + 1, :, :]
    phi_x = np.mean(np.abs(X) ** 2, axis=1)

    # Calculate spatial coherence of desired and generated signals
    for m in range(M - 1):
        if type_nf.lower() == 'spherical':
            sc_theory[m, :] = sinc(ww * (m + 1) * d / (c * np.pi))
        elif type_nf.lower() == 'cylindrical':
            sc_theory[m, :] = jv(0, ww * (m + 1) * d / c)

        # Compute cross-PSD of x_1 and x_(m+1)
        psi_x = np.mean(X[:, :, 0] * np.conj(X[:, :, m + 1]), axis=1)

        # Compute real-part of complex coherence between x_1 and x_(m+1)
        sc_generated[m, :] = np.real(psi_x / np.sqrt(phi_x[:, 0] * phi_x[:, m + 1]))

    # Calculate normalized mean square error
    NMSE = np.zeros(M)
    for m in range(M - 1):
        NMSE[m] = 10 * np.log10(np.sum((sc_theory[m, :] - sc_generated[m, :]) ** 2) / np.sum(sc_theory[m, :] ** 2))

    # Plot spatial coherence of two sensor pairs
    fig, axs = plt.subplots(min(2, M - 1), 1, figsize=(10, 8))
    Freqs = np.arange(0, Fs / 2 + Fs / K, Fs / K)
    for m in range(min(2, M - 1)):
        ax = axs[m] if min(2, M - 1) > 1 else axs
        ax.plot(Freqs / 1000, sc_theory[m, :], '-k', linewidth=1.5, label='Theory')
        ax.plot(Freqs / 1000, sc_generated[m, :], '-.b', linewidth=1.5,
                label=f'Proposed Method (NMSE = {NMSE[m]:.1f} dB)')
        ax.set_xlabel('Frequency [kHz]')
        ax.set_ylabel('Real(Spatial Coherence)')
        ax.set_title(f'Inter sensor distance {(m + 1) * d:.2f} m')
        ax.legend()
        ax.grid(True)

    plt.tight_layout()
    plt.show()

    return
