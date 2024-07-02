import numpy as np
from scipy.signal import stft, istft
from scipy.signal.windows import hann


def mix_signals3(n, DC, method):
    M = n.shape[1]  # Number of sensors
    K = (DC.shape[2] - 1) * 2  # Number of frequency bins

    # Compute short-time Fourier transform (STFT) of all input signals
    n = np.concatenate((np.zeros((K // 2, M)), n, np.zeros((K // 2, M))))
    N = stft(n, window=hann(K), noverlap=int(0.75 * K), nperseg=K)  #, centered=False)

    # Generate output signal in the STFT domain for each frequency bin k
    X = np.zeros_like(N, dtype=complex)
    for k in range(2, K // 2 + 1):
        if method.lower() == 'cholesky':
            C = np.linalg.cholesky(DC[:, :, k])
        elif method.lower() == 'eigen':
            V, D = np.linalg.eig(DC[:, :, k])
            C = np.sqrt(D) @ V.T
        else:
            raise ValueError('Unknown method specified.')

    X[k, :, :] = np.squeeze(N[k, :, :]) @ np.conj(C)

    X[K // 2 + 2:, :, :] = np.conj(X[K // 2 + 2::-1, :, :])

    # Compute inverse STFT
    x = istft(X, window=hann(K), noverlap=int(0.75 * K), nperseg=K)  # , centered=False, conjugatesymmetric=True)
    x = x[K // 2 + 1: -K // 2, :]

    return x
