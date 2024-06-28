import numpy as np
from scipy.linalg import cholesky, eigh
from scipy.signal import stft, istft
from scipy.signal.windows import hann


def my_mix_signals(n, DC, method):
    """ """
    M = n.shape[1]  # Number of sensors
    K = (DC.shape[2] - 1) * 2  # Number of frequency bins

    # Compute short-time Fourier transform (STFT) of all input signals
    n = np.vstack([np.zeros((K // 2, M)), n, np.zeros((K // 2, M))])
    f, t, N = stft(n.T, window=hann(K), noverlap=int(0.75 * K), nperseg=K, nfft=K, return_onesided=False)

    """
    DC = (M, M, K/2+1)
    N = (M, K, L)
    """

    # Generate output signal in the STFT domain for each frequency bin k
    X = np.zeros_like(N, dtype=complex)
    for k in range(1, K // 2 + 1):
        if method.lower() == 'cholesky':
            C = cholesky(DC[:, :, k])
        elif method.lower() == 'eigen':
            D, V = eigh(DC[:, :, k])
            V = np.asmatrix(V)
            C = np.sqrt(D) @ V.getH()
        else:
            raise ValueError('Unknown method specified.')

        X[:, k, :] = np.squeeze(N[:, k, :]) @ np.conj(C)

    X[K // 2 + 1:, :, :] = np.conj(X[1:K // 2, :, :][::-1])

    # Compute inverse STFT
    _, x = istft(X, window=hann(K), noverlap=int(0.75 * K), nperseg=K, nfft=K, input_onesided=False)
    x = x[:, K // 2: -K // 2].T

    return x
