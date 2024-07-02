import numpy as np
from scipy.signal import stft, istft
from scipy.linalg import cholesky


def mix_signals2(n, DC, method):
    M = n.shape[1]  # Number of sensors
    K = (DC.shape[2] - 1) * 2  # Number of frequency bins

    # Compute short-time Fourier transform (STFT) of all input signals
    n = np.vstack((np.zeros((K // 2, M)), n, np.zeros((K // 2, M))))
    f, t, N = stft(n, window='hann', nperseg=K, noverlap=0.75 * K, nfft=K, return_onesided=False)

    # Generate output signal in the STFT domain for each frequency bin k
    X = np.zeros_like(N)
    for k in range(1, K // 2 + 1):
        if method.lower() == 'cholesky':
            C = cholesky(DC[:, :, k])
        elif method.lower() == 'eigen':
            D, V = np.linalg.eig(DC[:, :, k])
            C = np.sqrt(np.diag(D)) @ V.T
        else:
            raise ValueError('Unknown method specified.')

        X[k, :, :] = N[k, :, :] @ np.conj(C).T

    X[K // 2 + 1:, :, :] = np.conj(X[K // 2 - 1:0:-1, :, :])

    # Compute inverse STFT
    _, x = istft(X, window='hann', nperseg=K, noverlap=0.75 * K, nfft=K, input_onesided=False)
    x = x[K // 2:-K // 2, :]

    return x
