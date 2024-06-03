import numpy as np

"""
" U(f,l) and V(f,l) are the stft of u(t) and v(t) respectively.
" f = freq_bin, l = time_frame
" Each STFT array must have dimensions (freq_bins, time_frames)
"""


def auto_and_cross_power_spectrum(U, V):
    UV_conj_product = U * np.conj(V)

    # Compute the expectation (mean) across the time axis to approximate E[UV*]
    CPS = np.mean(UV_conj_product, axis=-1)

    return CPS
