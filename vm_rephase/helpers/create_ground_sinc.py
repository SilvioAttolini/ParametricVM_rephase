import numpy as np

from helpers.plots import plot_sinc

"""
" unitary ground sinc
"""


def create_ground_sinc(freqs_gamma, params):
    c = 340
    d = params['d']
    k = 2 * np.pi * d / c

    if freqs_gamma[0] > 0:
        sinc = np.sin(k * freqs_gamma) / (k * freqs_gamma)
    else:
        sinc = np.sin(k * freqs_gamma[1:]) / (k * freqs_gamma[1:])
        sinc = np.insert(sinc, 0, 1)

    plot_sinc(freqs_gamma, sinc, False)

    return sinc
