import matlab.engine
eng = matlab.engine.start_matlab()
from fast_SC_test import fast_SC_test

def main() -> None:
    """
    I
        Noise field assumptions:
            1) Spatially homogeneous: the physical properties of the sound field do not depend on the absolute
                    positions of the sensors
            2) Isotropic: the physical properties of the sound field are the same in all directions
            3) Time invariant

        Algorithm overview:
            1) Generate M mutually independent noise signals (M = num of sensors)
            2) Filter and mix the noise signals in order to obtain the specified SC (spatial coherence)
            3) Get M sensors signals with the specified SC

    II                                     | 0 x1 ... xM |
        Define the sensors' positions: P = | 0 y1 ... yM |, with P0 as reference
                                           | 0 z1 ... zM |

        Let p, q in {1, ... , M}
        d_pq = distance(Pp, Pq)

        Let w = omega = 2 * pi * f

        Power Spectral Density (PSD):
            - AUTO: PHI_pp(w) = scipy.signal.welch(p, Fs, ...)
            - CROSS: PHI_pq(w) = scipy.signal.csd(p, q, Fs, ...)
            ref.: https://www.researchgate.net/profile/Fernando-Schlindwein/post/Analyzing-EEGs-with-FFTs-and-PSDs/attachment/59d64f8879197b80779a8a73/AS%3A498242888896512%401495801769748/download/Welch_The+use+of+fast+Fourier+transform+for+the+estimation+of+power+spectra_A+method+based+on+time+averaging+over+short%2C+modified+periodograms.pdf

        The homogeneity of the SF is represented by: PHI_pp(w) = PHI_qq(w) = PHI(w)

        Spatial Coherence definition:
            gamma_pq(w) = PHI_pq(w) / sqrt(PHI_pp(w)*PHI_qq(w))

        In 3D the expected (ground truth) SC is:
            let c = 343 m/s
            let k = w * d_pq * c
            Gamma_pq(w) = sinc(k) = sin(k) / k

        The STFT domain: (l, w_k), l time index, w_k in {1, ... K/2 -1} with K = frame length

                                       | Gamma_11(w) ... Gamma_1M(w) |
        Define the matrix GAMMA(w_k) = |    ...              ...     |  with SC ground truth
                                       | Gamma_M1(w) ... Gamma_MM(w) |
            in which:
                - we expect Gamma_pq(w) = gamma_pq(w)
                - we expect PSDp = PSDq

        Let the sensors signals: X(l, w_k) = [X_1(l, w_k), ..., X_M(l, w_k)].T
        Let the noise signals: N(l, w_k) = [N_1(l, w_k), ..., N_M(l, w_k)].T
        Let the mixing matrix C(w_k)

        In the STFT domain we calculate the instantaneous mixing coefficients from:
            X(l, w_k) = C(w_k).H * N(l, w_k)

    III
        Generate the M mutually independent noise signals, based on the homogeneity of the SF
        Assuming that the columns of the mixing matrix C have all equal norm
        If the short-term (aka inside a specific (l, w_k) frame) PSDs of the noise signals are
        equal, then the short-term PSDs of the sensors signals are equal.

        A) Perfectly homogeneous NF

            Let D_p(l, w_k) ~ U[-1, 1]
            Let PHI(l, w_k) be the same PSD for all the sensors

            The STFT coefficients of the p-th noise signal are:
                N_p(l, w_k) = sqrt(PHI(l, w_k)) * exp(i * pi * D_p(l, w_k))

            PHI(l, w_k) = PHI(w_k) when dealing with speech or factory noise, while with babble-speech only
            the short-term PSDs can coincide and the overall phase spectrum is destructed

        B) Approximately homogeneous NF

            Generate M mutually independent signals that consist of continuous (without periods of silence)
            babble speech or factory noise. These noise signals must have the same power: normalize them all with
            respect to the first one. In the STFT domain we do not expect the short-term PSDs of N_p(l, w_k)
            to be equal to that of N_q(l, w_k), therefore we can accept small fluctuations in the long-term PSDs.
            This represents the approximate homogeneity. Also, the signals resulting from the mixing operations
            will inevitably sound like a mixture of babble speeches or factory noises.

    IV
        Determine the mixing matrix C(w_k), which is responsible for imposing the desired SC.

        Prerequisites on C(w_k):
            1) the inner product of the columns p and q must give back Gamma_pq(w_k). E.g.: in a 2x2 C(w_k)
                    matrix, Cp · Cq = x_p*x_q + y_p*y_q
            2) the norm of each column vector of C(w_k) must be equal to 1. norm(C_p) = 1
            3) The SC between p and q must be Real and equal to Gamma_pq(w_k). C_p · C_q = Gamma_pq(w_k)

        Techniques to build C(w_k):
            1) using RIRs to obtain noise sources. This my_implementation is too computationally expensive and the SC
                    created depend on the physical positions of the sensors.
            2) using Cholesky decomposition: GAMMA(w_k) = C(w_k).H * C(w_k). This my_implementation creates an upper-right
                    triangular matrix, which creates unnatural results, since X_1 would depend only on N_1, while
                    X_M would depend on all the previous Ns. Also, this resulting matrix is only suitable
                    to omnidirectional sources.
            3) The general solution is given by the eigenvalue decomposition (EVD) of the matrix GAMMA(w_k):
                    GAMMA(w_k) = V(w_k) * D(w_k) * V(w_k).H
                               = V(w_k) * sqrt(D(w_k)) * sqrt(D(w_k)) * V(w_k).H
                               = V(w_k) * sqrt(D(w_k)) * C(w_k)
                   Even in this my_implementation all the noises to not contribute equally, but the result is
                   sufficiently satisfying.
                   In the single directional source, the EVD will build only the first row of C(w_k) with
                   positive values.

    V
        Algorithm summary:
            1) Define GAMMA(w_k) with the seeked SCs
            2) Calculate the EVD of GAMMA(w_k) to obtain C(w_k) = sqrt(D(w_k)) * V(w_k).H
            3) Generate the M mutually independent random signals N_p(l, w_k)
            4) For all sensors p
                    For all l in L
                            For all w_k in K
                                    X_p(l, w_k) = C(w_k).H * N_p(l, w_k)
            5) x_p(t) = istft(X_p(l, w_k))

            overall complexity = O(LKM^2 * log_2(K))

    VI
        Use more efficient techniques for the EDV
        Exploit the spectrum being conjugate symmetric
        Babble_speech_noise = Total_Signal - Direct_Signal
        Use K = 256
        Use d_pq = 0.2 m

    """

    # My implementation
    eng.eval("run('my_implementation/habets.m')", nargout=0)
    #fast_SC_test()

    # Original paper
    # eng.eval("run('ANF-Generator/gen_noisefield.m')", nargout=0)
    # eng.eval("run('ANF-Generator/gen_babble_speech.m')", nargout=0)


if __name__ == "__main__":
    main()

"""
~ = AltGr + ì
· = AltGr + .
"""
