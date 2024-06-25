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

    """


if __name__ == "__main__":
    main()