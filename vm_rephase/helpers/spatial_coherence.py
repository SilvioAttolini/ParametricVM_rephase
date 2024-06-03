import h5py
from scipy.io.wavfile import write
from scipy.signal import stft

from helpers.auto_and_cross_spectrum import auto_and_cross_power_spectrum
from helpers.create_ground_sinc import create_ground_sinc
from helpers.get_params import get_params
from helpers.plots import *


def spatial_coherence(reading_folder):
    print("Rebuilding data...")

    # parameters
    params = get_params()
    save_vm_audio = params['save_vm_audio']
    save_phases = params['save_phases']
    save_signal_time = params['save_signal_time']
    save_signal_stft = params['save_signal_stft']
    save_spatial_coherence = params['save_spatial_coherence']
    save_stereo = params['save_stereo']

    # import the time domain complete signals
    with h5py.File(f'{reading_folder}/completeEstimate.mat', 'r') as file:
        complete_vms_time = np.array(file[f'completeEstimate'])
    n_Vms, _ = np.shape(complete_vms_time)

    # get axis dimensions
    fAx_stft, tAx_stft, _ = stft(complete_vms_time[0, :], fs=params['Fs'], window='hann', nperseg=256, noverlap=128)

    # create ground sinc
    sinc = create_ground_sinc(fAx_stft, params)

    # read signals
    complete_vms_stft = np.zeros((n_Vms, len(fAx_stft), len(tAx_stft)), dtype='complex')  # 10x128x689
    for vm in range(n_Vms):
        print(f"vm: {vm}")

        imported_rir_time = complete_vms_time[vm, :]

        # create output file
        signal_int16 = np.int16(imported_rir_time / np.max(np.abs(imported_rir_time)) * 32767)
        filename = f'{save_vm_audio}/complete_vm{vm + 1}.wav'
        write(filename, params['Fs'], signal_int16)

        # Calculate the STFT of the diffuse contribution of the vm rir
        fAx_stft, tAx_stft, stft_of_imported = stft(
            imported_rir_time, fs=params['Fs'], window='hann', nperseg=256, noverlap=128)

        complete_vms_stft[vm, :, :] = stft_of_imported

        # phase
        phase = np.angle(stft_of_imported)
        t = 0
        ph = phase[:, t]  # all freqs at time zero
        plot_complex(fAx_stft, np.real(ph), np.imag(ph), "reference", t, vm, save_phases, False)
        t = 500
        ph = phase[:, t]  # all freqs at time 500
        plot_complex(fAx_stft, np.real(ph), np.imag(ph), "reference", t, vm, save_phases, False)

        # plot
        plot_signal_time(params, imported_rir_time, vm, save_signal_time, False)
        plot_signal_stft(tAx_stft, fAx_stft, complete_vms_stft[vm, :, :], vm, save_signal_stft, False)

    print("Computing Spatial Coherence...")

    micA = 0
    micB = micA + 1
    i = 0
    real_spat_coh_of_all_pairs = np.zeros((n_Vms // 2, len(fAx_stft)))
    while micB < n_Vms:
        print(f"pair {micA + 1}-{micB + 1}")

        # create current signals
        A = complete_vms_stft[micA, :, :]
        B = complete_vms_stft[micB, :, :]

        # compute the complex spatial coherence
        phi_na_nb = auto_and_cross_power_spectrum(A, B)
        phi_n_n = auto_and_cross_power_spectrum(A, A)  # ASSUMPTION: IT IS EQUAL TO (B,B)
        Gamma_nab = phi_na_nb / phi_n_n

        # keep just the real part to plot
        real_Gamma_nab = np.real(Gamma_nab)

        plot_spatial_coherence(fAx_stft, real_Gamma_nab, sinc, micA, micB, save_spatial_coherence, False)

        real_spat_coh_of_all_pairs[i, :] = real_Gamma_nab

        # stereo audio of this pair
        signal_left = np.int16(complete_vms_time[micA, :] / np.max(np.abs(complete_vms_time[micA, :])) * 32767)
        signal_right = np.int16(complete_vms_time[micB, :] / np.max(np.abs(complete_vms_time[micB, :])) * 32767)
        stereo_signal = np.column_stack((signal_left, signal_right))
        write(f"{save_stereo}/stereo_signal_vms{micA + 1}_{micB + 1}", params['Fs'], stereo_signal)

        # loop
        micA = micB + 1
        micB = micA + 1
        i = i + 1

    # mean of the spatial coherence over all the pairs
    spatial_coherence_avg = np.mean(real_spat_coh_of_all_pairs, axis=0)

    # plot average spatial coherence
    plot_average_spatial_coherence(fAx_stft, spatial_coherence_avg, sinc, n_Vms, save_spatial_coherence, False)

    print("Done.")

    return
