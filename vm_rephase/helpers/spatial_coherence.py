import h5py
import numpy as np
from scipy.io.wavfile import write
from scipy.signal import stft
from helpers.spat_coher_1 import spat_coher_1
from helpers.spat_cohers_2 import spat_cohers_2
from helpers.auto_and_cross_spectrum import auto_and_cross_power_spectrum
from helpers.create_ground_sinc import create_ground_sinc
from helpers.get_params import get_params
from helpers.plots import *


def spatial_coherence(reading_folder):
    print("Rebuilding data...")

    SAVE_FOLDER = 'results/prephased'
    INPUT_MAT = 'completeEstimate_prephased'
    MAT_NAME = 'completeEstimate'  # completeEstimate  completeReference
    METHOD = "rephased"  # "rephased"  # "reference"

    # parameters
    params = get_params(SAVE_FOLDER)
    save_vm_audio = params['save_vm_audio']
    save_phases = params['save_phases']
    save_signal_time = params['save_signal_time']
    save_signal_stft = params['save_signal_stft']
    save_spatial_coherence = params['save_spatial_coherence']
    save_stereo = params['save_stereo']

    # import the time domain complete signals
    with h5py.File(f'{reading_folder}/{INPUT_MAT}.mat', 'r') as file:
        complete_vms_time = np.array(file[f'{MAT_NAME}'])
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
        plot_complex(fAx_stft, np.real(ph), np.imag(ph), METHOD, t, vm, save_phases, False)
        t = 500
        ph = phase[:, t]  # all freqs at time 500
        plot_complex(fAx_stft, np.real(ph), np.imag(ph), METHOD, t, vm, save_phases, False)

        # plot
        plot_signal_time(params, imported_rir_time, vm, save_signal_time, False)
        plot_signal_stft(tAx_stft, fAx_stft, complete_vms_stft[vm, :, :], vm, save_signal_stft, False)

    print("Computing Spatial Coherence...")

    spat_coher_1(complete_vms_time, complete_vms_stft, n_Vms, fAx_stft, sinc,
                 save_spatial_coherence, save_stereo, params)

    # spat_cohers_2(complete_vms_time, params, n_Vms)

    print("Done.")

    return
