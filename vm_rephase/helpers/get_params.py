import os

"""
" Create parameters and folders
"""


def get_params(save_dir):
    # create params
    params = {'Fs': 16000,
              'winLength': 4096,
              'sigLenSecs': 5,
              'd': 0.078520,
              'save_directory': save_dir
              }
    save_directory = params['save_directory']

    save_stereo = f'{save_directory}/stereo'
    os.makedirs(save_stereo, exist_ok=True)
    params['save_stereo'] = save_stereo

    save_signal_time = f'{save_directory}/signal_time'
    os.makedirs(save_signal_time, exist_ok=True)
    params['save_signal_time'] = save_signal_time

    save_signal_stft = f'{save_directory}/signal_stft'
    os.makedirs(save_signal_stft, exist_ok=True)
    params['save_signal_stft'] = save_signal_stft

    save_phases = f'{save_directory}/phases'
    os.makedirs(save_phases, exist_ok=True)
    params['save_phases'] = save_phases

    save_spatial_coherence = f'{save_directory}/spatial_coherence'
    os.makedirs(save_spatial_coherence, exist_ok=True)
    params['save_spatial_coherence'] = save_spatial_coherence

    save_vm_audio = f'{save_directory}/vm_audio'
    os.makedirs(save_vm_audio, exist_ok=True)
    params['save_vm_audio'] = save_vm_audio

    return params
