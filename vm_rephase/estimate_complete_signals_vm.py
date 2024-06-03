import matlab.engine

eng = matlab.engine.start_matlab()
from helpers.spatial_coherence import spatial_coherence


def main():
    eng.eval("run('matlab_files/quickload.m')", nargout=0)

    spatial_coherence("matlab_files/export")

    return


if __name__ == '__main__':
    main()
