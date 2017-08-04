import h5py

# Whole band
# TODO add channel bw
rta_modes = {'1': {'n_chan': 32768,
                   'top_freq': 899972534.18,
                   'base_freq': 0.0},
             '2': {'n_chan': 32768,
                   'top_freq': 1199981689.45,
                   'base_freq': 600000000.0},
             '3': {'n_chan': 32768,
                   'top_freq': 1709973907.47,
                   'base_freq': 855000000.0},
             '4': {'n_chan': 32768,
                   'top_freq': 2699972534.18,
                   'base_freq': 1800000000.0}}


class RfiRta(object):
    def __init__(self, path):
        self.path = path

    def _open_file(self):
        """
        Open RTA h5 file
        :return:
        """
        self.file = h5py.File(self.path, mode='r+')
