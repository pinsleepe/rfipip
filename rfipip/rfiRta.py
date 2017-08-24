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
        self.data = None
        self.mode = None

    def _open_file(self):
        """
        Open RTA h5 file
        :return:
        """
        self.file = h5py.File(self.path, mode='r+')



    def read_file(self):
        zero_data = self.file['spectra'][:]
        # strip data and choose frequency channels
        self.data = self._strip_zeros(zero_data)
        # [:, ch_start:ch_stop]
        # assuming that mode doesnt change in observation
        self.mode = self.file['mode'][:][self.time_vector][0]
        self._create_freqs()

