from sigpyproc.Readers import FilReader
import numpy as np
from rfipip import rfiUtils


class RfiFil(object):
    def __init__(self, path):
        self.path = path
        self.file = None
        self.fs = None
        self.freqs = None
        self.time = None
        self.time_series = None
        self.bandpass = None

    def _open_file(self):
        """
        Open filterbank file
        :return:
        """
        self.file = FilReader(self.path)
        self.fs = 1.0 / self.file.header.tsamp
        self._create_freqs()

    def _create_freqs(self):
        """
        Create vector with frequencies for each channel
        :return:
        """
        f0 = self.file.header['fch1']
        f_delt = self.file.header['foff']

        i_start, i_stop = 0, self.file.header['nchans']

        # calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx = np.int(i_stop)

        # create freq array
        if i_start < i_stop:
            i_vals = np.arange(chan_start_idx, chan_stop_idx)
        else:
            i_vals = np.arange(chan_stop_idx, chan_start_idx)

        self.freqs = f_delt * i_vals + f0

        # In filterbank channel 0 is the highest frequency
        # so invert
        if f_delt < 0:
            self.freqs = self.freqs[::-1]

    def read_file(self,
                  start_time=0.0,
                  duration=0.0):
        """

        :param start_time: in seconds
        :param duration: in seconds
        :return:
        """
        start_sample = start_time / self.file.header.tsamp
        if duration == 0:
            nsamples = self.file.header.nsamples - start_sample
        else:
            nsamples = duration / self.file.header.tsamp
        block = self.file.readBlock(start_sample, nsamples)
        self._create_time(start_time, duration, block)
        return block, nsamples

        # # TODO change to channel number
        # self.data = self.file.read_filterbank(f_start=f_start, f_stop=f_stop)
        # self.freqs = self.file.freqs

    def _create_time(self,
                     start_time,
                     duration,
                     block):
        """

        :param start_time:
        :param duration:
        :param block:
        :return:
        """  # TODO calc MJD for plotting and analysis
        self.time = np.linspace(start_time,
                                start_time + duration,
                                num=block.shape[1])

    def time_vector(self, vec_length):
        """

        :param vec_length: int
        :return: vecotr with start times and duration of a block 
        """
        start_vector = np.linspace(0,
                                   self.file.header.tobs,
                                   num=vec_length,
                                   endpoint=False,
                                   retstep=True)
        return start_vector[0], start_vector[1]

    def rfi_block(self, start_time, duration):
        """

        :param start_time: 
        :param duration: 
        :return: 
        """
        block, _ = self.read_time_freq(start_time, duration)
        return rfiUtils.rfi_threshold(block)

    def read_time_freq(self,
                       start_time,
                       duration):
        """

        :param start_time:
        :param duration:
        :return:
        """
        block, nsamples = self.read_file(start_time, duration)
        return block, nsamples

    def read_time_series(self):
        """

        :return:
        """
        raw = self.file.collapse()  # collapse freq channels
        self.time_series = raw

    def read_bandpass(self):
        """

        :param data:
        :return:
        """
        self.bandpass = self.file.bandpass()

