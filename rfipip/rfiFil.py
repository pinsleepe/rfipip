from sigpyproc.Readers import FilReader
import numpy as np
from rfipip import rfiUtils
from astropy.time import Time
# from datetime import datetime

# no info about antennas if fil files ??


class RfiFil(object):
    def __init__(self, path,
                 observer):
        """

        :param path:
        :param observer: ephem.Observer() object
        """
        self.path = path
        self.file = None
        self.observer = observer
        self.metadata = {}

        # self.freqs_vector = None
        # self.time_vector = None
        # self.time_series = None
        # self.bandpass = None

    def init(self):
        self.open_file()
        self.fill_metadata()

    def fill_metadata(self):
        """
        Fill the metadata with info from H5
        :return:
        """
        self.freq_metadata()
        self.other_metadata()

    def fill_observer(self):
        """

        :return:
        """
        self.observer.horizon = np.radians(15)

    def open_file(self):
        """
        Open filterbank file
        :return:
        """
        self.file = FilReader(self.path)

    def freq_metadata(self):
        """
        Frequency specific metadata
        :return:
        """
        self.metadata['fs'] = 1.0 / self.file.header.tsamp
        self.metadata['bandwidth'] = self.file.header.bandwidth * 10e5
        self.metadata['f0'] = self.file.header['fch1']
        self.metadata['f_delt'] = self.file.header['foff']
        self.metadata['nchannels'] = self.file.header['nchans']
        self.metadata['cenfreq'] = self.file.header['fcenter'] + self.metadata['f_delt']/2

    def other_metadata(self):
        self.metadata['nsamples'] = self.file.header.nsamples
        self.metadata['tsamp'] = self.file.header['tsamp']
        # time stamp (MJD) of first sample
        self.metadata['tstart'] = self.file.header['tstart']

    def target_metadata(self):
        """
        antenna specific metadata
        :return:
        """
        ant = 'antenna'  # as we don't know the antenna
        self.metadata[ant] = {}
        self.metadata[ant]['target'] = self.file.header.source_name
        # TODO has to monitor those numbers with future updates
        self.metadata[ant]['az'] = self.file.header.src_raj
        self.metadata[ant]['el'] = self.file.header.src_dej
        t = Time(self.metadata['tstart'], format='mjd')
        self.observer.date = t.datetime
        ra, dec = self.observer.radec_of(np.radians(self.metadata[ant]['az']),
                                         np.radians(self.metadata[ant]['el']))
        self.metadata[ant]['ra'] = str(ra)
        self.metadata[ant]['dec'] = str(dec)

    def create_freqs(self):
        """
        Create vector with centre frequencies for each channel
        :return:
        """
        i_start, i_stop = 0, self.metadata['nchannels']

        # calculate closest true index value
        chan_start_idx = np.int(i_start)
        chan_stop_idx = np.int(i_stop)

        # create freq array
        if i_start < i_stop:
            i_vals = np.arange(chan_start_idx, chan_stop_idx)
        else:
            i_vals = np.arange(chan_stop_idx, chan_start_idx)

        freqs_vector = self.metadata['f_delt'] * i_vals + self.metadata['f0']

        # In filterbank channel 0 is the highest frequency
        # so invert
        if self.metadata['f_delt'] < 0:
            freqs_vector = freqs_vector[::-1]
        return freqs_vector

    def read_file(self,
                  start_time=0.0,
                  duration=0.0):
        """

        :param start_time: in seconds
        :param duration: in seconds
        :return:
        """
        start_sample = start_time / self.metadata['tsamp']
        if duration == 0:
            nsamples = self.metadata['nsamples'] - start_sample
        else:
            nsamples = duration / self.metadata['tsamp']
        block = self.file.readBlock(start_sample, nsamples)
        # self._create_time(start_time, duration, block)
        return block, nsamples

        # # TODO change to channel number
        # self.data = self.file.read_filterbank(f_start=f_start, f_stop=f_stop)
        # self.freqs_vector = self.file.freqs_vector

    # def _create_time(self,
    #                  start_time,
    #                  duration,
    #                  block):
    #     """
    #
    #     :param start_time:
    #     :param duration:
    #     :param block:
    #     :return:
    #     """  # TODO calc MJD for plotting and analysis
    #     self.time_vector = np.linspace(start_time,
    #                                    start_time + duration,
    #                                    num=block.shape[1])
    #
    # def time_vector(self, vec_length):
    #     """
    #
    #     :param vec_length: int
    #     :return: vecotr with start times and duration of a block
    #     """
    #     start_vector = np.linspace(0,
    #                                self.file.header.tobs,
    #                                num=vec_length,
    #                                endpoint=False,
    #                                retstep=True)
    #     return start_vector[0], start_vector[1]
    #
    # def rfi_block(self, start_time, duration):
    #     """
    #
    #     :param start_time:
    #     :param duration:
    #     :return:
    #     """
    #     block, _ = self.read_time_freq(start_time, duration)
    #     return rfiUtils.rfi_threshold(block)
    #
    # def read_time_freq(self,
    #                    start_time,
    #                    duration):
    #     """
    #
    #     :param start_time:
    #     :param duration:
    #     :return:
    #     """
    #     block, nsamples = self.read_file(start_time, duration)
    #     return block, nsamples
    #
    # def read_time_series(self):
    #     """
    #
    #     :return:
    #     """
    #     raw = self.file.collapse()  # collapse freq channels
    #     self.time_series = raw
    #
    # def read_bandpass(self):
    #     """
    #
    #     :param data:
    #     :return:
    #     """
    #     self.bandpass = self.file.bandpass()
    #
