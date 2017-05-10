#!/usr/bin/env python

# This module uses snippets of code taken from the rfDB2 (https://github.com/ska-sa/rfDB2).
# This is because rfDB2 hasn't been maintained for 3 years now, and many packages it uses are
# also obsolete.

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from rfiUtils import detect_peaks
from skimage import filters
from statistics import median
from rfipip import rfiDatabase, rfiUtils, rfiEvent
import pandas as pd
from scipy.ndimage import measurements

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

training_set = pd.DataFrame(columns=('event',
                                     'c_freq',
                                     'bw',
                                     't_start',
                                     'duration',
                                     'culprit',
                                     'description',
                                     'band'))


class RfiObservation(object):
    def __init__(self,
                 path='',
                 h5_file=False,
                 fil_file=False,
                 rta_file=False):
        """

        :param path:
        :param h5_file:
        :param fil_file:
        :param rta_file:
        """
        self.path = path

        self._h5_file = h5_file
        self._fil_file = fil_file
        self._rta_file = rta_file

        self.bandpass = None
        self.file = None
        self.freqs = None
        self.time = None
        self.data = None

        if self._rta_file:
            self.mode = None
            self.num_channels = None
            self.channel_bw = None

        if self._fil_file:
            self.fs = None
            self.time_series = None

        self.database = None
        self.threshold = None
        self.corrupted_samples = 0
        self.events = None

        # Initialise
        self._init()

    def _init(self):
        """
        Initialise the rfi observation i.e. open and read the data
        :return:
        """
        self._open_file()
        if self._rta_file:
            self.read_file()

    def _open_h5_file(self):
        """
        Open h5 file
        :return:
        """
        if self._h5_file:
            import h5py
            self.file = h5py.File(self.path, mode='r+')

    def _open_rta_file(self):
        """
        Open ratty h5 file
        :return:
        """
        if self._rta_file:
            import h5py
            self.file = h5py.File(self.path, mode='r+')

    def _open_fil_file(self):
        """
        Open filterbank file
        :return:
        """
        if self._fil_file:
            from sigpyproc.Readers import FilReader
            self.file = FilReader(self.path)
            self.fs = 1.0 / self.file.header.tsamp
            self._create_freqs()

    def _strip_zeros(self, data):
        """
        Strip zero columns from data
        :param data:
        :return:
        """
        if self._rta_file:
            # ratty samples every 2 sec so get rid of the zeros
            num_sec = data.shape[0]  # 3600 seconds
            # [0...3599]
            full_sec_idx = range(num_sec)
            # find idx with zeros
            zero_idx = list(np.where(~data.any(axis=1))[0])
            # delete zero idx from full idx list
            data_idx = np.delete(full_sec_idx, zero_idx, axis=0)
            stripped_data = data[data_idx]
            self.time = data_idx
            return stripped_data
        if self._fil_file:
            num_sec = data.shape[1]  # time samples
            full_sec_idx = range(num_sec)
            # find idx with zeros
            zero_idx = list(np.where(~data.any(axis=0))[0])
            # delete zero idx from full idx list
            data_idx = np.delete(full_sec_idx, zero_idx, axis=0)
            stripped_data = data[:, data_idx]
            self.time = data_idx
            return stripped_data

    def _open_file(self):
        """
        Open file according to type
        :return:
        """
        if self._h5_file:
            self._open_h5_file()  # hdf5 object
        if self._fil_file:
            self._open_fil_file()
        if self._rta_file:
            self._open_rta_file()

    def read_file(self, start_time=0, duration=0):
        """

        :param start_time: in seconds
        :param duration: in seconds
        :return:
        """
        if self.file:
            if self._fil_file:
                start_sample = long(start_time / self.file.header.tsamp)
                if duration == 0:
                    nsamples = self.file.header.nsamples - start_sample
                else:
                    nsamples = long(duration / self.file.header.tsamp)
                block = self.file.readBlock(start_sample, nsamples)
                self._create_time(start_time, duration, block)
                return block, nsamples

                # # TODO change to channel number
                # self.data = self.file.read_filterbank(f_start=f_start, f_stop=f_stop)
                # self.freqs = self.file.freqs
            if self._rta_file:
                zero_data = self.file['spectra'][:]
                # strip data and choose frequency channels
                self.data = self._strip_zeros(zero_data)
                # [:, ch_start:ch_stop]
                # assuming that mode doesnt change in observation
                self.mode = self.file['mode'][:][self.time][0]
                self._create_freqs()

    def _create_freqs(self):
        """
        Create vector with frequencies for each channel
        :return:
        """
        if self._rta_file:
            mode = str(self.mode)
            self.freqs = np.linspace(rta_modes[mode]['base_freq'],
                                     rta_modes[mode]['top_freq'],
                                     rta_modes[mode]['n_chan'])
            self.channel_bw = self.freqs[1] - self.freqs[0]
        if self._fil_file:
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

    def _create_time(self,
                     start_time,
                     duration,
                     block):
        """

        :param start_time:
        :param duration:
        :param block:
        :return:
        """
        if self._fil_file:
            # TODO calc MJD for plotting and analysis
            self.time = np.linspace(start_time,
                                    start_time+duration,
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

    def rfi_threshold(self, data_array):
        """
        Calculate threshold based on Yen algorithm
        :param data_array: numpy array
        :return: 
        """
        if self._fil_file:
            return filters.threshold_yen(data_array)

    def rfi_block(self, start_time, duration):
        """
        
        :param start_time: 
        :param duration: 
        :return: 
        """
        block, _ = self.read_time_freq(start_time, duration)
        return self.rfi_threshold(block)

    def rfi_median(self, vec_length):
        """
        Run through entire file and give median value for threshold
        :param vec_length: int, how many block there is in a file,
                            need some logic...
        :return: 
        """
        start_vector, duration = self.time_vector(vec_length)
        obs_val = [self.rfi_block(st, duration) for _, st in enumerate(start_vector)]
        self.threshold = median(obs_val)
        # return median(obs_val)

    def sum_dimension(self,
                      data,
                      axis=0):
        """

        :param data: 2D array
        :param axis:
        :return:
        """
        return data.sum(axis=axis)

    def read_bandpass(self):
        """

        :param data:
        :return:
        """
        if self._rta_file:
            if self.data is not None:
                self.bandpass = self.data.sum(axis=0)
        if self._fil_file:
            self.bandpass = self.file.bandpass()

    def plot_bandpass(self,
                      data=None,
                      ch_start=0,
                      ch_stop=0):
        """
        Works for summing 2D data along Y axis
        :param data:
        :return:
        """
        if self._rta_file:
            fig, ax = plt.subplots()
            ax.plot(self.freqs[ch_start:ch_stop]/1e9, data)
            ax.set_xlabel('Frequency')
            ax.set_ylabel('Power')
            ax.set_aspect('auto')
            plt.show()
        if self._fil_file:
            if self.bandpass is not None:
                bp = self.bandpass  # collapse freq channels
            else:
                self.read_bandpass()
                bp = self.bandpass
            plt.figure()
            ax = plt.subplot(1, 1, 1)
            ax.plot(bp)
            ax.set_xlabel('Freq channel')
            ax.set_title('Bandpass')
            ax.set_aspect('auto')
            plt.show()

    # TODO axis ranges incorrect - fix
    def plot_spectrum(self,
                      data,
                      ch_start,
                      ch_stop):
        """

        :param data:
        :param ch_start:
        :param ch_stop:
        :return:
        """
        fig, ax = plt.subplots()
        if self._rta_file:
            ax.imshow(data,
                      extent=[self.freqs[ch_start] / 1e9,
                              self.freqs[ch_stop] / 1e9,
                              self.time[-1],
                              self.time[0]],
                      cmap=cm.Blues)
            ax.set_xlabel('Frequency [MHz]')
            ax.set_ylabel('Time (seconds)')
        if self._fil_file:
            ax.imshow(data,
                      extent=[self.time[0],
                              self.time[-1],
                              self.freqs[ch_start],
                              self.freqs[ch_stop]],
                      cmap=cm.Blues)
            ax.set_ylabel('Frequency [MHz]')
            ax.set_xlabel('Time (seconds)')

        ax.set_aspect('auto')
        plt.show()

    def plot_mask(self,
                  mask,
                  start_time):
        start_sample = long(start_time / self.file.header.tsamp)

        plt.figure()
        ax = plt.subplot(1, 1, 1)
        plt.imshow(mask, extent=[start_sample * self.file.header.tsamp,
                                 (start_sample + np.shape(mask)[1]) * self.file.header.tsamp,
                                 np.shape(mask)[0], 0])
        ax.set_aspect("auto")
        ax.set_xlabel('observation time (secs)')
        ax.set_ylabel('freq channel')
        ax.set_title('mask time-freq plot')
        ax.set_aspect('auto')
        plt.colorbar()
        plt.show()

    def read_time_freq(self,
                       start_time,
                       duration):
        """

        :param start_time:
        :param duration:
        :return:
        """
        if self._fil_file:
            block, nsamples = self.read_file(start_time, duration)
            return block, nsamples

    def plot_time_freq(self,
                       start_time,
                       duration,
                       zeros=False):
        """

        :param start_time:
        :param duration:
        :return:
        """
        if self._fil_file:
            start_sample = long(start_time / self.file.header.tsamp)
            block, nsamples = self.read_time_freq(start_time, duration)
            if zeros:
                new_block = self._strip_zeros(block)
                block = new_block

            plt.figure()
            ax = plt.subplot(1, 1, 1)
            plt.imshow(block, extent=[start_sample * self.file.header.tsamp,
                                      (start_sample + np.shape(block)[1]) * self.file.header.tsamp,
                                      np.shape(block)[0], 0])
            ax.set_aspect("auto")
            ax.set_xlabel('observation time (secs)')
            ax.set_ylabel('freq channel')
            ax.set_title('time-freq plot')
            ax.set_aspect('auto')
            plt.colorbar()
            plt.show()

    def read_time_series(self):
        """

        :return:
        """
        if self._fil_file:
            raw = self.file.collapse()  # collapse freq channels
            self.time_series = raw

    def plot_time_series(self):
        """

        :return:
        """
        if self._fil_file:
            if self.time_series is not None:
                raw = self.time_series  # collapse freq channels
            else:
                self.read_time_series()
                raw = self.time_series

            plt.figure()
            ax = plt.subplot(2, 1, 1)
            ax.plot(np.linspace(0, np.size(raw) / self.fs, np.size(raw)), raw)
            ax.set_xlabel("observation time (sec)")
            ax.set_title("raw data (full dataset)")
            ax.set_aspect('auto')
            plt.show()

    def read_database(self, path):
        """
        
        :param path: 
        :return: 
        """
        # TODO download from google docs
        rfiDb = rfiDatabase.RfiDatabase()
        rfiDb.write_dict([path])
        # bands = rfiDb.dictionary.keys()
        int_bands = [rfiDb.dictionary[k]['band'] for k in rfiDb.dictionary.keys()]
        int_dict = dict(zip(int_bands, rfiDb.dictionary.keys()))
        self.database = rfiDb
        return int_dict

    def _count_corrupted(self, close_img):
        """
        
        :param corrupted_samples: fil_rfiObs
        :param close_img: 
        :return: int_dict
        """
        corrupt_block = np.count_nonzero(close_img == 1)
        self.corrupted_samples += corrupt_block

    def apply_threshold(self, block):
        """
        
        :param block: 
        :return: 
        """
        return block < self.threshold

    def clean_mask(self, mask):
        """
        
        :param mask: 
        :return: 
        """
        open_img = rfiUtils.open_blob(mask)
        return rfiUtils.close_blob(open_img)

    def mask_events(self, close_img):
        """
        
        :param close_img: 
        :return: 
        """
        labeled_array, num_features = measurements.label(close_img)
        self._count_corrupted(close_img)
        non_zero_array = labeled_array.nonzero()
        return labeled_array, num_features, non_zero_array

    def find_rfi_events(self, block):
        """
        
        :param block: 
        :return: 
        """
        mask = self.apply_threshold(block)
        close_img = self.clean_mask(mask)
        labeled_array, num_features, non_zero_array = self.mask_events(close_img)
        return labeled_array, num_features, non_zero_array

    def block_events(self, block, int_dict):
        """
        
        :param block: 
        :param int_dict: 
        :return: 
        """
        labeled_array, num_features, non_zero_array = self.find_rfi_events(block)
        feature_range = np.linspace(1, num_features, num=num_features)
        rfi_evs = [rfiEvent.RfiEvent(ev, labeled_array, non_zero_array)
                   for ev in feature_range]
        t_df = self.time[1] - self.time[0]
        [ev.finetune_attr(self.file.header.foff,
                          self.freqs,
                          t_df, self.time) for ev in rfi_evs]
        [ev.find_bands(int_dict) for ev in rfi_evs]
        [ev.find_culprit(self.database.dictionary, int_dict) for ev in rfi_evs]
        return rfi_evs

    def find_obs_event(self, start_time, duration, int_dict):
        """
        
        :param start_time: 
        :param duration: 
        :return: 
        """
        block, num_sam = self.read_time_freq(start_time,
                                             duration)
        return self.block_events(block, int_dict)

    def obs_event(self, vec_length, int_dict):
        """
        
        :param vec_length: 
        :param int_dict: 
        :return: 
        """
        start_vector, duration = self.time_vector(vec_length)
        self.events = [self.find_obs_events(start_vector[sv],
                                            duration,
                                            int_dict) for sv in range(vec_length)]

        #
    # # from rfDB2
    # def get_rfi(self, data, sigma=4):
    #     """
    #     OBSOLETE
    #     Get the rfi for the middle of the window
    #     Data is assumed to be a 2D array (y-axis Time, x-axis Frequency)
    #     and sigma is the standard deviation that is used
    #     :param data: np.array
    #     :param sigma:
    #     :return:
    #     """
    #     mid = data.shape[0] // 2 + 1  # mid point of window
    #     if data.ndim == 2:
    #         med = np.median(data[:, :])
    #     else:
    #         med = np.median(data)
    #     mad = np.median(np.abs(data - med), axis=0)
    #     mad_limit = sigma / 1.4826  # see relation to standard deviation
    #     mid_rfi = (data[mid] > mad_limit * mad + med) + \
    #               (data[mid] < -mad_limit * mad + med)
    #     return mid_rfi, mad_limit
    #
    # # from rfDB2
    # def median_filter(self,
    #                   f_data,
    #                   window=10,
    #                   sigma=4):
    #     """
    #     OBSOLETE
    #     :param f_data:
    #     :param window:
    #     :param sigma:
    #     :return:
    #     """
    #     data = None
    #     if self._rta_file:
    #         data = np.transpose(f_data)
    #     if self._fil_file:
    #         data = f_data
    #     rfi = np.zeros(data.shape)
    #     for t in range(window, data.shape[0] - window):
    #         rfi[t, :] = self.get_rfi(data[t - window:t + window + 1],
    #                                  sigma=sigma)[0]
    #     if self._rta_file:
    #         return np.transpose(rfi)
    #     if self._fil_file:
    #         return rfi
