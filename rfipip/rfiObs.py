#!/usr/bin/env python

# Copyright [2017] [SKA SA]
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from rfipip import rfiDatabase, rfiUtils, rfiEvent, rfiH5
import pandas as pd
from scipy.ndimage import measurements
import ephem


class RfiObservation(object):
    def __init__(self,
                 path,
                 h5_file=False,
                 fil_file=False,
                 rta_file=False):
        """
        
        :param path: 
        :param h5_file: H5 observation file
        :param fil_file: filterbank observation file
        :param rta_file: RTA system observation file
        """
        self.path = path
        self._h5_file = h5_file

        if h5_file:
            self.beamformer = {}

        # MeerKAT coords
        observer = ephem.Observer()
        observer.lon = '21:26:38.0'
        observer.lat = '-30:42:47.4'
        observer.elevation = 1060.0
        observer.epoch = ephem.J2000
        self.observer = observer

        self._fil_file = fil_file
        self._rta_file = rta_file

        self.file = None
        # self.freqs = None
        # self.time = None
        # self.data = None
        # self.bandpass = None

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

    def _open_file(self):
        """
        Open file according to type
        :return:
        """
        if self._h5_file:
            from rfiH5 import RfiH5
            self.file = RfiH5(self.path, self.observer)
        if self._fil_file:
            from rfiFil import RfiFil
            self.file = RfiFil(self.path, self.observer)
        self.file.init()

    def create_bf(self):
        """

        :return:
        """
        for idx, filename in enumerate(self.path):
            h5_obj = rfiH5.RfiH5(filename, self.observer)
            h5_obj.init()

            key = 'pol_%d' % idx
            print('Reading file %s into %s' % (filename, key))
            self.beamformer[key] = {
                                    'h5file': filename,
                                    'metadata': h5_obj.metadata
                                    }

    def verify_h5(self):
        """
        verify parameters that should be the same for observations on both pol files
        basic verification between files of the 2 polarisations
        :return:
        """
        ants = self.beamformer['pol_0']['metadata']['ants']
        if len(ants) > 1:
            raise SystemExit('More than 1 antenna: '
                             'Complete implementation before proceding')
        source_name = self.beamformer['pol_0']['metadata'][ants[0]]['target']
        right_ascension = self.beamformer['pol_0']['metadata'][ants[0]]['ra']
        declination = self.beamformer['pol_0']['metadata'][ants[0]]['dec']
        if self.verbose:
            print('sourceName: %s' % source_name)
            print('rightAscension: %s' % right_ascension)
            print('declination: %s' % declination)

        if self.beamformer['pol_0']['metadata']['fs'] != \
                self.beamformer['pol_1']['metadata']['fs']:
            raise RuntimeError('Files have different sample rates')
        sample_rate = self.beamformer['pol_0']['metadata']['fs']  # sps
        if self.verbose:
            print('samplingClock: %.f Hz' % sample_rate)

        if self.beamformer['pol_0']['metadata']['nchannels'] != \
                self.beamformer['pol_1']['metadata']['nchannels']:
            print('Number of channels differs between the polarizations.')
            raise RuntimeError('channelNumberPol0 %d != channelNumberPol1 %d' %
                               (self.beamformer['pol_0']['metadata']['nchannels'],
                                self.beamformer['pol_1']['metadata']['nchannels']))
        self.beamformer['metadata'] = {'ants': ants,
                                       'source_name': source_name,
                                       'ra': right_ascension,
                                       'dec': declination,
                                       'sample_rate': sample_rate}

    def find_utc(self):
        """
        find UTC sync time for data files
        :return:
        """
        if self.beamformer['pol_0']['metadata']['clk_sync'] != \
                self.beamformer['pol_1']['metadata']['clk_sync']:
            print('System sync timestamp differ between the polarizations.')
            raise RuntimeError()
        sync_ts = beamformer['pol_0']['metadata']['clk_sync']
        if self.verbose:
            print('syncTime: %d' % sync_ts)

        if self.beamformer['pol_0']['metadata']['nsamples'] != \
                self.beamformer['pol_1']['metadata']['nsamples']:
            print('Number of channels differs between the polarizations.')
            raise RuntimeError('channelNumberPol0 %d != channelNumberPol1 %d' %
                               (beamformer['pol_0']['metadata']['nsamples'],
                                beamformer['pol_1']['metadata']['nsamples']))
        nsamples = self.beamformer['pol_0']['metadata']['nsamples']
        if self.verbose:
            print('ADCsnapblock: %d' % nsamples)

        sample_interval = nsamples / self.beamformer['pol_0']['metadata']['fs']
        if self.verbose:
            print('Nyquist sample_interval: %.20f s' % sample_interval)
            print('samplingTime: %.20f ms' % (sample_interval * 1e6))
        self.beamformer['metadata']['syncTime'] = sync_ts
        self.beamformer['metadata']['nsamples'] = nsamples
        self.beamformer['metadata']['sample_interval'] = sample_interval

    def overlap_start(self):
        """
        calculating where both files start overlapping
        :return:
        """
        start_sync_ts = numpy.max([self.beamformer['pol_0']['adc_clks'][0],
                                   self.beamformer['pol_1']['adc_clks'][0]])
        start_sync_ts_pol_a = numpy.where(self.beamformer['pol_0']['adc_clks'] == start_sync_ts)[0][0]
        start_sync_ts_pol_b = numpy.where(self.beamformer['pol_1']['adc_clks'] == start_sync_ts)[0][0]
        if self.verbose:
            print('countADCPol0[0]: %d' % self.beamformer['pol_0']['adc_clks'][0])
            print('countADCPol1[0]: %d' % self.beamformer['pol_1']['adc_clks'][0])
            print('Sync to start timestamp: %d' % start_sync_ts)
            print('\tStart timestamp polA[%d]: %d' % (start_sync_ts_pol_a,
                                                      self.beamformer['pol_0']['adc_clks'][start_sync_ts_pol_a]))
            print('\tStart timestamp polB[%d]: %d' % (start_sync_ts_pol_b,
                                                      self.beamformer['pol_1']['adc_clks'][start_sync_ts_pol_b]))
            print('countADCPol0[%d]: %d' % (start_sync_ts_pol_a,
                                            self.beamformer['pol_0']['adc_clks'][start_sync_ts_pol_a]))
            print('countADCPol1[%d]: %d' % (start_sync_ts_pol_b,
                                            self.beamformer['pol_1']['adc_clks'][start_sync_ts_pol_b]))
            print('startSyncADC: %d' % start_sync_ts)
        return {'start_sync': start_sync_ts,
                'start_sync_pol_a': start_sync_ts_pol_a,
                'start_sync_pol_b': start_sync_ts_pol_b}

    def overlap_end(self):
        """
        calculating where both files end overlaping
        :return:
        """
        end_sync_ts = numpy.min([self.beamformer['pol_0']['adc_clks'][-1],
                                 self.beamformer['pol_1']['adc_clks'][-1]])
        end_sync_ts_pol_a = numpy.where(self.beamformer['pol_0']['adc_clks'] == end_sync_ts)[0][0]
        end_sync_ts_pol_b = numpy.where(self.beamformer['pol_1']['adc_clks'] == end_sync_ts)[0][0]
        if verbose:
            print('countADCPol0[-1]: %d' % self.beamformer['pol_0']['adc_clks'][-1])
            print('countADCPol1[-1]: %d' % self.beamformer['pol_1']['adc_clks'][-1])
            print('Sync to end timestamp: %d' % end_sync_ts)
            print(
                '\tEnd timestamp polA[%d]: %d' % (end_sync_ts_pol_a,
                                                  self.beamformer['pol_0']['adc_clks'][end_sync_ts_pol_a]))
            print(
                '\tEnd timestamp polB[%d]: %d' % (end_sync_ts_pol_b,
                                                  self.beamformer['pol_1']['adc_clks'][end_sync_ts_pol_b]))
            print('countADCPol0[%d]: %d' % (end_sync_ts_pol_a,
                                            self.beamformer['pol_0']['adc_clks'][end_sync_ts_pol_a]))
            print('countADCPol1[%d]: %d' % (end_sync_ts_pol_b,
                                            self.beamformer['pol_1']['adc_clks'][end_sync_ts_pol_b]))
            print('endSyncADC: %d' % end_sync_ts)
        return {'end_sync': end_sync_ts,
                'end_sync_pol_a': end_sync_ts_pol_a,
                'end_sync_pol_b': end_sync_ts_pol_b}

    def sync_h5(self):
        """

        :return:
        """
        self.beamformer['metadata']['syncADC']['start'] = self.overlap_start()
        self.beamformer['metadata']['syncADC']['end'] = self.overlap_end()


    # def read_file(self, start_time=0.0, duration=0.0, polarisation=0):
    #     """
    #
    #     :param start_time: in seconds
    #     :param duration: in seconds
    #     :return:
    #     """
    #     if self.file:
    #         if self._fil_file:
    #             start_sample = start_time / self.file.header.tsamp
    #             if duration == 0:
    #                 nsamples = self.file.header.nsamples - start_sample
    #             else:
    #                 nsamples = duration / self.file.header.tsamp
    #             block = self.file.readBlock(start_sample, nsamples)
    #             self._create_time(start_time, duration, block)
    #             return block, nsamples
    #
    #             # TODO change to channel number
    #             self.data = self.file.read_filterbank(f_start=f_start, f_stop=f_stop)
    #             self.freqs_vector = self.file.freqs_vector
    #         if self._rta_file:
    #             zero_data = self.file['spectra'][:]
    #             # strip data and choose frequency channels
    #             self.data = self._strip_zeros(zero_data)
    #             # [:, ch_start:ch_stop]
    #             # assuming that mode doesnt change in observation
    #             self.mode = self.file['mode'][:][self.time_vector][0]
    #             self._create_freqs()
    #         if self._h5_file:
    #             start_sample = start_time / self.file.sampling_time
    #             if duration == 0:
    #                 nsamples = self.file.header.nsamples - start_sample
    #             else:
    #                 nsamples = duration / self.file.header.tsamp

    # def rfi_block(self, start_time, duration):
    #     """
    #
    #     :param start_time:
    #     :param duration:
    #     :return:
    #     """
    #     # if self._fil_file:
    #     #     block, _ = self.read_time_freq(start_time, duration)
    #     #     return rfiUtils.rfi_threshold(block)
    #
    # def rfi_median(self, vec_length):
    #     """
    #     Run through entire file and give median value for threshold
    #     :param vec_length: int, how many block there is in a file,
    #                         need some logic...
    #     :return:
    #     """
    #     if self._fil_file:
    #         start_vector, duration = self.time_vector(vec_length)
    #         obs_val = [self.rfi_block(st, duration) for _, st in enumerate(start_vector)]
    #         obs_val2array = np.array(obs_val)
    #         self.threshold = np.median(obs_val2array)
    #
    # def sum_dimension(self,
    #                   data,
    #                   axis=0):
    #     """
    #
    #     :param data: 2D array
    #     :param axis:
    #     :return:
    #     """
    #     return data.sum(axis=axis)
    #
    # def read_bandpass(self):
    #     """
    #
    #     :param data:
    #     :return:
    #     """
    #     # if self._rta_file:
    #     #     if self.data is not None:
    #     #         self.bandpass = self.data.sum(axis=0)
    #     # if self._fil_file:
    #     #     self.bandpass = self.file.bandpass()
    #
    # def read_database(self, path):
    #     """
    #
    #     :param path:
    #     :return:
    #     """
    #     # TODO download from google docs
    #     rfi_db = rfiDatabase.RfiDatabase()
    #     rfi_db.write_dict([path])
    #     # bands = rfiDb.dictionary.keys()
    #     int_bands = [rfi_db.dictionary[k]['band'] for k in rfi_db.dictionary.keys()]
    #     int_dict = dict(zip(int_bands, rfi_db.dictionary.keys()))
    #     self.database = rfi_db
    #     return int_dict
    #
    # def _count_corrupted(self, close_img):
    #     """
    #
    #     :param corrupted_samples: fil_rfiObs
    #     :param close_img:
    #     :return: int_dict
    #     """
    #     if self._fil_file:
    #         corrupt_block = np.count_nonzero(close_img == 1)
    #         self.corrupted_samples += corrupt_block
    #
    # def apply_threshold(self, block):
    #     """
    #
    #     :param block:
    #     :return:
    #     """
    #     return block < self.threshold
    #
    # def clean_mask(self, mask):
    #     """
    #
    #     :param mask:
    #     :return:
    #     """
    #     open_img = rfiUtils.open_blob(mask)
    #     return rfiUtils.close_blob(open_img)
    #
    # def mask_events(self, close_img):
    #     """
    #
    #     :param close_img:
    #     :return:
    #     """
    #     if self._fil_file:
    #         labeled_array, num_features = measurements.label(close_img)
    #         self._count_corrupted(close_img)
    #         non_zero_array = labeled_array.nonzero()
    #         return labeled_array, num_features, non_zero_array
    #
    # def percentage_rfi(self, vec_length):
    #     """
    #
    #     :param vec_length:
    #     :return:
    #     """
    #     if self._fil_file:
    #         _, duration = self.time_vector(vec_length)
    #         num_sam = long(duration / self.file.header.tsamp)
    #         file_sam = vec_length * num_sam * self.file.header.nchans
    #         return rfiUtils.percentage(self.corrupted_samples, file_sam)
    #
    # def find_rfi_events(self, block):
    #     """
    #
    #     :param block:
    #     :return:
    #     """
    #     if self._fil_file:
    #         mask = self.apply_threshold(block)
    #         close_img = self.clean_mask(mask)
    #         labeled_array, num_features, non_zero_array = self.mask_events(close_img)
    #         return labeled_array, num_features, non_zero_array
    #
    # def block_events(self, block, int_dict):
    #     """
    #
    #     :param block:
    #     :param int_dict:
    #     :return:
    #     """
    #     if self._fil_file:
    #         labeled_array, num_features, non_zero_array = self.find_rfi_events(block)
    #         feature_range = np.linspace(1, num_features, num=num_features)
    #         rfi_evs = [rfiEvent.RfiEvent(ev, labeled_array, non_zero_array)
    #                    for ev in feature_range]
    #         t_df = self.time[1] - self.time[0]
    #         [ev.finetune_attr(self.file.header.foff,
    #                           self.freqs,
    #                           t_df, self.time) for ev in rfi_evs]
    #         [ev.find_bands(int_dict) for ev in rfi_evs]
    #         [ev.find_culprit(self.database.dictionary, int_dict) for ev in rfi_evs]
    #         return rfi_evs
    #
    # def find_obs_event(self, start_time, duration, int_dict):
    #     """
    #
    #     :param start_time:
    #     :param duration:
    #     :return:
    #     """
    #     if self._fil_file:
    #         block, num_sam = self.read_time_freq(start_time,
    #                                              duration)
    #         return self.block_events(block, int_dict)
    #
    # def obs_events(self, vec_length, int_dict):
    #     """
    #
    #     :param vec_length:
    #     :param int_dict:
    #     :return:
    #     """
    #     if self._fil_file:
    #         start_vector, duration = self.time_vector(vec_length)
    #         self.events = [self.find_obs_event(start_vector[sv],
    #                                            duration,
    #                                            int_dict) for sv in range(vec_length)]
    #
    # def write2csv(self,
    #               h5_name='rfi_measurements.h5',
    #               return_h5=False):
    #     """
    #
    #     :param h5_name:
    #     :param return_h5:
    #     :return:
    #     """
    #     csv_name = self.file.header.source_name + '_' + \
    #                str(self.file.header.src_dej) + '_' + \
    #                str(self.file.header.src_raj) + '.csv'
    #     # TODO get rid of for loops
    #     # TODO reset idx in csv file
    #     columns = ('event',
    #                'c_freq',
    #                'bw',
    #                't_start',
    #                'duration',
    #                'culprit',
    #                'description',
    #                'band')
    #     # check if h5 file exist
    #     # TODO if exist, stop, delete or change name
    #     with pd.HDFStore(h5_name) as store:
    #         for i, x in enumerate(self.events):
    #             training_set = pd.DataFrame(columns=columns)
    #             pd_idx = 0
    #             for eb in x:
    #                 if type(eb.culprit) is list:
    #                     for b in range(len(eb.culprit)):
    #                         training_set.loc[pd_idx] = [eb.event,
    #                                                     eb.c_freq,
    #                                                     eb.bw,
    #                                                     eb.t_start,
    #                                                     eb.duration,
    #                                                     eb.culprit[b],
    #                                                     eb.description[b],
    #                                                     eb.band]
    #                         pd_idx += 1
    #
    #                 else:
    #                     training_set.loc[pd_idx] = [eb.event,
    #                                                 eb.c_freq,
    #                                                 eb.bw,
    #                                                 eb.t_start,
    #                                                 eb.duration,
    #                                                 eb.culprit,
    #                                                 eb.description,
    #                                                 eb.band]
    #                     pd_idx += 1
    #             store.put('test',
    #                       training_set,
    #                       format='table',
    #                       data_columns=True,
    #                       append=True,
    #                       min_itemsize={'band': 150,
    #                                     'description': 500})
    #         store['test'].to_csv(csv_name)
    #         if return_h5:
    #             return store['test']
    #

    # def _strip_zeros(self, data):
    #     """
    #     Strip zero columns from data
    #     :param data:
    #     :return:
    #     """
    #     if self._rta_file:
    #         # ratty samples every 2 sec so get rid of the zeros
    #         num_sec = data.shape[0]  # 3600 seconds
    #         # [0...3599]
    #         full_sec_idx = range(num_sec)
    #         # find idx with zeros
    #         zero_idx = list(np.where(~data.any(axis=1))[0])
    #         # delete zero idx from full idx list
    #         data_idx = np.delete(full_sec_idx, zero_idx, axis=0)
    #         stripped_data = data[data_idx]
    #         self.time_vector = data_idx
    #         return stripped_data
    #     if self._fil_file:
    #         num_sec = data.shape[1]  # time_vector samples
    #         full_sec_idx = range(num_sec)
    #         # find idx with zeros
    #         zero_idx = list(np.where(~data.any(axis=0))[0])
    #         # delete zero idx from full idx list
    #         data_idx = np.delete(full_sec_idx, zero_idx, axis=0)
    #         stripped_data = data[:, data_idx]
    #         self.time_vector = data_idx
    #         return stripped_data

    # def _open_rta_file(self):
    #     """
    #     Open RTA h5 file
    #     :return:
    #     """
    #     if self._rta_file:
    #         import h5py
    #         self.file = h5py.File(self.path, mode='r+')

    # def plot_bandpass(self,
    #                   data=None,
    #                   ch_start=0,
    #                   ch_stop=0):
    #     """
    #     Works for summing 2D data along Y axis
    #     :param data:
    #     :return:
    #     """
    #     if self._rta_file:
    #         fig, ax = plt.subplots()
    #         ax.plot(self.freqs_vector[ch_start:ch_stop]/1e9, data)
    #         ax.set_xlabel('Frequency')
    #         ax.set_ylabel('Power')
    #         ax.set_aspect('auto')
    #         plt.show()
    #     if self._fil_file:
    #         if self.bandpass is not None:
    #             bp = self.bandpass  # collapse freq channels
    #         else:
    #             self.read_bandpass()
    #             bp = self.bandpass
    #         plt.figure()
    #         ax = plt.subplot(1, 1, 1)
    #         ax.plot(bp)
    #         ax.set_xlabel('Freq channel')
    #         ax.set_title('Bandpass')
    #         ax.set_aspect('auto')
    #         plt.show()
    #
    # # TODO axis ranges incorrect - fix
    # def plot_spectrum(self,
    #                   data,
    #                   ch_start,
    #                   ch_stop):
    #     """
    #
    #     :param data:
    #     :param ch_start:
    #     :param ch_stop:
    #     :return:
    #     """
    #     fig, ax = plt.subplots()
    #     if self._rta_file:
    #         ax.imshow(data,
    #                   extent=[self.freqs_vector[ch_start] / 1e9,
    #                           self.freqs_vector[ch_stop] / 1e9,
    #                           self.time_vector[-1],
    #                           self.time_vector[0]],
    #                   cmap=cm.Blues)
    #         ax.set_xlabel('Frequency [MHz]')
    #         ax.set_ylabel('Time (seconds)')
    #     if self._fil_file:
    #         ax.imshow(data,
    #                   extent=[self.time_vector[0],
    #                           self.time_vector[-1],
    #                           self.freqs_vector[ch_start],
    #                           self.freqs_vector[ch_stop]],
    #                   cmap=cm.Blues)
    #         ax.set_ylabel('Frequency [MHz]')
    #         ax.set_xlabel('Time (seconds)')
    #
    #     ax.set_aspect('auto')
    #     plt.show()
    #
    # def plot_mask(self,
    #               mask,
    #               start_time):
    #     start_sample = long(start_time / self.file.header.tsamp)
    #
    #     plt.figure()
    #     ax = plt.subplot(1, 1, 1)
    #     plt.imshow(mask, extent=[start_sample * self.file.header.tsamp,
    #                              (start_sample + np.shape(mask)[1]) * self.file.header.tsamp,
    #                              np.shape(mask)[0], 0])
    #     ax.set_aspect("auto")
    #     ax.set_xlabel('observation time_vector (secs)')
    #     ax.set_ylabel('freq channel')
    #     ax.set_title('mask time_vector-freq plot')
    #     ax.set_aspect('auto')
    #     plt.colorbar()
    #     plt.show()
    #
    # def plot_time_freq(self,
    #                    start_time,
    #                    duration,
    #                    zeros=False):
    #     """
    #
    #     :param start_time:
    #     :param duration:
    #     :return:
    #     """
    #     if self._fil_file:
    #         start_sample = long(start_time / self.file.header.tsamp)
    #         block, nsamples = self.read_time_freq(start_time, duration)
    #         if zeros:
    #             new_block = self._strip_zeros(block)
    #             block = new_block
    #
    #         plt.figure()
    #         ax = plt.subplot(1, 1, 1)
    #         plt.imshow(block, extent=[start_sample * self.file.header.tsamp,
    #                                   (start_sample + np.shape(block)[1]) * self.file.header.tsamp,
    #                                   np.shape(block)[0], 0])
    #         ax.set_aspect("auto")
    #         ax.set_xlabel('observation time_vector (secs)')
    #         ax.set_ylabel('freq channel')
    #         ax.set_title('time_vector-freq plot')
    #         ax.set_aspect('auto')
    #         plt.colorbar()
    #         plt.show()
    #
    #
    # def plot_time_series(self):
    #     """
    #
    #     :return:
    #     """
    #     if self._fil_file:
    #         if self.time_series is not None:
    #             raw = self.time_series  # collapse freq channels
    #         else:
    #             self.read_time_series()
    #             raw = self.time_series
    #
    #         plt.figure()
    #         ax = plt.subplot(2, 1, 1)
    #         ax.plot(np.linspace(0, np.size(raw) / self.fs, np.size(raw)), raw)
    #         ax.set_xlabel("observation time_vector (sec)")
    #         ax.set_title("raw data (full dataset)")
    #         ax.set_aspect('auto')
    #         plt.show()