# Author: Maciej Serylak
# https://github.com/mserylak/scripts/blob/master/fastH5.py
# Modified by Bruce Merry
# Modified by Ruby van Rooyen
# Modified by Monika Obrocka (Aug 2017)

import h5py
import numpy as np
from datetime import datetime


class RfiH5(object):
    def __init__(self, path,
                 observer):
        """
        Initialise object 
        :param path: 
        :param observer: ephem.Observer() object
        """

        self.path = path
        self.file = None
        self.observer = observer
        self.verbose = False
        self.metadata = {}
        self.clockcounts = None

    def init(self):
        """
        Init BF observation
        :return: 
        """
        self.open_file(self.path)
        self.fill_metadata()
        # do stuff
        self.file.close()

    def open_file(self, filename):
        """
        Loading the files.
        :return:
        """
        self.file = h5py.File(filename, 'r')

    def fill_metadata(self):
        """
        Fill the metadata with info from H5
        :return: 
        """
        self.read_timestamps()
        self.adc_counts_diff()
        self.read_data()
        self.freq_metadata()
        self.fill_observer()
        self.ant_metadata()

    def read_timestamps(self):
        """
        Read timestamp type
        :return: 
        """
        self.metadata['ts_type'] = self.file['Data/timestamps'].attrs['timestamp_type']
        if self.verbose:
            print('Get %s timestamps from file' % self.metadata['ts_type'])
        self.clockcounts = self.file['Data/timestamps'][:]
        if self.verbose:
            print('Read %d timestamps' % len(self.clockcounts))

    def adc_counts_diff(self):
        """
        difference between adc_counts to check for missing packets
        :return: 
        """
        nts_adc_snapblock = np.diff(self.clockcounts)[0]  # all should be 8192
        self.metadata['nsamples'] = nts_adc_snapblock
        if self.verbose:
            print('Nr time samples per ADC snap block = %d' % nts_adc_snapblock)
        if (np.where(np.diff(self.clockcounts) != nts_adc_snapblock)[0].size) > 0:
            print('Missing spectra in: %s' % self.path)
            print(np.where(np.diff(self.clockcounts) != nts_adc_snapblock)[0])
            raise RuntimeError('Missing spectra')

    def read_data(self):
        """
        read data, assume layout (n_chans, n_timestamps, n_complex)
        :return: 
        """
        self.metadata['rawdata'] = 'Data/bf_raw'
        rawdata = self.file['Data/bf_raw']
        if self.verbose:
            print('Reading raw data from file')
        [n_chans, n_ts, n_c] = rawdata.shape
        self.metadata['nchannels'] = n_chans
        self.metadata['nspectra'] = n_ts
        if abs(n_ts - len(self.clockcounts)) > 0:
            raise RuntimeError('Not expected data format: '
                               'Nr timestamps do not match')  # noqa
        if abs(n_chans - int(self.file['TelescopeModel/cbf'].attrs['n_chans'])) > 0:
            raise RuntimeError('Not expected data format: '
                               'Nr channels do not match')  # noqa

    def freq_metadata(self):
        """
        Frequency specific metadata
        :return: 
        """
        scale_factor_timestamp = float(self.file['TelescopeModel/cbf'].attrs['scale_factor_timestamp'])  # noqa
        tot_ts = self.clockcounts[-1] - self.clockcounts[0]
        n_sec = float(tot_ts) / scale_factor_timestamp
        sps = int(tot_ts / n_sec)
        self.metadata['clks_per_sec'] = self.metadata['nspectra'] / n_sec
        self.metadata['fs'] = sps
        self.metadata['clk_sync'] = self.file['/TelescopeModel/cbf'].attrs['sync_time']
        self.metadata['cenfreq'] = self.file['TelescopeModel/cbf'].attrs['center_freq']
        self.metadata['bandwidth'] = self.file['TelescopeModel/cbf'].attrs['bandwidth']

    def fill_observer(self):
        """
        
        :return: 
        """
        self.observer.temp = self.file['TelescopeModel/anc/air_temperature'][0]['value']
        self.observer.pressure = self.file['TelescopeModel/anc/air_pressure'][0]['value']
        self.observer.horizon = np.radians(15)

    def ant_metadata(self):
        """
        antenna specific metadata
        :return: 
        """
        ants = [str(key) for key in self.file['TelescopeModel'].keys() if 'm0' in key]
        self.metadata['ants'] = ants
        for ant in ants:
            self.metadata[ant] = {}
            self.metadata[ant]['target'] = self.file['TelescopeModel'][ant]['target'][0]['value']
            self.metadata[ant]['az'] = self.file['TelescopeModel'][ant]['pos_actual_scan_azim'][0]['value']
            self.metadata[ant]['el'] = self.file['TelescopeModel'][ant]['pos_actual_scan_elev'][0]['value']
            self.observer.date = datetime.utcfromtimestamp(self.file['TelescopeModel'][ant]['pos_actual_scan_azim'][0]['timestamp'])
            ra, dec = self.observer.radec_of(np.radians(self.metadata[ant]['az']), np.radians(self.metadata[ant]['el']))
            self.metadata[ant]['ra'] = str(ra)
            self.metadata[ant]['dec'] = str(dec)

    #     """
    #     Getting number of channels from each file.
    #     :return:
    #     """
    #     self.p0_nchans = self.p0_data["Data/bf_raw"].shape[0]
    #     self.p1_nchans = self.p1_data["Data/bf_raw"].shape[0]
    #     self.nchans_equal = self._nchans_equal()
    #     if self.nchans_equal:
    #         self.nchans = self.p0_nchans
    #
    # def _nchans_equal(self):
    #     """
    #     Checking if number of channels is the same in each file.
    #     :return:
    #     """
    #     equal = True
    #     if self.p0_nchans != self.p1_nchans:
    #         equal = False
    #     return equal
    #
    # def spectra(self):
    #     """
    #     Getting number of spectra from each file.
    #     :return:
    #     """
    #     self.p0_spectra = self.p0_data["Data/bf_raw"].shape[1]
    #     self.p1_spectra = self.p1_data["Data/bf_raw"].shape[1]
    #
    # def adc_count(self):
    #     """
    #     Getting ADC counts from each file.
    #     :return:
    #     """
    #     self.p0_adc_count = self.p0_data["Data/timestamps"][:]
    #     self.p1_adc_count = self.p1_data["Data/timestamps"][:]
    #
    # def find_subset(self):
    #     """
    #     Calculating where both files are overlaping.
    #     :return:
    #     """
    #     p0_shape = self.p0_adc_count.shape
    #     p1_shape = self.p1_adc_count.shape
    #
    #     if p0_shape > p1_shape:
    #         intersect = np.intersect1d(self.p0_adc_count,
    #                                    self.p1_adc_count,
    #                                    assume_unique=True)
    #     elif p0_shape < p1_shape:
    #         intersect = np.intersect1d(self.p1_adc_count,
    #                                    self.p0_adc_count,
    #                                    assume_unique=True)
    #
    #     else:
    #         intersect = np.intersect1d(self.p1_adc_count,
    #                                    self.p0_adc_count,
    #                                    assume_unique=True)
    #
    #     adc0_start = np.where(self.p0_adc_count == intersect[0])[0][0]
    #     adc0_stop = np.where(self.p0_adc_count == intersect[-1])[0][0]
    #     adc1_start = np.where(self.p1_adc_count == intersect[0])[0][0]
    #     adc1_stop = np.where(self.p1_adc_count == intersect[-1])[0][0]
    #
    #     self.p0_adc_idx = np.linspace(adc0_start,
    #                                   adc0_stop,
    #                                   num=intersect.shape[0],
    #                                   dtype=np.int)
    #     self.p1_adc_idx = np.linspace(adc1_start,
    #                                   adc1_stop,
    #                                   num=intersect.shape[0],
    #                                   dtype=np.int)
    #
    # def missing_packets(self):
    #     """
    #     Getting difference between countADC to check for missing packets.
    #     8192 ??
    #     :return:
    #     """
    #     differencePol0 = np.diff(self.p0_adc_count[self.p0_adc_idx[0]:])
    #     differencePol1 = np.diff(self.p1_adc_count[self.p1_adc_idx[0]:])
    #     breaksPol0 = np.where(differencePol0 != self.obs_bw/1e6)[0]
    #     breaksPol1 = np.where(differencePol1 != self.obs_bw/1e6)[0]
    #     self.p0_missing_spectra = breaksPol0
    #     self.p1_missing_spectra = breaksPol1
    #
    # def find_sync_time(self):
    #     """
    #     Find sync time in the data files.
    #     :return: unix time in seconds
    #     """
    #     try:
    #         syncTime = self.p0_data["/TelescopeModel/cbf"].attrs['sync_time']
    #         self.sync_time = syncTime
    #         syncFactor = self.p0_data["/TelescopeModel/cbf"].attrs['scale_factor_timestamp']
    #         self.sync_time_factor = syncFactor
    #     except KeyError:
    #         print "Data does not have sync time in the header! " \
    #               "Specify it manually in the script!"
    #
    # def find_obs_start_time(self):
    #     """
    #     Calculate observation start time,
    #     :return: observation start time in seconds
    #     """
    #     self.obs_start_time = self.p0_adc_count[self.p0_adc_idx[0]] / self.sync_time_factor
    #
    # def find_unix_time(self):
    #     """
    #     Convert observation start time to unix time
    #     :return:
    #     """
    #     self.unix_time = float(self.sync_time) + self.obs_start_time
    #
    # def unix2mjd(self):
    #     """
    #     Convert unix time to MJD
    #     :return:
    #     """
    #     startTimeMJD = katpoint.Timestamp(self.unix_time)
    #     self.mjd = startTimeMJD.to_mjd()
    #
    # def find_freqs(self):
    #     """
    #     Find centre, top and bottom frequencies from the file.
    #     :return:
    #     """
    #     freqCent = float(self.p0_data["/TelescopeState/obs_params"][:][6][1].split()[1].strip("'"))
    #     self.centre_freq = freqCent
    #
    #     freqTop = freqCent + (((self.nchans / 2) - 1) * self.channel_bw)
    #     self.top_freq = freqTop
    #
    #     freqBottom = freqCent - ((self.nchans / 2) * self.channel_bw)
    #     self.bottom_freq = freqBottom
    #
    # def find_radec(self):
    #     """
    #     Find Ra, Dec coordinates for az, el pointing
    #     :return:
    #     """
    #     az = self.p0_data["/TelescopeState/obs_params"][:][15][1].split()[3][:-3]
    #     el = self.p0_data["/TelescopeState/obs_params"][:][15][1].split()[2][:-1]
    #
    #     self.azimuth = ephem.degrees(az)  # radians
    #     self.elevation = ephem.degrees(el)  # radians
    #
    #     # katpoint.construct_azel_target(az, el).radec(timestamp, antenna)
    #
    # def check_block_goodness(self, voltagefile, min_percent=10):
    #     """
    #     Check how many zero channels is in the block from H5 file
    #     :param voltagefile: numpy array
    #     :param min_percent: min percent ob bad channels that are tolerated
    #     :return:
    #     """
    #     all_zeros = [not np.any(voltagefile[ch, :]) for ch in range(self.nchans)]
    #     good_channels = all_zeros.count(False)
    #     min_goodness = int(self.nchans - (self.nchans / 100.0) * min_percent)
    #     skip = False
    #     if good_channels < min_goodness:
    #         skip = True
    #     return skip
