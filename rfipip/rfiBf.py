# Author: Maciej Serylak
# https://github.com/mserylak/scripts/blob/master/fastH5.py
# Modified by Bruce Merry
# Modified by Ruby van Rooyen
# Modified by Monika Obrocka (Aug 2017)

import numpy as np
import ephem
from datetime import datetime


class RfiBf(object):
    def __init__(self,
                 pol0,
                 pol1,
                 verbose=False,
                 debug=False):
        """
         Class to handle BF eng mode observations
        :param pol0: RfiH5 object
        :param pol1: RfiH5 object
        :param verbose: Show various parameters and results
        :param debug: Show very verbose output for debugging
        """
        self.pol0 = pol0
        self.pol1 = pol1
        self.beamformer = None
        self.verbose = verbose
        self.debug = debug

        self._init()

    def _init(self):
        """
        
        :return: 
        """
        beamformer = {}
        for idx, filename in enumerate([self.pol0, self.pol1]):
            key = 'pol_%d' % idx
            print('Reading file %s into %s' % (filename, key))
            [metadata, adc_clks] = filename.metadata, filename.clockcounts
            beamformer[key] = {
                'h5file': filename,
                'metadata': metadata,
                'adc_clks': adc_clks,
            }
        self.beamformer = beamformer

        self.verify_parameters()
        self.obs_metadata()

    def verify_parameters(self):
        """
        basic verification between files of the 2 polarisations
        :return: 
        """
        ants = self.check_ants()
        nchannels = self.check_nchannels()
        sample_rate = self.check_sample_rate()
        sync_ts = self.sync_time()
        nsamples = self.check_nsamples()
        info = {'ants': ants,
                'sample_rate': sample_rate,
                'nchannels': nchannels,
                'sync_time': sync_ts,
                'nsamples': nsamples}
        self.beamformer['info'] = info
        start_sync = self.files_overlap_start()
        end_sync = self.files_overlap_stop()
        self.beamformer['info']['sync'] = {'start': start_sync,
                                           'stop': end_sync}

    def check_sample_rate(self):
        """
        
        :return: 
        """
        if self.beamformer['pol_0']['metadata']['fs'] \
                != self.beamformer['pol_1']['metadata']['fs']:
            raise RuntimeError('Files have different sample rates')
        sample_rate = self.beamformer['pol_0']['metadata']['fs']  # sps
        if self.verbose:
            print('samplingClock: %.f Hz' % sample_rate)
        return sample_rate

    def check_nchannels(self):
        """
        
        :return: 
        """
        if self.beamformer['pol_0']['metadata']['nchannels'] \
                != self.beamformer['pol_1']['metadata']['nchannels']:
            print('Number of channels differs between the polarizations.')
            raise RuntimeError('channelNumberPol0 %d != channelNumberPol1 %d' %
                               (self.beamformer['pol_0']['metadata']['nchannels'],
                                self.beamformer['pol_1']['metadata']['nchannels']))
        nchannels = self.beamformer['pol_0']['metadata']['nchannels']
        if self.verbose:
            print('complexChannels: %d' % nchannels)
        return nchannels

    def check_ants(self):
        """
        
        :return: 
        """
        ants = self.beamformer['pol_0']['metadata']['ants']
        if len(ants) > 1:
            raise SystemExit('More than 1 antenna: '
                             'Complete implementation before proceding')
        return ants

    def sync_time(self):
        """
        find UTC sync time_vector for data files
        :return: 
        """
        if self.beamformer['pol_0']['metadata']['clk_sync'] \
                != self.beamformer['pol_1']['metadata']['clk_sync']:
            print('System sync timestamp differ between the polarizations.')
            raise RuntimeError()
        sync_ts = self.beamformer['pol_0']['metadata']['clk_sync']
        if self.verbose:
            print('syncTime: %d' % sync_ts)
        return sync_ts

    def check_nsamples(self):
        """
        
        :return: 
        """
        if self.beamformer['pol_0']['metadata']['nsamples'] \
                != self.beamformer['pol_1']['metadata']['nsamples']:
            print('Number of channels differs between the polarizations.')
            raise RuntimeError('channelNumberPol0 %d != channelNumberPol1 %d' %
                               (self.beamformer['pol_0']['metadata']['nsamples'],
                                self.beamformer['pol_1']['metadata']['nsamples']))
        nsamples = self.beamformer['pol_0']['metadata']['nsamples']
        if self.verbose:
            print('ADCsnapblock: %d' % nsamples)
        return nsamples

    def files_overlap_start(self):
        """
        calculating where both files start overlapping
        :return: 
        """
        start_sync_ts = np.max([self.beamformer['pol_0']['adc_clks'][0], self.beamformer['pol_1']['adc_clks'][0]])
        start_sync_ts_pol_a = np.where(self.beamformer['pol_0']['adc_clks'] == start_sync_ts)[0][0]
        start_sync_ts_pol_b = np.where(self.beamformer['pol_1']['adc_clks'] == start_sync_ts)[0][0]
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
        sync = {'sync_ts': start_sync_ts,
                'sync_polA': start_sync_ts_pol_a,
                'sync_polB': start_sync_ts_pol_b}
        return sync

    def files_overlap_stop(self):
        """
        calculating where both files end overlaping
        :return: 
        """
        end_sync_ts = np.min([self.beamformer['pol_0']['adc_clks'][-1],
                              self.beamformer['pol_1']['adc_clks'][-1]])
        end_sync_ts_pol_a = np.where(self.beamformer['pol_0']['adc_clks'] ==
                                    end_sync_ts)[0][0]
        end_sync_ts_pol_b = np.where(self.beamformer['pol_1']['adc_clks'] ==
                                     end_sync_ts)[0][0]
        if self.verbose:
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
        sync = {'sync_ts': end_sync_ts,
                'sync_polA': end_sync_ts_pol_a,
                'sync_polB': end_sync_ts_pol_b}
        return sync

    def obs_metadata(self):
        """
        
        :return: 
        """
        source_name = self.beamformer['pol_0']['metadata']['ants'][0]['target']
        right_ascension = self.beamformer['pol_0']['metadata']['ants'][0]['ra']
        declination = self.beamformer['pol_0']['metadata']['ants'][0]['dec']
        self.beamformer['info']['target'] = {'source_name': source_name,
                                             'right_ascension': right_ascension,
                                             'declination': declination}
        if self.verbose:
            print('sourceName: %s' % source_name)
            print('rightAscension: %s' % right_ascension)
            print('declination: %s' % declination)

        start_ts = self.beamformer['info']['sync']['start']['sync_ts'] / self.beamformer['info']['sample_rate']
        unix_start_ts = float(self.beamformer['info']['sync_time']) + start_ts
        unix_dt = datetime.utcfromtimestamp(unix_start_ts).strftime('%Y-%m-%d %H:%M:%S.%f')
        start_mjd = ephem.julian_date(ephem.Date(unix_dt)) - 2400000.5  # convert to MJD
        self.beamformer['info']['time_info'] = {'obsStartTime': start_ts,
                                                'unix_start_ts': unix_start_ts,
                                                'obsSyncDate': unix_dt,
                                                'startTimeMJD': start_mjd}
        cen_freq = self.beamformer['pol_0']['metadata']['cenfreq']
        bandwidth = self.beamformer['pol_0']['metadata']['bandwidth']
        channel_bw = bandwidth / self.beamformer['info']['nchannels']
        if self.verbose:
            print('obsStartTime: %.12f' % start_ts)
            print('obsSyncDate: %s' % unix_dt)
            print('startTimeMJD: %.12f' % start_mjd)
            print('freqCent: %f' % cen_freq)
            print('channelBW: %.10f MHz' % (channel_bw * 1e-6))
        upper_freq = cen_freq + bandwidth / 2.
        lower_freq = cen_freq - bandwidth / 2.
        self.beamformer['info']['freq_info'] = {'cen_freq': cen_freq,
                                                'bandwidth': bandwidth,
                                                'channel_bw': channel_bw,
                                                'freqTop': upper_freq,
                                                'freqBottom': lower_freq}
        # Getting number of spectra from each file.
        nspectra_pol_a = self.beamformer['pol_0']['metadata']['nspectra']
        nspectra_pol_b = self.beamformer['pol_1']['metadata']['nspectra']
        if self.verbose:
            print('freqTop: %.10f MHz' % (upper_freq * 1e-6))
            print('freqBottom: %.10f MHz' % (lower_freq * 1e-6))
            print('spectraNumberPol0: %d' % nspectra_pol_a)
            print('spectraNumberPol1: %d' % nspectra_pol_b)
        self.beamformer['info']['spectra'] = {'spectraNumberPol0': nspectra_pol_a,
                                              'spectraNumberPol1': nspectra_pol_b}

