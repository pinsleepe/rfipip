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

import numpy as np


class RfiEvent(object):
    """

    """
    def __init__(self,
                 event,
                 labeled_array,
                 non_zero_arr):

        self.event = event
        self.c_freq = None
        self.bw = None
        self.t_start = None
        self.duration = None
        self.channel = None
        self.culprit = []
        self.description = []
        self.band = []

        self.init_info = None

        self.init(labeled_array, non_zero_arr)

    def init(self, labeled_array, non_zero_arr):
        """
        
        :param labeled_array: 
        :return: 
        """
        # non_zero_arr = labeled_array.nonzero()
        f, t = self.look_4_event(labeled_array, non_zero_arr)
        self.assign_attr(f, t)

    def look_4_event(self, labeled_array, non_zero_array):
        """
        
        :return: 
        """
        ev = np.where(labeled_array[non_zero_array] == self.event)
        x = non_zero_array[0][ev]
        y = non_zero_array[1][ev]
        return x, y

    def assign_attr(self, x, y):
        """
        
        :param x: 
        :param y: 
        :return: 
        """
        # c_freq 0
        # channel 1
        # t_start 2
        # duration 3
        self.init_info = [x[0],
                          x.max() - x[0],
                          y[0],
                          y.max() - y[0]]
        self.channel = self.init_info[1]

    def finetune_attr(self, foff,
                      freqs_v,
                      t_df,
                      time_v):
        """
        
        :param foff: header.foff
        :param freqs_vector: fil_rfiObs.freqs_vector
        :param t_df: fil_rfiObs.time_vector[1] - fil_rfiObs.time_vector[0]
        :param time_v: fil_rfiObs.time_vector
        :return: 
        """
        # check if the event occupies more than one channel
        # yes
        if self.init_info[1] > 0:
            # freq channels times BW
            temp_bw = self.init_info[1] * foff
            # freq of middle channel
            temp_freq = freqs_v[self.init_info[0]] + temp_bw / 2.0
        # no
        else:
            temp_freq = freqs_v[self.init_info[0]]
            temp_bw = foff
        # duration
        if self.init_info[3] == 0:
            temp_dur = t_df
        else:
            temp_dur = self.init_info[3] * t_df
        self.c_freq = temp_freq
        self.duration = temp_dur
        self.bw = temp_bw
        self.t_start = time_v[self.init_info[2]]

    def find_bands(self, range_dict):
        """
        
        :param range_dict: dict {0: '47.0-87.5'}
        :return: 
        """

        self.band = [key
                     for key in range_dict.keys()
                     if float(range_dict[key].split('-')[0]) <=
                     self.c_freq <=
                     float(range_dict[key].split('-')[1])]

    def _find_culprit_freq(self,
                           full_dict,
                           range_dict,
                           band,
                           f_range=None):
        """
        
        :param full_dict: 
        :param range_dict: 
        :param band: 
        :return: 
        """
        if f_range:
            f_range = f_range
        else:
            f_range = 1
        if full_dict[range_dict[band]]['frequencies'] is not []:
            culprit = [(1, full_dict[range_dict[band]]['frequencies'][av_f])
                       for av_f in full_dict[range_dict[band]]['frequencies'].keys()
                       if (float(av_f) - f_range) <= self.c_freq <= (float(av_f) + f_range)]
            if culprit == []:
                culprit = [(0, 'Unknown')]
        else:
            culprit = [(0, 'Unknown')]
        return culprit[0]

    def find_culprit(self, full_dict, range_dict):
        """
        
        :param full_dict: 
        :param range_dict: 
        :param f_range: 
        :return: 
        """
        culprit_info = [self._find_culprit_freq(full_dict, range_dict, band)
                        for band in self.band]
        self.culprit = [cul[0] for cul in culprit_info]
        self.description = [cul[1] for cul in culprit_info]
        self.clean_culprit_info()

    def clean_culprit_info(self):
        """
        
        :return: 
        """
        # clean up descriptions
        known_idx = [i for i, x in enumerate(self.culprit) if x == 1]
        # more than one culprit identified
        if len(known_idx) > 1:
            # check if all are the same
            list_ones = [self.description[x] for i, x in enumerate(known_idx)]
            # yes
            if list_ones.count(list_ones[0]) == len(list_ones):
                self.culprit = 1
                self.description = list_ones[0]
            # no
            else:
                # add new label?
                unique = set(list_ones)
                self.description = list(unique)
                self.culprit = [1]*len(unique)
        # one culprit identified
        if len(known_idx) == 1:
            self.culprit = 1
            clean_d = self.description[known_idx[0]]
            self.description = clean_d
        # no culprit identified
        if len(known_idx) == 0:
            self.culprit = 0
            self.description = 'Unknown'

        # integers to a string for data frame
        clean_band = '-'.join(str(x) for x in self.band)
        self.band = clean_band

        return
