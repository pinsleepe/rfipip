import numpy as np

class RfiEvent(object):
    """

    """
    def __init__(self,
                 event,
                 labeled_array):

        self.event = event
        self.c_freq = None
        self.bw = None
        self.t_start = None
        self.duration = None
        self.channel = None
        self.culprit = []
        self.description = []
        self.band = []

        self.init(labeled_array)

    def init(self, labeled_array):
        """
        
        :param labeled_array: 
        :return: 
        """
        f, t = self.look_4_event(labeled_array)
        self.assign_attr(f, t)

    def look_4_event(self, labeled_array):
        """
        
        :return: 
        """
        x, y = np.where(labeled_array == self.event)
        return x, y

    def assign_attr(self, x, y):
        """
        
        :param x: 
        :param y: 
        :return: 
        """
        self.c_freq = x[0]
        self.channel = x.max() - x[0]
        self.t_start = y[0]
        self.duration = y.max() - y[0]

    def finetune_attr(self, foff,
                      freqs_v,
                      t_df,
                      time_v):
        """
        
        :param foff: header.foff
        :param freqs: fil_rfiObs.freqs
        :param t_df: fil_rfiObs.time[1] - fil_rfiObs.time[0]
        :param time_v: fil_rfiObs.time
        :return: 
        """
        if self.channel > 0:
            # freq channels times BW
            temp_bw = self.channel * foff
            # freq of middle channel
            temp_freq = freqs_v[self.c_freq] + temp_bw / 2.0
            # duration
            temp_dur = self.duration * t_df
        # no
        else:
            temp_freq = freqs_v[self.c_freq]
            temp_bw = foff
            temp_dur = t_df

        self.c_freq = temp_freq
        self.duration = temp_dur
        self.bw = temp_bw
        self.t_start = time_v

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
        return culprit

    def find_culprit(self, full_dict, range_dict):
        """
        
        :param full_dict: 
        :param range_dict: 
        :param f_range: 
        :return: 
        """
        culprit_info = [self._find_culprit_freq(full_dict, range_dict, band)
                        for band in self.band]
        self.culprit = [culprit_info[cul][0] for cul in culprit_info]
        self.description = [culprit_info[cul][1] for cul in culprit_info]
        # return dict(zip(self.band, culprit_info))

    # def chan_to_freq(self, chan, freq_vector):
    #     """
    #     OBSOLETE
    #     Returns the channel number where a given frequency is to be found.
    #     Frequency is in Hz.
    #     :param freq_vector:
    #     :param chan:
    #     :return:
    #     """
    #     self.peak_freq = freq_vector[chan]
    #
    # def grab_data(self, arr):
    #     """
    #     OBSOLETE
    #     :return:
    #     """
    #     num_samples = 50  # TODO hardcoded, choose 50 samples around event
    #     self.data = arr[:, self.peak_channel - num_samples: self.peak_channel + num_samples]
