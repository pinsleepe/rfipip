# Author: Monika Obrocka

import numpy as np


class RfiBlock(object):
    def __init__(self, data, times=None, freqs=None, tsamp=1.0):
        self.data = data
        self.channels = self.data.shape[1]
        # self.num_features = None
        self.nsamp = self.data.shape[0]
        if times is None:
            self.times = np.arange(self.nsamp, dtype=float)
        if freqs is None:
            self.freqs = np.arange(self.channels, dtype=float)
        self.tsamp = tsamp
        self.channelBw = self.freqs[1] - self.freqs[0]
