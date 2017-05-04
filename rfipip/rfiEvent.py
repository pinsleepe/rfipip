
class RfiEvent(object):
    """

    """
    def __init__(self,
                 mode=0,
                 peak_channel=None,
                 freq_vector=None,
                 arr_data=None):

        self.time_occupancy = None
        self.peak_channel = peak_channel
        self.peak_freq = None
        self.bandwidth = None
        self.mode = mode  # not used atm
        self.data = None

        self.init(freq_vector, arr_data)

    def init(self, freq_vector, arr_data):
        """
        
        :param freq_vector:
        :param arr_data:
        :return:
        """
        return

    def chan_to_freq(self, chan, freq_vector):
        """
        OBSOLETE
        Returns the channel number where a given frequency is to be found.
        Frequency is in Hz.
        :param freq_vector:
        :param chan:
        :return:
        """
        self.peak_freq = freq_vector[chan]

    def grab_data(self, arr):
        """
        OBSOLETE
        :return:
        """
        num_samples = 50  # TODO hardcoded, choose 50 samples around event
        self.data = arr[:, self.peak_channel - num_samples: self.peak_channel + num_samples]
