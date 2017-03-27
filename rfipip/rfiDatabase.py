import numpy as np


class RfiDatabase(object):
    """

    This class is going to check in which frequency band does the maximum frequency on the spectrum
    correspond to, and will print out all the possible RFI candidates.

    """

    def __init__(self):
        self.dictionary = None

    def read_csv(self, csv_path):
        """

        :param csv_path:
        :return:
        """
        my_data = np.genfromtxt(csv_path,
                                delimiter=',',
                                dtype=[('dfloat', 'f8'),
                                       ('ufloat', 'f8'),
                                       ('string0', 'S50'),
                                       ('string1', 'S50')],
                                skip_header=2,
                                usecols=(7, 8, 9, 10))
        # if start freq is None delete row
        to_delete = []
        for row in range(my_data.shape[0]):
            if not my_data[row][0] > 0:
                to_delete.append(row)
        pruned_data = np.delete(my_data, to_delete, axis=0)
        return pruned_data

    def write_dict(self, file_list):
        """

        :param file_list:
        :return:
        """
        csvArray = []
        for file in file_list:
            csvArray.append(self.read_csv(file))

        # separate for bands and frequencies in a dictionary
        # find the band:
        rfiSpectrum = {}
        band = 0
        for arr in csvArray:
            for row in range(arr.shape[0]):
                # look for not empty fields
                if arr[row][1] > 0:
                    dict_str = str(arr[row][0]) + '-' + str(arr[row][1])
                    description = arr[row][2]
                    rfiSpectrum[dict_str] = {'description': description}
                    rfiSpectrum[dict_str]['band'] = band
                    band += 1
                    newIdx = row + 1
                    frequencies = {}
                    # look for non-empty freq-end
                    while newIdx < arr.shape[0] and not arr[newIdx][1] > 0:
                        frequencies[str(arr[newIdx][0])] = arr[newIdx][-1]
                        newIdx += 1
                rfiSpectrum[dict_str]['frequencies'] = frequencies
        self.dictionary = rfiSpectrum

