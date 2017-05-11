import datetime
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


class RfiReport(object):

    def __init__(self, report_name,
                 report_path=None):
        self.name = report_name
        self.path = None
        if report_path:
            self.path = report_path

    def plot_corrupted(self, rfi_sam):
        """
        
        :param rfi_sam: 
        :return: 
        """
        labels = 'Good', 'RFI'
        rfi_perc = rfi_sam * 100
        good_perc = 100 - rfi_perc
        sizes = [good_perc, rfi_perc]

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.pie(sizes,
               labels=labels,
               autopct='%1.1f%%',
               shadow=True,
               startangle=90)
        ax.axis('equal')
        plt.title('Cleanliness of the observation')
        return fig, ax

    def plot_culprit_classyfication(self, training_set):
        """
        
        :param training_set: 
        :return: 
        """
        fig, ax = plt.subplots(figsize=(5, 5))
        training_set.groupby('culprit')['event'].nunique().plot.pie(autopct='%1.1f%%',
                                                                    labels=['unknown',
                                                                            'known'])
        plt.ylabel('')
        plt.title('Culprit classification')
        return fig, ax

    def plot_culprit_occurences(self, training_set):
        """
        
        :param training_set: 
        :return: 
        """
        fig, ax = plt.subplots(figsize=(5, 5))
        training_set.loc[lambda df: df.culprit > 0, :].groupby('description')['event'].nunique().plot.pie(
            autopct='%1.1f%%')
        plt.ylabel('')
        plt.title('Known culprit occurences')
        return fig, ax

    def plot_culprit_time_occupancy(self, training_set):
        """
        
        :param training_set: 
        :return: 
        """
        fig, ax = plt.subplots(figsize=(5, 5))
        training_set.loc[lambda df: df.culprit > 0, :].groupby('description')['duration'].sum().plot.pie(
            autopct='%1.1f%%')
        plt.ylabel('')
        plt.title('Known culprit time occupancy')
        return fig, ax

    def write_report(self, training_set,
                     rfi_sam):
        """
        
        :param training_set: 
        :param rfi_sam: 
        :return: 
        """
        with PdfPages(self.name) as pdf:

            f, _ = self.plot_corrupted(rfi_sam)
            pdf.savefig()

            f, _ = self.plot_culprit_classyfication(training_set)
            pdf.savefig()

            f, _ = self.plot_culprit_occurences(training_set)
            pdf.savefig()

            f, _ = self.plot_culprit_time_occupancy(training_set)
            pdf.savefig()

            # We can also set the file's metadata via the PdfPages object:
            d = pdf.infodict()
            d['Title'] = 'RFI Report'
            d['ModDate'] = datetime.datetime.today()
