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

    def write_report(self, training_set, rfi_sam):
        with PdfPages('rfi_report.pdf') as pdf:
            f, _ = self.plot_corrupted(rfi_sam)
            pdf.savefig()
            f.close()
            #
            # plt.rc('text', usetex=False)
            # plt.figure(figsize=(5, 5))
            # labels = 'Good', 'RFI'
            # rfi_perc = rfi_sam * 100
            # good_perc = 100 - rfi_perc
            # sizes = [good_perc, rfi_perc]
            #
            # fig1, ax1 = plt.subplots()
            # ax1.pie(sizes,
            #         labels=labels,
            #         autopct='%1.1f%%',
            #         shadow=True,
            #         startangle=90)
            # ax1.axis('equal')
            # plt.title('Cleanliness of the observation')
            # pdf.savefig()  # saves the current figure into a pdf page
            # plt.close()

            plt.rc('text', usetex=False)
            plt.figure(figsize=(5, 5))
            training_set.groupby('culprit')['event'].nunique().plot.pie(autopct='%1.1f%%',
                                                                        labels=['unknown',
                                                                                'known'])
            plt.ylabel('')
            plt.title('Culprit classyfication')
            pdf.savefig()
            plt.close()

            plt.rc('text', usetex=False)
            fig = plt.figure(figsize=(5, 5))
            training_set.loc[lambda df: df.culprit > 0, :].groupby('description')['event'].nunique().plot.pie(autopct='%1.1f%%')
            plt.ylabel('')
            plt.title('Known culprit occurences')
            pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
            plt.close()

            plt.rc('text', usetex=False)
            fig = plt.figure(figsize=(5, 5))
            training_set.loc[lambda df: df.culprit > 0, :].groupby('description')['duration'].sum().plot.pie(autopct='%1.1f%%')
            plt.ylabel('')
            plt.title('Known culprit time occupancy')
            pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
            plt.close()

            # We can also set the file's metadata via the PdfPages object:
            d = pdf.infodict()
            d['Title'] = 'RFI Report'
            d['ModDate'] = datetime.datetime.today()
