#!/usr/bin/env python
import Tkinter as tk
from PIL import Image, ImageTk
import os
from os.path import join, isfile
from os import listdir
import numpy as np


class Viewer:
    def __init__(self,
                 master,
                 filelist=None,
                 figurelist=None,
                 writePath=None):
        """

        :param master:
        :param filelist:
        :param figurelist: list with figures and paths
        """
        # TODO hardcoded file name
        if writePath:
            self.labels = join(writePath,
                               'labels.txt')
        else:
            self.labels = 'labels.txt'
        self.top = master
        if filelist:
            self.fileType = 0
        else:
            self.fileType = 1
        if filelist:
            self.files = filelist
            self.paths = filelist
        else:
            newfigurelist = []
            for fig in figurelist[0]:
                newfigurelist.append(self.fig2img(fig))
            self.files = newfigurelist
            self.paths = figurelist[1]
        self.index = 0
        # display first image
        if filelist:
            filename = filelist[0]
            if not os.path.exists(filename):
                print "Unable to find %s" % filename
                self.top.quit()
        else:
            filename = figurelist[0][0]
            chackPath = figurelist[1][0]
            if not os.path.exists(chackPath):
                print "Unable to find %s" % chackPath
                self.top.quit()
        if filelist:
            self.title = tk.Label(text=os.path.basename(filename))
        else:
            chackPath = figurelist[1][0]
            self.title = tk.Label(text=os.path.basename(chackPath))
        self.title.pack()

        if filelist:
            im = Image.open(filename)
        else:
            im = self.fig2img(filename)
        if im.format == "SPIDER":
            im = im.convert2byte()
        self.size = im.size
        self.tkimage = ImageTk.PhotoImage(im,
                                          palette=256)

        self.lbl = tk.Label(master,
                            image=self.tkimage)
        self.lbl.pack(side='top')

        # the button frame
        fr = tk.Frame(master)
        fr.pack(side='top',
                expand=1,
                fill='x')
        back = tk.Button(fr,
                         text="back",
                         command=lambda: self.nextframe(-1))
        back.grid(row=0,
                  column=0,
                  sticky="w",
                  padx=4,
                  pady=4)

        ilabel = tk.Label(fr,
                          text="image number:")
        ilabel.grid(row=0,
                    column=1,
                    sticky="e",
                    pady=4)

        self.evar = tk.IntVar()
        self.evar.set(1)
        entry = tk.Entry(fr,
                         textvariable=self.evar)
        entry.grid(row=0,
                   column=2,
                   sticky="w",
                   pady=4)
        entry.bind('<Return>',
                   self.getimgnum)

        # nextF = self.nextframe(1)
        next = tk.Button(fr,
                         text="next",
                         command=lambda: self.nextframe(1))
        next.grid(row=0,
                  column=3,
                  sticky="e",
                  padx=4,
                  pady=4)

        # Show the user feedback options:
        # 1. Interesting
        # 2. Maybe
        # 3. Uninteresting

        # TODO hardcoded paths
        self.thumbsup = ImageTk.PhotoImage(Image.open('/home/monika/Transcend/katai_home/python_projects/scripts/up.png'))
        self.thumbsdown = ImageTk.PhotoImage(Image.open('/home/monika/Transcend/katai_home/python_projects/scripts/down.png'))
        self.interButton = tk.Button(fr,
                                     text='Interesting',
                                     compound=tk.LEFT,
                                     image=self.thumbsup,
                                     command=self.chooseInteresting)
        self.interButton.grid(row=1,
                              column=0)
        self.maybeButton = tk.Button(fr,
                                     text='Maybe',
                                     command=self.chooseMaybe)
        self.maybeButton.grid(row=1,
                              column=1)
        self.unintButton = tk.Button(fr,
                                     text='Uninteresting',
                                     compound=tk.RIGHT,
                                     image=self.thumbsdown,
                                     command=self.chooseUninteresting)
        self.unintButton.grid(row=1,
                              column=2)

        # Create the quit button
        self.quitButton = tk.Button(fr,
                                    text='Quit',
                                    command=self.top.quit,
                                    background='#fcc')
        self.quitButton.grid(ipadx=30,
                             pady=10)

    def getImage(self, filename):
        """

        :param filename:
        :return:
        """
        if self.fileType == 0:
            im = Image.open(filename)
        else:
            im = filename
        if im.format == "SPIDER":
            im = im.convert2byte()
        if im.size != self.size:
            print "all images must be same dimensions:"
            f1 = os.path.basename(self.files[0])
            f2 = os.path.basename(filename)
            print "%s: %s, %s : %s" % (f1, str(self.size), f2, str(im.size))
            self.top.quit()
        return im

    def nextframe(self,
                  i=1,
                  imgnum=-1):
        """

        :param i:
        :param imgnum:
        :return:
        """
        if imgnum == -1:
            self.index += i
        else:
            self.index = imgnum - 1
        if self.index >= len(self.files):
            self.index = 0
        elif self.index < 0:
            self.index = len(self.files) - 1
        filename = self.files[self.index]
        if self.fileType == 0:
            if not os.path.exists(filename):
                print "Unable to find %s" % filename
                self.top.quit()
            self.title.configure(text=os.path.basename(filename))
        self.evar.set(self.index + 1)

        im = self.getImage(filename)
        self.tkimage.paste(im)

    def getimgnum(self):
        """

        :param event:
        :return:
        """
        self.nextframe(imgnum=self.evar.get())

    def chooseInteresting(self):
        string = self.paths[self.index] + ' 1' + '\n'
        self.append2file(string)

    def chooseMaybe(self):
        string = self.paths[self.index] + ' 2' + '\n'
        self.append2file(string)

    def chooseUninteresting(self):
        string = self.paths[self.index] + ' 0' + '\n'
        self.append2file(string)

    # taken from http://www.icare.univ-lille1.fr/node/1141
    def fig2data(self, fig):
        """
        @brief Convert a Matplotlib figure to a 4D numpy array 
        with RGBA channels and return it
        @param fig a matplotlib figure
        @return a numpy 3D array of RGBA values
        """
        # draw the renderer
        fig.canvas.draw()

        # Get the RGBA buffer from the figure
        w, h = fig.canvas.get_width_height()
        buf = np.fromstring(fig.canvas.tostring_argb(), dtype=np.uint8)
        buf.shape = (w, h, 4)

        # canvas.tostring_argb give pixmap in ARGB mode.
        #  Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll(buf, 3, axis=2)
        return buf

    def fig2img(self, fig):
        """
        @brief Convert a Matplotlib figure to a PIL Image 
        in RGBA format and return it
        @param fig a matplotlib figure
        @return a Python Imaging Library ( PIL ) image
        """
        # put the figure pixmap into a numpy array
        buf = self.fig2data(fig)
        w, h, d = buf.shape
        return Image.frombytes("RGBA", (w, h), buf.tostring())

    def append2file(self, string):
        with open(self.labels, "a") as myfile:
            myfile.write(string)

# --------------------------------------------------------------------
if __name__ == "__main__":

    path = '/home/monika/Transcend/katai_home/python_projects/scripts/pngs'
    filelist = [os.path.join(path, f) for f in listdir(path) if isfile(join(path, f))]

    root = tk.Tk()
    app = Viewer(root, filelist)
    root.mainloop()
