#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import pyBigWig
import numpy as np
from coolbox.plots.track.bigwig import PlotBigWig
from coolbox.plots.track.base import TrackPlot

def plot_y_axis(y_ax, ymin, ymax, size):
    """
    Plot the scale of the y axis with respect to the plot_axis
    plot something that looks like this:
    ymax ┐
            │
            │
    ymin ┘
    Parameters
    ----------
    y_ax : matplotlib.axes.Axes
        Axis to use to plot the scale
    """
    def value_to_str(value):
        if value % 1 == 0:
            str_value = str(int(value))
        else:
            if value < 0.01:
                str_value = "{:.4f}".format(value)
            else:
                str_value = "{:.2f}".format(value)
        return str_value

    ymax_str = value_to_str(ymax)
    ymin_str = value_to_str(ymin)
    x_pos = [0, 0.5, 0.5, 0]
    y_pos = [0.01, 0.01, 0.99, 0.99]
    y_ax.plot(x_pos, y_pos, color='black', linewidth=1, transform=y_ax.transAxes)
    y_ax.text(-0.2, -0.01, ymin_str, verticalalignment='bottom', horizontalalignment='right',
                transform=y_ax.transAxes, fontsize=size)
    y_ax.text(-0.2, 1, ymax_str, verticalalignment='top', horizontalalignment='right',
                transform=y_ax.transAxes, fontsize=size)
    y_ax.patch.set_visible(False)

class SigTrack(object):
    """
    Map a reference signal track in .bigwig to a local assembly.
    """
    def __init__(self, file_, bounds, orients, res=10000, nBins=500, multiply=1,
        default_chrom='chrN'):

        self.factor = int(np.ceil(nBins / (bounds[-1][0]+1)))
        self.arr = np.zeros((bounds[-1][0]+1) * self.factor)
        self.multiply = multiply

        blocks = []
        for i in range(0, len(bounds), 2):
            b1 = bounds[i][1]
            b2 = bounds[i+1][1]
            if b1[1] > b2[1]:
                b1, b2 = b2, b1
            o = orients[i//2]
            si = bounds[i][0] * self.factor
            ei = (bounds[i+1][0] + 1) * self.factor
            blocks.append((si, ei, b1[0], b1[1], b2[1], o))
        self.blocks = blocks
        self.res = res
        self.chromsizes = {default_chrom: bounds[-1][0]*res}

        self.bw = pyBigWig.open(file_)


    def chroms(self):

        return self.chromsizes
    
    def stats(self, chrom, start, end, nBins):
        '''
        The input parameters doesn't mean too much, just to make a consistent interface with CoolBox.
        '''
        for si, ei, c, s, e, o in self.blocks:
            part = np.array(self.bw.stats(c, s, e, nBins=ei-si)).astype(float)
            if o == '+':
                self.arr[si:ei] = part
            else:
                self.arr[si:ei] = part[::-1]
        
        self.arr = self.arr * self.multiply

        return self.arr

class plotSignal(PlotBigWig):

    def __init__(self, bw, **kwargs):

        properties_dict = {
            'file': 'inMemory.bw',
            'signal': bw,
            'height': 3,
            'color': '#dfccde',
            'style': 'fill',
        }

        properties_dict.update(kwargs)
        properties_dict['type'] = properties_dict['style']

        self.properties = properties_dict
        self.bw = bw
    
    def plot_y_axis(self, y_ax, size):
        """
        Plot the scale of the y axis with respect to the plot_axis
        plot something that looks like this:
        ymax ┐
             │
             │
        ymin ┘
        Parameters
        ----------
        y_ax : matplotlib.axes.Axes
            Axis to use to plot the scale
        """
        plot_axis = self.ax
        ymin, ymax = plot_axis.get_ylim()

        plot_y_axis(y_ax, ymin, ymax, size)
    
    def plot_range_text(self, ymin, ymax, fontsize):

        genome_range = self.genome_range
        ydelta = ymax - ymin

        # set min max
        format_lim = lambda lim: int(lim) if float(lim) %1 == 0 else "{:.2f}".format(lim)
        ymax_print = format_lim(ymax)
        ymin_print = format_lim(ymin)
        small_x = 0.01 * genome_range.length
        # by default show the data range
        self.ax.text(genome_range.start - small_x, ymax - ydelta * 0.2,
                     "[ {} ~ {} ]".format(ymin_print, ymax_print),
                     horizontalalignment='left',
                     verticalalignment='top',
                     fontsize=fontsize)


        