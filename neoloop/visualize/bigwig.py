#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import pyBigWig
import numpy as np
from neoloop.visualize.base import TrackPlot, GenomeRange, get_logger

log = get_logger(__name__)

# coolbox.track.bigwig
class PlotBigWig(TrackPlot):

    DEFAULT_COLOR = "#33a02c"

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])
        if 'color' not in self.properties:
            self.properties['color'] = PlotBigWig.DEFAULT_COLOR
        if 'data_range_style' not in self.properties:
            # choices: {'text', 'y-axis'}
            # default 'text'
            # 'y-axis' style need set the .y_ax attribute
            self.properties['data_range_style'] = 'text'

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax

        genome_range = GenomeRange(chrom_region, start_region, end_region)
        log.debug("plotting {}".format(self.properties['file']))

        num_bins = self.__get_bins_num()
        self.__check_chrom_name(genome_range)
        scores_per_bin = self.__get_scores_per_bin(genome_range, num_bins)

        x_values = np.linspace(genome_range.start, genome_range.end, num_bins)

        if 'type' in self.properties and self.properties['type'] != 'fill':
            self.__plot_line_or_points(scores_per_bin, x_values)
        else:
            self.__plot_fill(scores_per_bin, x_values)

        ymin, ymax = self.__adjust_plot(genome_range)

        if "show_data_range" in self.properties and self.properties["show_data_range"] == 'no':
            pass
        else:
            self.genome_range = genome_range
            self.plot_data_range(ymin, ymax, self.properties['data_range_style'])

        self.plot_label()

        return self.ax

    def __get_bins_num(self, default_num=700):
        num_bins = default_num
        if 'number_of_bins' in self.properties:
            try:
                num_bins = int(self.properties['number_of_bins'])
            except TypeError:
                num_bins = default_num
                log.warning("'number_of_bins' value: {} for bigwig file {} "
                            "is not valid. Using default value (700)".format(self.properties['number_of_bins'],
                                                                             self.properties['file']))
        return num_bins

    def __get_scores_per_bin(self, genome_range, num_bins, max_try_nums=5):
        # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        # of the memory. This only occurs when calling trackPlot from different
        # processors. Reloading the file solves the problem.
        num_tries = 0
        while num_tries < max_try_nums:
            num_tries += 1
            try:
                scores_per_bin = np.array(self.bw.stats(genome_range.chrom, genome_range.start,
                                                        genome_range.end, nBins=num_bins)).astype(float)
            except Exception as e:
                import pyBigWig
                self.bw = pyBigWig.open(self.properties['file'])

                log.warning("error found while reading bigwig scores ({}).\nTrying again. Iter num: {}".
                            format(e, num_tries))
                pass
            else:
                if num_tries > 1:
                    log.warning("After {} the scores could be computed".format(num_tries))
                break
        return scores_per_bin

    def __check_chrom_name(self, genome_range):
        if genome_range.chrom not in self.bw.chroms().keys():
            genome_range.change_chrom_names()

        if genome_range.chrom not in self.bw.chroms().keys():
            log.warning("Can not read region {} from bigwig file:\n\n"
                        "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                        "and that the region is valid".format(str(genome_range), self.properties['file']))

    def __plot_line_or_points(self, scores_per_bin, x_values):
        if self.properties['type'].find(":") > 0:
            plot_type, size = self.properties['type'].split(":")
            try:
                size = float(size)
            except ValueError:
                raise ValueError("Invalid value: 'type = {}' in Track: {}\n"
                                 "A number was expected and found '{}'".format(
                    self.properties['type'],
                    self.properties['name'],
                    size))
        else:
            plot_type = self.properties['type']
            size = None

        if plot_type == 'line':
            self.ax.plot(x_values, scores_per_bin, '-', linewidth=size, color=self.properties['color'])

        elif plot_type == 'points':
            self.ax.plot(x_values, scores_per_bin, '.', markersize=size, color=self.properties['color'])

        else:
            raise ValueError("Invalid: 'type = {}' in Track: {}\n".format(self.properties['type'],
                                                                          self.properties['name']))

    def __plot_fill(self, scores_per_bin, x_values):
        if 'positive_color' not in self.properties or 'negative_color' not in self.properties:
            self.ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                                 color=self.properties['color'],
                                 facecolor=self.properties['color'],
                                 rasterized=True)
        else:
            self.ax.fill_between(x_values, 0, scores_per_bin, where=(scores_per_bin > 0),
                                 linewidth=0.1, color=self.properties['positive_color'],
                                 facecolor=self.properties['positive_color'],
                                 rasterized=True)
            self.ax.fill_between(x_values, 0, scores_per_bin, where=(scores_per_bin < 0),
                                 linewidth=0.1, color=self.properties['negative_color'],
                                 facecolor=self.properties['negative_color'],
                                 rasterized=True)

    def __adjust_plot(self, genome_range):
        self.ax.set_xlim(genome_range.start, genome_range.end)
        ymin, ymax = self.ax.get_ylim()
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            self.ax.set_ylim(ymax, ymin)
        else:
            self.ax.set_ylim(ymin, ymax)
        return ymin, ymax

    def plot_data_range(self, ymin, ymax, data_range_style):

        if data_range_style == 'text':
            assert hasattr(self, 'genome_range'), \
                "If use text style data range must, must set the .genome_range attribute"
            self.__plot_range_text(ymin, ymax)

        else:  # 'y-axis' style
            try:
                y_ax = self.y_ax
                self.plot_y_axis(y_ax)
            except AttributeError as e:
                log.exception(e)
                msg = "If use y-axis style data range must, must set the .y_ax attribute, switch to text style."
                log.warn(msg)
                self.plot_data_range(ymin, ymax, data_range_style='text')

    def __plot_range_text(self, ymin, ymax):
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
                     verticalalignment='top')

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


        