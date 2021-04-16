from coolbox.core.track import Bed
from coolbox.plots.track import PlotBed
from coolbox.utilities import ReadBed, Interval, IntervalTree
from pyensembl import EnsemblRelease
import io, logging, bisect

log = logging.getLogger(__name__)

def get_blocks(bounds, orients, res):

    blocks = []
    for i in range(0, len(bounds), 2):
        b1 = bounds[i][1]
        b2 = bounds[i+1][1]
        if b1[1] > b2[1]:
            b1, b2 = b2, b1
        o = orients[i//2]
        s = bounds[i][0] * res
        blocks.append((s, b1[0].lstrip('chr'), b1[1], b2[1], o))
    
    return blocks


class Elements(object):
    '''
    Map cis regulatory elements to a new assembly.
    '''
    def __init__(self, bedfil, bounds, orients, res=10000, default_chrom='chrN'):

        self.blocks = get_blocks(bounds, orients, res)
        self.chrom_name = default_chrom
        self.elements = self.parse_bed(bedfil)
        self.map_elements()
    
    def map_elements(self):

        pool = set()
        for p, c, s, e, o in self.blocks:
            eles = self.extract_by_region('chr'+c, int(s), int(e))
            block_bounds = [p, p+e-s]
            for g in eles:
                if o == '-':
                    start = p + e - g[1]
                    end = p + e - g[2]
                else:
                    start = p + g[1] - s
                    end = p + g[2] - s

                if start > end:
                    start, end = end, start
                
                # truncate genes if it's beyond the original chromosome bounds
                if start < block_bounds[0]:
                    start = block_bounds[0]
                if end > block_bounds[1]:
                    end = block_bounds[1]
                
                if end <= start:
                    continue
                
                if len(g) > 3:
                    extra = list(g[3:])
                    if extra[0] == '+':
                        if o == '-':
                            extra[0] = '-'
                    elif extra[0] == '-':
                        if o == '-':
                            extra[0] = '+'
                    extra = tuple(extra)
                else:
                    extra = tuple()

                record = (self.chrom_name, start, end) + extra
                pool.add(record)
        
        self.cREs = sorted(pool)

    
    def parse_bed(self, bedfil):

        D =  {}
        with open(bedfil, 'r') as source:
            for line in source:
                p = line.rstrip().split()
                if not p[0] in D:
                    D[p[0]] = []
                D[p[0]].append([int(p[1]), int(p[2])] + p[3:])
        
        for c in D:
            D[c].sort()
        
        return D
    
    def extract_by_region(self, chrom, start, end):

        cache = set()
        List = self.elements[chrom]
        idx = max(0, bisect.bisect(List, [start, end])-1)
        for q in List[idx:]:
            if q[1] <= start:
                continue
            if q[0] >= end:
                break
            
            cache.add((chrom,) + tuple(q))
        
        return cache

class Genes(object):
    """
    Map reference genes onto a new assembly. 
    """
    def __init__(self, bounds, orients, res=10000, default_chrom='chrN',
        species='human', release=97, filter_=None):
        
        ref = EnsemblRelease(release, species=species)
        ref.download()
        ref.index()
        self.ref = ref
        if not filter_ is None:
            # a pre-defined gene list, such as cancer-related genes
            self.filter_list = set([g.upper() for g in filter_])
        else:
            self.filter_list = None

        self.blocks = get_blocks(bounds, orients, res)

        self.chrom_name = default_chrom

        self.map_genes()
        self.file_handler = io.StringIO('\n'.join(['\t'.join(list(map(str, g))) for g in self.genes]))
    
    def _check_overlap(self, i1, i2):

        if i1[0] != i2[0]:
            return 0, 0
        
        t1 = sorted(i1[1:])
        t2 = sorted(i2[1:])

        if (t1[1]<=t2[0]) or (t2[1]<=t1[0]):
            return 0, 0
        
        mi = t1 + t2
        mi.sort()
        OR = (mi[2] - mi[1]) / (mi[3] - mi[0]) # intersect / union

        return OR, mi[2]-mi[1]
    
    def map_genes(self):

        pool = []
        for p, c, s, e, o in self.blocks:
            genes = self.ref.genes_at_locus(c, int(s), int(e))
            block_bounds = [p, p+e-s]
            for g in genes:
                name = g.gene_name
                if not self.filter_list is None:
                    if not name.upper() in self.filter_list:
                        continue

                if not g.strand in ['+', '-']:
                    continue
                # gene body
                if o == '-':
                    start = p + e - g.start
                    end = p + e - g.end
                else:
                    start = p + g.start - s
                    end = p + g.end - s

                if start > end:
                    start, end = end, start
                
                # truncate genes if it's beyond the original chromosome bounds
                if start < block_bounds[0]:
                    start = block_bounds[0]
                if end > block_bounds[1]:
                    end = block_bounds[1]
                
                if end <= start:
                    continue

                if o == '+':
                    strand = g.strand
                else:
                    if g.strand == '+':
                        strand = '-'
                    else:
                        strand = '+'
                
                if len(g.exons):
                    # collect exons
                    exons = []
                    for en in g.exons:
                        size = en.end - en.start
                        if o == '-':
                            e_st = p + e - en.start
                            e_ed = p + e - en.end
                        else:
                            e_st = p + en.start - s
                            e_ed = p + en.end - s
                        if e_st > e_ed:
                            e_st, e_ed = e_ed, e_st
                        
                        # truncate exons
                        if e_st < start:
                            e_st = start
                        if e_ed > end:
                            e_ed = end
                        
                        if e_ed <= e_st:
                            continue

                        if e_st - start > 0:
                            exons.append((e_st-start, size)) # start relative to gene start
                    exons.sort()
                    exons = self._remove_redundant_exons(exons)
                    starts = ','.join([str(en[0]) for en in exons])+','
                    sizes = ','.join([str(en[1]) for en in exons])+','
                    exon_num = len(exons)
                else:
                    exon_num = 1
                    starts = '{0},'.format(0)
                    sizes = '{0},'.format(end-start)

                record = [self.chrom_name, start, end, name, 0, strand, start, end,
                        0, exon_num, sizes, starts]
                pool.append(record)
        
        self.genes = self._remove_redundant_genes(pool)
    
    def _remove_redundant_genes(self, genes):

        dedup = {g[3]:[0,0,0] for g in genes}
        for g in genes:
            Len = g[2] - g[1]
            cL = dedup[g[3]][2] - dedup[g[3]][1]
            if Len > cL:
                dedup[g[3]] = g
        
        genes = sorted([dedup[g] for g in dedup])

        return genes
    
    def _remove_redundant_exons(self, exons):

        pre_tree = {}
        for e in exons:
            if not e[0] in pre_tree:
                pre_tree[e[0]] = []
            pre_tree[e[0]].append(e[1])
        
        pre_exon = {e:max(pre_tree[e]) for e in pre_tree}

        nondup = []
        cur = 0
        for e in sorted(pre_exon):
            if e < cur:
                continue
            nondup.append((e, pre_exon[e]))
            cur = e + pre_exon[e]
        
        return nondup


class plotGenes(PlotBed):

    DEFAULT_COLOR = "#1f78b4"

    def __init__(self, file_, **kwargs):

        # here, file_ is an opened file or io.StringIO
        properties_dict = {
            'file_': file_,
            'file': 'inMemory.bed',
            'height': Bed.DEFAULT_HEIGHT,
            'color': 'bed_rgb',
            'border_color': 'black',
            'fontsize': Bed.DEFAULT_FONTSIZE,
            'title': '',
            'labels': 'auto',
            'style': 'flybase',
            'display': 'stacked',
            'global_max_row': False,
        }

        properties_dict.update(kwargs)

        labels = properties_dict.get('labels')
        if labels == 'auto':
            properties_dict['labels'] = 'auto'
        elif labels is True:
            properties_dict['labels'] = 'on'
        else:
            properties_dict['labels'] = 'off'
        
        self.properties = properties_dict
        
        self.bed_type = None  # once the bed file is read, this is bed3, bed6 or bed12
        self.len_w = None  # this is the length of the letter 'w' given the font size
        self.interval_tree = {}  # interval tree of the bed regions

        from matplotlib import font_manager

        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])

        if 'interval_height' not in self.properties:
            self.properties['interval_height'] = 100

        if self.properties['labels'] != 'on':
            self.is_draw_labels = False
        else:
            self.is_draw_labels = True

        self.colormap = None

        # to set the distance between rows
        self.row_scale = self.properties['interval_height'] * 2.3

        self.interval_tree, min_score, max_score = self._process_bed()
    

    def _process_bed(self):

        bed_file_h = ReadBed(self.properties['file_'])
        self.bed_type = bed_file_h.file_type

        if 'color' in self.properties and self.properties['color'] == 'bed_rgb' and \
                self.bed_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(PlotBed.DEFAULT_COLOR))
            self.properties['color'] = PlotBed.DEFAULT_COLOR

        valid_intervals = 0
        interval_tree = {}

        max_score = float('-inf')
        min_score = float('inf')
        for bed in bed_file_h:
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            if bed.chromosome not in interval_tree:
                interval_tree[bed.chromosome] = IntervalTree()

            interval_tree[bed.chromosome].add(Interval(bed.start, bed.end, bed))
            valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file_name']))

        return interval_tree, min_score, max_score
    
    def _get_length_w(self, fig_width, region_start, region_end):
        '''
        to improve the visualization of the genes it is good to have an estimation of the label
        length. In the following code I try to get the length of a 'W' in base pairs.
        '''
        if self.is_draw_labels:
            # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
            inches_per_pt = 1.0 / 72.27
            font_in_inches = self.properties['fontsize'] * inches_per_pt
            region_len = region_end - region_start
            bp_per_inch = region_len / fig_width
            font_in_bp = font_in_inches * bp_per_inch
            self.len_w = font_in_bp
            log.debug("len of w set to: {} bp".format(self.len_w))
        else:
            self.len_w = 1

        return self.len_w
    
    def _get_y_pos(self, free_row):
        """
        The y_pos is set such that regions to be plotted do not overlap (stacked). To override this
        the properties['collapsed'] needs to be set.
        The algorithm uses a interval tree (self.region_interval) to check the overlaps
        and a sort of coverage vector 'rows used' to identify the row in which to plot
        Return
        ------
        ypos : int
            y position.
        """
        # if the domain directive is given, ypos simply oscilates between 0 and 100
        if self.properties['display'] == 'interlaced':
            ypos = self.properties['interval_height'] if self.counter % 2 == 0 else 1

        elif self.properties['display'] == 'collapsed':
            ypos = 0

        else:
            ypos = free_row * self.row_scale
        return ypos
    
    
    def plot(self, ax, chrom_region, start_region, end_region, label_align={}):
        self.counter = 0
        self.small_relative = 0.004 * (end_region - start_region)
        self._get_length_w(ax.get_figure().get_figwidth(), start_region, end_region)

        genes_overlap = sorted(self.interval_tree[chrom_region][start_region:end_region])

        if self.properties['labels'] == 'auto':
            if len(genes_overlap) > 30:
                # turn labels off when too many intervals are visible.
                self.is_draw_labels = False
            else:
                self.is_draw_labels = True

        max_num_row_local = 1
        max_ypos = 0

        row_last_position = []
        for region in genes_overlap:
            
            self.counter += 1
            bed = region.data

            if self.is_draw_labels:
                num_name_characters = len(bed.name) + 2  # +2 to account for an space before and after the name
                bed_extended_end = int(bed.end + (num_name_characters * self.len_w))
            else:
                bed_extended_end = (bed.end + 2 * self.small_relative)

            # get smallest free row
            if len(row_last_position) == 0:
                free_row = 0
                row_last_position.append(bed_extended_end)
            else:
                # get list of rows that are less than bed.start, then take the min
                idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                if len(idx_list):
                    free_row = min(idx_list)
                    row_last_position[free_row] = bed_extended_end
                else:
                    free_row = len(row_last_position)
                    row_last_position.append(bed_extended_end)

            rgb, edgecolor = self.get_rgb_and_edge_color(bed)

            ypos = self._get_y_pos(free_row)

            # do not plot if the maximum interval rows to plot is reached
            if 'gene_rows' in self.properties and free_row >= int(self.properties['gene_rows']):
                continue

            if free_row > max_num_row_local:
                max_num_row_local = free_row
            if ypos > max_ypos:
                max_ypos = ypos

            if self.properties['style'] == 'flybase':
                self._draw_gene_with_introns_flybase_style(ax, bed, ypos, rgb, edgecolor)
            else:
                self._draw_gene_simple(ax, bed, ypos, rgb, edgecolor)

            if not self.is_draw_labels:
                pass
            elif bed.end > start_region and bed.start < end_region:
                if bed.name in label_align:
                    if label_align[bed.name] == 'left':
                        ax.text(bed.end + self.small_relative, ypos + (float(self.properties['interval_height']) / 2),
                                bed.name, horizontalalignment='left',
                                verticalalignment='center', fontproperties=self.fp)
                    else:
                        ax.text(bed.start - self.small_relative, ypos + (float(self.properties['interval_height']) / 2),
                                bed.name, horizontalalignment='right',
                                verticalalignment='center', fontproperties=self.fp)
                else:
                    ax.text(bed.end + self.small_relative, ypos + (float(self.properties['interval_height']) / 2),
                                bed.name, horizontalalignment='left',
                                verticalalignment='center', fontproperties=self.fp)

        if self.counter == 0:
            log.warning("*Warning* No intervals were found for file {} "
                        "in Track '{}' for the interval plotted ({}:{}-{}).\n".
                        format(self.properties['file'], self.properties['name'], chrom_region, start_region, end_region))

        ymax = 0

        if 'global_max_row' in self.properties and self.properties['global_max_row'] == 'yes':
            ymin = self.max_num_row[chrom_region] * self.row_scale

        elif 'gene_rows' in self.properties:
            ymin = int(self.properties['gene_rows']) * self.row_scale
        else:
            ymin = max_ypos + self.properties['interval_height']

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

        if 'display' in self.properties:
            if self.properties['display'] == 'domain':
                ax.set_ylim(-5, 205)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(-5, 105)

        ax.set_xlim(start_region, end_region)

        self.plot_label()
    
    def _draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws an interval with direction (if given)
        """
        from matplotlib.patches import Polygon

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos), bed.end - bed.start, self.properties['interval_height'],
                                   edgecolor=edgecolor, facecolor=rgb, linewidth=0.5))
        else:
            vertices = self._draw_arrow(ax, bed.start, bed.end, bed.strand, ypos)
            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=rgb,
                                 linewidth=0.5))
    
    def _draw_arrow(self, ax, start, end, strand, ypos):

        half_height = float(self.properties['interval_height']) / 2
        if strand == '+':
            x0 = start
            x1 = end  # - self.small_relative
            y0 = ypos
            y1 = ypos + self.properties['interval_height']

            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.
            -----------------\
            ---------------- /
            """

            vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + self.small_relative, y0 + half_height), (x1, y0)]

        else:
            x0 = start  # + self.small_relative
            x1 = end
            y0 = ypos
            y1 = ypos + self.properties['interval_height']

            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.
            /-----------------
            \-----------------
            """

            vertices = [(x0, y0), (x0 - self.small_relative, y0 + half_height), (x0, y1), (x1, y1), (x1, y0)]

        return vertices
    
    def _draw_gene_with_introns_flybase_style(self, ax, bed, ypos, rgb, edgecolor):
        """
        draws a gene using different styles
        """
        from matplotlib.patches import Polygon
        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self._draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = float(self.properties['interval_height']) / 2
        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=0.5, zorder=-1)

        # get start, end of all the blocks
        positions = []
        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x0 < bed.thick_start < x1:
                positions.append((x0, bed.thick_start, 'UTR'))
                positions.append((bed.thick_start, x1, 'coding'))

            elif x0 < bed.thick_end < x1:
                positions.append((x0, bed.thick_end, 'coding'))
                positions.append((bed.thick_end, x1, 'UTR'))

            else:
                if x1 < bed.thick_start or x0 > bed.thick_end:
                    type = 'UTR'
                else:
                    type = 'coding'

                positions.append((x0, x1, type))

        # plot all blocks as rectangles except the last if the strand is + or
        # the first is the strand is -, which are drawn as arrows.
        if bed.strand == '-':
            positions = positions[::-1]

        first_pos = positions.pop()
        if first_pos[2] == 'UTR':
            _rgb = 'grey'
        else:
            _rgb = rgb

        vertices = self._draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=0.5))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = 'grey'
            else:
                _rgb = rgb
            vertices = [(start_pos, ypos), (start_pos, ypos + self.properties['interval_height']),
                        (end_pos, ypos + self.properties['interval_height']), (end_pos, ypos)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=0.5))

    

    