#cython: language_level=3
#cython: boundscheck=False
#cython: cdivision=True

import itertools, matplotlib, neoloop, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.gridspec import GridSpec
from neoloop.visualize.bed import Genes, plotGenes, Elements
from neoloop.visualize.bigwig import plotSignal, SigTrack
from neoloop.visualize.loops import Loops
from neoloop.assembly import complexSV
from neoloop.callers import Peakachu
from scipy import sparse

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

matplotlib.rcParams.update(new_rc_params)

class Triangle(Peakachu):

    track_colors = {
        'CNV': '#666666',
        'RNA+': '#E31A1C',
        'RNA-': '#E31A1C',
        'RNA': '#E31A1C',
        'H3K27ac': '#6A3D9A',
        'DNase': '#6A3D9A',
        '4C': '#33A02C'
    }

    chrom_colors = {'chr1': '#B15928',
            'chr2': '#6A3D9A',
            'chr3': '#CAB2D6',
            'chr4': '#FDBF6F',
            'chr5': '#A6CEE3',
            'chr6': '#B2DF8A',
            'chr7': '#1B9E77',
            'chr8': '#A6CEE3',
            'chr9': '#33A02C',
            'chr10':'#E6AB02',
            'chr11':'#E58606',
            'chr12':'#5D69B1',
            'chr13':'#FF7F00',
            'chr14':'#52BCA3',
            'chr15':'#99C945',
            'chr16':'#CC61B0',
            'chr17':'#88CCEE',
            'chr18':'#1D6996',
            'chr19':'#117733',
            'chr20':'#661100',
            'chr21':'#882255',
            'chr22':'#1F78B4',
            'chrX':'#666666',
            'chrY':'#5F4690'
            }

    def __init__(self, clr, assembly, span=None, protocol='insitu', slopes={},
        n_rows=3, track_partition=[4,0.1,0.3], space=0.04, correct='sweight', figsize=(7, 4)):
        """
        Recommended parameter settings.

        Set1: (Contact map + genes + chromosome bar)
            figsize=(7, 4)
            track_partition=[5, 0.2, 0.5]
            n_rows=3
            space=0.02
        
        Set2: (Contact map + genes + 1 signal track + chromosome bar)
            figsize=(7, 4.6)
            track_partition=[5, 0.2, 0.8, 0.5]
            n_rows=4
            genes: fontsize=7
            signal: data_range_size=7
            chromosome bar: name_size=7
        
        Set2_2: (Contact map + genes + cREs + chromosome bar)
            figsize=(7, 4.4)
            track_partition=[5, 0.3, 0.3, 0.5]
            n_rows=4
            genes: fontsize=7
            signal: data_range_size=7
            chromosome bar: name_size=7
        
        Set3: (Contact map + genes + 2 signal tracks + chromosome bar)
            figsize=(7, 5.2)
            track_partition=[5, 0.2, 0.8, 0.8, 0.5]
            n_rows=5
        
        Set4: (Contact map + genes + 3 signal tracks + chromosome bar)
            figsize=(7, 5.8)
            track_partition=[5, 0.2, 0.8, 0.8, 0.8, 0.5]
            n_rows=6
        
        Set5: (Contact map + genes + 4 signal tracks + chromosome bar)
            figsize=(7, 6.4)
            track_partition=[5, 0.2, 0.8, 0.8, 0.8, 0.8, 0.5]
            n_rows=7 
        
        Set6: (Contact map + genes + 5 signal tracks + chromosome bar)
            figsize=(7,7)
            track_partition=[5, 0.2, 0.8, 0.8, 0.8, 0.8, 0.8, 0.5]
            n_rows=8
        
        Set7: (Contact map + Arcs + 4C + genes + chromosome bar)
            figsize=(7, 5.2)
            track_partition=[5, 0.8, 0.8, 0.2, 0.5]
            n_rows=5

        """
        parse = assembly.split()
        self.assembly_ID = parse[0]
        if not span is None:
            assembly = complexSV(clr, '\t'.join(parse[1:]), col=correct, protocol=protocol,
                                 span=span, flexible=False, slopes=slopes)
        else:
            assembly = complexSV(clr, '\t'.join(parse[1:]), col=correct, protocol=protocol,
                                slopes=slopes)
        assembly.reorganize_matrix()
        assembly.correct_heterozygosity()
        self.assembly = assembly

        fig = plt.figure(figsize=figsize)
        self.fig = fig
        self.grid = GridSpec(n_rows, 1, figure=fig, left=0.1, right=0.9,
                    bottom=0.1, top=0.9, hspace=space, height_ratios=track_partition)
        self.track_count = 0
        self.w_h = figsize[0] / figsize[1]

        self.matrix = assembly.symm_hcm # should be a symmetric matrix
        self.bounds = assembly.bounds
        self.orients = assembly.orients
        self.Map = assembly.Map
        self.res = clr.binsize
        self.sig_data = {}
        self.sig_tracks = {}

        # define my colormap (traditional w --> r)
        self.cmap = LinearSegmentedColormap.from_list('interaction',
                ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])
    
    def print_coordinate(self, pos):

        if pos % 1000000 == 0:
            return '{0}M'.format(pos//1000000)
        else:
            return '{0:.2f}M'.format(pos/1000000)
    
    def matrix_plot(self, colormap='traditional', vmin=None, vmax=None, log=False,
        cbr_width=0.02, cbr_height=0.1, cbr_fontsize=7, no_colorbar=False):

        h_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        heatmap_pos = h_ax.get_position().bounds

        M = self.matrix
        n = M.shape[0]

        # Create the rotation matrix
        t = np.array([[1,0.5], [-1,0.5]])
        A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)

        if colormap=='traditional':
            cmap = self.cmap
        else:
            cmap = colormap
        
        # Plot the Heatmap ...
        x = A[:,1].reshape(n+1, n+1)
        y = A[:,0].reshape(n+1, n+1)
        y[y<0] = -y[y<0]

        if vmax is None:
            vmax = np.percentile(M[M.nonzero()], 95)
        if vmin is None:
            vmin = M.min()
        
        if log:
            vmin = M[np.nonzero(M)].min()
            vmax = M.max()
            sc = h_ax.pcolormesh(x, y, np.flipud(M), cmap=cmap,
                        edgecolor='none', snap=True, linewidth=.001,
                        norm=LogNorm(vmin, vmax), rasterized=True)
        else:
            sc = h_ax.pcolormesh(x, y, np.flipud(M), vmin=vmin, vmax=vmax, cmap=cmap,
                        edgecolor='none', snap=True, linewidth=.001, rasterized=True)
        
        h_ax.axis('off')
        self.heatmap_ax = h_ax
        self.hx = x
        self.hy = y
        
        # colorbar
        if not no_colorbar:
            c_ax = self.fig.add_axes([heatmap_pos[0]-0.02,
                                 (heatmap_pos[1]+0.9)/2,
                                 cbr_width,
                                 cbr_height])
            cbar = self.fig.colorbar(sc, cax=c_ax, ticks=[vmin, vmax], format='%.3g')
            c_ax.tick_params(labelsize=cbr_fontsize)
            self.cbar_ax = c_ax
        
    
    def plot_chromosome_bounds(self, line_color='k', linewidth=1.5, linestype=':'):

        n = self.matrix.shape[0]

        bounds = []
        for i in range(2, len(self.bounds), 2):
            bounds.append(self.bounds[i][0])
        
        for si in bounds:
            ei = n - 1
            x = [self.hx[:-1, :-1][n-1-si, si],
                 self.hx[:-1, :-1][n-1-si, ei]]
            y = [self.hy[:-1, :-1][n-1-si, si] - 1,
                 self.hy[:-1, :-1][n-1-si, ei] + 1]
            self.heatmap_ax.plot(x, y, color=line_color, linestyle=linestype,
                linewidth=linewidth)
        
        for ei in bounds:
            si = 0
            x = [self.hx[:-1, :-1][n-1-si, ei],
                 self.hx[:-1, :-1][n-1-ei, ei]]
            y = [self.hy[:-1, :-1][n-1-si, ei] + 1,
                 self.hy[:-1, :-1][n-1-ei, ei] - 1]
            self.heatmap_ax.plot(x, y, color=line_color, linestyle=linestype,
                linewidth=linewidth)
        
        self.heatmap_ax.set_xlim(self.hx.min(), self.hx.max())
        self.heatmap_ax.set_ylim(self.hy.min(), self.hy.max())
    
    def plot_neoTAD(self, ws=500000, color='#60636A'):

        from neoloop.tadtool.core import TADcaller
        import joblib

        hmm_folder = os.path.join(os.path.split(neoloop.__file__)[0], 'data')
        hmm = joblib.load(os.path.join(hmm_folder, 'HMM-model.pkl'))

        work = TADcaller(self.matrix, self.res, hmm, window_size=ws)
        work.callDomains()
        self.tads = work.loop_like()

        pairs = []
        for st, et, _ in self.tads:
            pairs.append((st+0.5, et+0.5))
            
        pairs.sort()

        tad_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        tad_ax.set_xlim(self.hx.min(), self.hx.max())
        bottoms = np.r_[[i%2 for i in range(len(pairs))]]
        widths = np.r_[[d[1]-d[0] for d in pairs]]
        lefts = np.r_[[d[0] for d in pairs]]
        tad_ax.barh(bottoms, widths, 0.8, lefts, color=color,
                    edgecolor='none')
        tad_ax.set_ylim(-0.5, 2.1)
        tad_ax.set_axis_off()
        self.tad_ax = tad_ax
    
    def plot_loops(self, loop_fil, marker_size=50, face_color='#1F78B4', edgecolors='#1F78B4', 
        marker_type='o', marker_alpha=1, cluster=True, onlyneo=False, filter_by_res=False):

        self.M = sparse.csr_matrix(self.matrix)
        self.r = self.res

        loops = Loops(loop_fil)
        self.loops = loops.loops[self.assembly_ID]
        
        Donuts = {}
        Bool = np.zeros(self.matrix.shape, dtype=bool)
        for l1, l2, label in self.loops:
            # Lodate the loop pixel at given resolution
            xs, xe = l1[1:]
            ys, ye = l2[1:]
            if filter_by_res:
                if xe - xs != self.res:
                    continue
            s_l = range(xs//self.res-1, int(np.ceil(xe/float(self.res)))+1)
            e_l = range(ys//self.res-1, int(np.ceil(ye/float(self.res)))+1)
            si, ei = None, None
            for i in s_l:
                for j in e_l:
                    loci1 = (l1[0], i*self.res)
                    loci2 = (l2[0], j*self.res)
                    if (loci1 in self.Map) and (loci2 in self.Map):
                        st, et = self.Map[loci1], self.Map[loci2]
                        if st > et:
                            st, et = et, st
                        if label:
                            if self.assembly.chains[st] == self.assembly.chains[et]:
                                continue
                        else:
                            if self.assembly.chains[st] != self.assembly.chains[et]:
                                continue
                        if si is None:
                            si, ei = st, et
                        else:
                            if self.matrix[st,et] > self.matrix[si,ei]:
                                si, ei = st, et
            if not si is None:
                if onlyneo:
                    if label:
                        Donuts[(si, ei)] = self.matrix[si, ei]
                else:
                    Donuts[(si, ei)] = self.matrix[si, ei]
        
        if cluster:
            Donuts = self.local_clustering(Donuts, min_count=1, r=20000, chains=self.assembly.chains)
        
        for si, ei in Donuts:
            Bool[si, ei] = 1
        
        lx = self.hx[:-1,:-1][np.flipud(Bool)]
        ly = self.hy[:-1,:-1][np.flipud(Bool)] + 1
        self.Bool = Bool

        if lx.size > 0:
            self.heatmap_ax.scatter(lx, ly, s=marker_size, c=face_color, marker=marker_type,
                edgecolors=edgecolors, alpha=marker_alpha)
        
        self.heatmap_ax.set_xlim(self.hx.min(), self.hx.max())
        self.heatmap_ax.set_ylim(self.hy.min(), self.hy.max())
    
    def plot_DI(self, ws=2000000, pos_color='#FB9A99', neg_color='#A6CEE3', y_axis_offset=0.01,
        data_range_size=7):

        from tadlib.hitad.chromLev import Chrom
        from scipy.sparse import csr_matrix
        from neoloop.visualize.bigwig import plot_y_axis

        di_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        
        hicdata = csr_matrix(self.matrix, shape=self.matrix.shape)
        tad_core = Chrom('chrN', self.res, hicdata, 'pseudo')
        tad_core._dw = min(tad_core._dw, hicdata.shape[0]-1)
        tad_core.windows = np.ones(hicdata.shape[0], dtype=int) * (ws // self.res)
        tad_core.calDI(tad_core.windows, 0)
        arr = tad_core.DIs

        pos_mask = arr >= 0
        neg_mask = arr < 0
        di_ax.fill_between(np.arange(arr.size), arr, where=pos_mask, color=pos_color,
                           edgecolor='none')
        di_ax.fill_between(np.arange(arr.size), arr, where=neg_mask, color=neg_color,
                           edgecolor='none')
        di_ax.set_xlim(self.hx.min(), self.hx.max())

        ymin, ymax = di_ax.get_ylim()
        ax_pos = di_ax.get_position().bounds
        y_ax = self.fig.add_axes([ax_pos[0]-y_axis_offset, ax_pos[1],
                                y_axis_offset, ax_pos[3]])
        plot_y_axis(y_ax, ymin, ymax, size=data_range_size)
        self.clear_frame(y_ax)
        self.clear_frame(di_ax)

        self.di_ax = di_ax
        self.DI_arr = arr


    def plot_arcs(self, h_ratio=0.3, arc_color='#386CB0', arc_alpha=1,
        lw=1.5, cutoff='bottom', species='human', release=97, gene_filter=[]):
        """
        Both plot_loops and plot_genes need to be called before this method.
        """
        from matplotlib.patches import Arc

        arc_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        arc_ax.set_xlim(self.hx.min(), self.hx.max())

        xs, ys = np.where(self.Bool)
        if not len(gene_filter):
            x, y = xs, ys
        else:
            # filter loops by genes
            refgenes = Genes(self.bounds, self.orients, self.res, species=species, release=release,
                         filter_=gene_filter).genes
            pos = set()
            for g in gene_filter:
                for t in refgenes:
                    if t[3] == g:
                        for i in range(t[1]//self.res, t[2]//self.res+1):
                            pos.add(i)
            x = []; y = []
            for i, j in zip(xs, ys):
                if (i in pos) or (j in pos):
                    x.append(i+0.5)
                    y.append(j+0.5)
        
        miny = maxy = 0
        for i, j in zip(x, y):
            mid = (i + j) / 2
            a, b = j - i, (j - i) * h_ratio
            if cutoff=='bottom':
                ty = b / 2 + 0.5
                if ty > maxy:
                    maxy = ty
                arc_ax.add_patch(Arc((mid, 0), a, b, 
                             theta1=0.0, theta2=180.0, edgecolor=arc_color, lw=lw, alpha=arc_alpha))
            else:
                ty = -b / 2 - 0.5
                if ty < miny:
                    miny = ty
                arc_ax.add_patch(Arc((mid, 0), a, b, 
                             theta1=180.0, theta2=360.0, edgecolor=arc_color, lw=lw, alpha=arc_alpha))
        
        arc_ax.set_ylim(miny, maxy)
        arc_ax.set_axis_off()
        self.arc_ax = arc_ax
        self.arc_x = x
        self.arc_y = y
            
    def plot_chromosome_bar(self, coord_ypos=0.3, name_ypos=0, coord_size=3.5, name_size=7,
        width=3, headwidth=8, remove_coord=False, color_by_order=[]):

        chrom_ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1
        chrom_colors = self.chrom_colors

        chrom_ax.set_xlim(self.hx.min(), self.hx.max())
        chrom_ax.set_ylim(0, 1)
        # chromosome bar
        n = self.matrix.shape[0]
        for i in range(0, len(self.bounds), 2):
            s = self.bounds[i][0]
            e = self.bounds[i+1][0]
            si = self.hx[:-1, :-1][n-1-s, s]
            ei = self.hx[:-1, :-1][n-1-e, e]
            o = self.orients[i//2]
            c = 'chr'+self.bounds[i][1][0].lstrip('chr')
            bar_color = chrom_colors[c]
            if len(color_by_order):
                bar_color = color_by_order[i//2]
            if o=='+':
                chrom_ax.annotate('', xy=(ei, 0.7), xytext=(si, 0.7),
                        xycoords='data', arrowprops=dict(color=bar_color, shrink=0, width=width, headwidth=headwidth),
                    )
            else:
                chrom_ax.annotate('', xy=(si, 0.7), xytext=(ei, 0.7),
                        xycoords='data', arrowprops=dict(color=bar_color, shrink=0, width=width, headwidth=headwidth),
                    )
            
            if not remove_coord:
                chrom_ax.text(si, coord_ypos, self.print_coordinate(self.bounds[i][1][1]),
                        ha='left', va='top', fontsize=coord_size)
                chrom_ax.text(ei, coord_ypos, self.print_coordinate(self.bounds[i+1][1][1]),
                        ha='right', va='top', fontsize=coord_size)
        
        # concatenate chromosome labels
        chroms = []
        for b in self.bounds:
            if not len(chroms):
                chroms.append([b[1][0], b[0], b[0]])
            else:
                if b[1][0] == chroms[-1][0]:
                    chroms[-1][-1] = b[0]
                else:
                    chroms.append([b[1][0], b[0], b[0]])
        
        unique_chrom_names = set([c[0] for c in chroms])
        for c in chroms:
            si = self.hx[:-1, :-1][n-1-c[1], c[1]]
            ei = self.hx[:-1, :-1][n-1-c[2], c[2]]
            if len(unique_chrom_names) > 1:
                chrom_ax.text((si+ei)/2, name_ypos, 'chr'+c[0].lstrip('chr'),
                        ha='center', va='top', fontsize=name_size)
            else:
                chrom_ax.text((si+ei)/2, name_ypos-0.15, 'chr'+c[0].lstrip('chr'),
                        ha='center', va='top', fontsize=name_size)
        
        chrom_ax.axis('off')
    
    def plot_genes(self, species='human', release=97, filter_=None,
        color='#999999', border_color='#999999', fontsize=7, labels='auto',
        style='flybase', global_max_row=False, label_aligns={}):

        genes = Genes(self.bounds, self.orients, self.res, species=species, release=release,
                filter_=filter_)
        _wk = plotGenes(genes.file_handler, color=color, fontsize=fontsize, labels=labels,
                style=style, global_max_row=global_max_row)
        
        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        _wk.plot(ax, 'chrN', 0, self.matrix.shape[0]*self.res, label_aligns)

        ax.set_xlim(-self.res/2, self.matrix.shape[0]*self.res+self.res/2)
        ax.axis('off')

        self.gene_ax = ax
        self.genes = genes
    
    def plot_cREs(self, bedfil, color='#E41A1C', alpha=0.5):

        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        eles = Elements(bedfil, self.bounds, self.orients, self.res)
        cREs = eles.cREs

        bottoms = np.zeros(len(cREs))
        widths = np.r_[[d[2]-d[1] for d in cREs]]
        lefts = np.r_[[d[1]+self.res/2 for d in cREs]]
        ax.barh(bottoms, widths, 0.8, lefts, color=color,
                edgecolor=color, alpha=alpha, linewidth=0)
        ax.set_ylim(-0.2, 1)
        ax.set_axis_off()
        ax.set_xlim(0, self.matrix.shape[0]*self.res+self.res)

        self.cRE_ax = ax
    
    def plot_motif(self, bedfil, plus_color='#444444', minus_color='#999999', ypos=0.5,
        marker_size=10, subset='+'):

        ax = self.fig.add_subplot(self.grid[self.track_count])
        self.track_count += 1

        eles = Elements(bedfil, self.bounds, self.orients, self.res)
        motifs = eles.cREs
        for m in motifs:
            x = (m[1] + m[2]) // 2
            strand = m[3]
            if strand != subset:
                continue
            if strand == '+':
                ax.scatter(x, ypos, s=marker_size, c=plus_color, marker='>')
            else:
                ax.scatter(x, ypos, s=marker_size, c=minus_color, marker='<')

        ax.set_ylim(0, 1)
        ax.set_axis_off()
        ax.set_xlim(0, self.matrix.shape[0]*self.res+self.res)
        self.motif_ax = ax
        
    
    def clear_frame(self, ax):

        for spine in ax.spines:
            ax.spines[spine].set_visible(False)
        
        ax.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
            labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    
    def plot_signal(self, track_name, bw_fil, factor=1, ax=None, color='auto', show_data_range=True,
        data_range_style='y-axis', data_range_size=7, max_value='auto', min_value='auto',
        y_axis_offset=0.01, show_label=True, label_size=7, label_pad=2, nBins=500,
        style_params={'type':'fill'}):
        '''
        Choices for data_range_style: ['y-axis', 'text'].
        '''
        sigs = SigTrack(bw_fil, self.bounds, self.orients, res=self.res, nBins=nBins,
                        multiply=factor)
        arr = sigs.stats('chrN', 0, self.matrix.shape[0]*self.res, 100)

        if color=='auto':
            if track_name in self.track_colors:
                color = self.track_colors[track_name]
            else:
                color = '#dfccde'

        if max_value=='auto':
            max_value = np.nanmax(arr)
        
        if min_value=='auto':
            min_value = 0

        _wk = plotSignal(sigs, color=color, show_data_range='no',
                        max_value=max_value, min_value=min_value,
                        number_of_bins=self.matrix.shape[0] * sigs.factor)
        
        if ax is None:
            ax = self.fig.add_subplot(self.grid[self.track_count])
            self.track_count += 1
        else:
            ax = ax

        _wk.properties.update(style_params)
        _wk.plot(ax, 'chrN', 0, self.matrix.shape[0] * self.res)
        self.clear_frame(ax)

        if show_data_range:
            if data_range_style=='y-axis':
                ax_pos = ax.get_position().bounds
                y_ax = self.fig.add_axes([ax_pos[0]-y_axis_offset, ax_pos[1],
                                        y_axis_offset, ax_pos[3]])
                _wk.plot_y_axis(y_ax, size=data_range_size)
                self.clear_frame(y_ax)
            else:
                _wk.plot_range_text(min_value, max_value, data_range_size)
        
        if show_label:
            if show_data_range and (data_range_style=='y-axis'):
                y_ax.set_ylabel(track_name, rotation=0, va='center', ha='right',
                    fontsize=label_size, labelpad=label_pad)
            else:
                ax.set_ylabel(track_name, rotation=0, va='center', ha='right', 
                    fontsize=label_size, labelpad=label_pad)
        
        ax.set_xlim(-self.res/2, self.matrix.shape[0]*self.res+self.res/2)

        self.sig_data[track_name] = sigs
        self.sig_tracks[track_name] = ax


    def outfig(self, outfile, dpi=200, bbox_inches='tight'):

        self.fig.savefig(outfile, dpi=dpi, bbox_inches=bbox_inches)
    
    def show(self):

        self.fig.show()

def plot_exp(exp, res, outfig, dpi=200):

    x = []
    y = []
    for i in sorted(exp):
        x.append(i*res)
        y.append(exp[i])
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    l1, = ax.plot(x, y)
    ax.set_xlabel('Genomic distance (bp)')
    ax.set_ylabel('Average contact frequency')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)

    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)

    ax.ticklabel_format(scilimits=(-3,3))

    plt.savefig(outfig, dpi=dpi, bbox_inches='tight')
    plt.close()
