.. note:: the reference data for *calculate-cnv* have been migrated to dropbox,
   if you encounter errors when you run *calculate-cnv*, consider upgrading your
   package by running ``pip install -U neoloop``.

Neo-loop Finder
***************
.. image:: https://codeocean.com/codeocean-assets/badge/open-in-code-ocean.svg
   :target: https://codeocean.com/capsule/8407443/tree/v1
.. image:: https://static.pepy.tech/personalized-badge/neoloop?period=total&units=international_system&left_color=black&right_color=orange&left_text=Downloads
   :target: https://pepy.tech/project/neoloop

Although recent efforts have shown that structural variations (SVs) can disrupt the 3D genome organization and induce enhancer-hijacking, no computational tools exist to detect such events from chromatin interaction data, such as Hi-C. Here, we develop NeoLoopFinder, a computational framework to identify the chromatin interactions induced by SVs, such as inter-chromosomal translocations, large deletions, and inversions. Our framework can automatically reconstruct local Hi-C maps surrounding the breakpoints, normalize copy number variation and allele effects, and capture local optimal signals. We applied NeoLoopFinder in Hi-C data from 50 cancer cell lines and primary tumors and identified tens of recurrent genes associated with enhancer-hijacking in different cancer types. To validate the algorithm, we deleted hijacked enhancers by CRISPR/Cas9 and showed that the deletions resulted in the reduction of the target oncogene expression. In summary, NeoLoopFinder is a novel tool for identifying potential tumorigenic mechanisms and suggesting new diagnostic and therapeutic targets.

Citation
========
Wang, X., Xu, J., Zhang, B., Hou, Y., Song, F., Lyu, H., Yue, F. Genome-wide detection of
enhancer-hijacking events from chromatin interaction data in re-arranged genomes. Nat Methods. 2021.


Installation
============
NeoLoopFinder and all the dependencies can be installed using `conda <https://conda.io/miniconda.html>`_
or `pip <https://pypi.org/project/pip/>`_::

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda create -n neoloop python=3.7.1 cython=0.29.13 cooler=0.8.6 numpy=1.17.2 scipy=1.3.1 joblib=0.13.2 scikit-learn=0.20.2 networkx=1.11 pyensembl=1.8.0 matplotlib=3.1.1 pybigwig=0.3.17 pomegranate=0.10.0
    $ conda activate neoloop
    $ conda install -c r r=3.5.1 rpy2=2.9.4 r-mgcv=1.8_23
    $ pip install neoloop TADLib==0.4.2 coolbox==0.1.7

Overview
========
neoloop-finder is distributed with 8 scripts. You can learn the basic usage of each script
by typing ``command [-h]`` in a terminal window, where "command" is one of the following
script names:

- calculate-cnv

  Calculate the copy number variation profile from Hi-C map using a generalized additive model with the Poisson link function

- segment-cnv

  Perform HMM segmentation on a pre-calculated copy number variation profile.

- plot-cnv
  
  Plot genome-wide CNV profiles and segments.

- correct-cnv

  Remove copy number variation effects from cancer Hi-C.

- simulate-cnv

  Simulate CNV effects on a normal Hi-C. The inputs are the Hi-C matrix of a normal cell in .cool format,
  the Hi-C matrix of a cancer cell in .cool format, and the CNV segmentation file of the same cancer cell
  in bedGraph format.

- assemble-complexSVs

  Assemble complex SVs. The inputs are a list of simple SVs and the Hi-C matrix of the same sample.

- neoloop-caller

  Identify neo-loops across SV breakpoints. The required inputs are the output SV assemblies from
  ``assemble-complexSVs`` and the corresponding Hi-C map in .cool format.

- neotad-caller

  Identify neo-TADs. The inputs are the same as ``neoloop-caller``.

- searchSVbyGene

  Search SV assemblies by gene name.

CNV normalization
=================
As copy number variations (CNVs) can distort Hi-C signals in cancer cells, we proposed a modified
matrix balancing algorithm to remove such effects along with other systematic biases including mappability,
GC content, and restriction fragment sizes. In our implementation, you can perform this CNV normalization by
sequentially running ``calculate-cnv``, ``segment-cnv``, and ``correct-cnv``. The Hi-C map in
`.cool <https://github.com/open2c/cooler>`_ format is the only required input to this pipeline, and the
bias vector returned by this algorithm will be stored in the "sweight" column in the `bins <https://cooler.readthedocs.io/en/latest/datamodel.html#bins>`_
table of the cool file.

By default, ``assemble-complexSVs``, ``neoloop-caller``, and ``neotad-caller`` will use the "sweight" column to
normalize the Hi-C matrix. However, you can change this option to ICE normalization by specifying ``--balance-type ICE``.

.. note:: if your .cool files were transformed from .hic files, please make sure to add the "chr" prefix to your cool files using `add_prefix_to_cool.py <https://raw.githubusercontent.com/XiaoTaoWang/NeoLoopFinder/master/scripts/add_prefix_to_cool.py>`_ before your run ``calculate-cnv`` (`issue #1 <https://github.com/XiaoTaoWang/NeoLoopFinder/issues/1>`_). Also make sure you have run ``cooler balance`` on your cool files before ``correct-cnv`` (`issue #8 <https://github.com/XiaoTaoWang/NeoLoopFinder/issues/8>`_).


Format of the input SV list
===========================
The input SV file to the command ``assemble-complexSVs`` should contain following 6 columns separated by tab::

    chr7    chr14   ++      14000000        37500000        translocation
    chr7    chr14   --      7901149 37573191        translocation

1. **chrA**: The chromosome name of the 1st breakpoint.
2. **chrB**: The chromosome name of the 2nd breakpoint.
3. **orientation**: The orientation type of the fusion, one of ++, +-, -+, or --.
4. **b1**: The position of the 1st breakpoint on *chrA*.
5. **b2**: The position of the 2nd breakpoint on *chrB*.
6. **type**: SV type. Allowable choices are: *deletion*, *inversion*, *duplication*, and *translocation*.


Tutorial
========
This tutorial will cover the basic usage of ``assemble-complexSVs``, ``neoloop-caller`` and the
visualization module.

First, change your current working directory to the *test* folder and download the Hi-C contact map in K562::

    $ cd test
    $ wget -O K562-MboI-allReps-hg38.10K.cool https://www.dropbox.com/s/z3z5bye1tuywf18/K562-MboI-allReps-hg38.10K.cool?dl=0

To detect and assemble complex SVs in K562, submit the command below::

    $ assemble-complexSVs -O K562 -B K562-test-SVs.txt -H K562-MboI-allReps-hg38.10K.cool

The job should be finished within 1 minute, and all candidate local assemblies will be reported into
a TXT file named "K562.assemblies.txt"::

    A0	translocation,22,23290555,+,9,130731760,-	translocation,9,131280137,+,13,108009063,+	deletion,13,107848624,-,13,93371155,+	22,22300000	13,93200000
    A1	translocation,9,131280000,+,13,93252000,-	deletion,13,93371155,+,13,107848624,-	9,130720000	13,108030000
    A2	translocation,22,23290555,+,9,130731760,-	translocation,9,131280000,+,13,93252000,-	22,22300000	13,93480000
    A3	translocation,22,23290555,+,9,130731760,-	translocation,9,131199197,+,22,16819349,+	22,22300000	22,16240000
    C0	deletion,13,93371155,+,13,107848624,-	13,93200000	13,108030000
    C1	translocation,22,16819349,+,9,131199197,+	22,16240000	9,130710000
    C2	translocation,22,23290555,+,9,130731760,-	22,22300000	9,131290000
    C3	translocation,9,131280000,+,13,93252000,-	9,130720000	13,93480000
    C4	translocation,9,131280137,+,13,108009063,+	9,130720000	13,107810000

Then you can detect neo-loops on each assembly using the ``neoloop-caller`` command::

    $ neoloop-caller -O K562.neo-loops.txt -H K562-MboI-allReps-hg38.10K.cool --assembly K562.assemblies.txt --no-clustering --prob 0.95

Wait ~1 minute... The loop coordinates in both shuffled (neo-loops) and undisrupted regions near SV breakpoints will be
reported into "K562.neo-loops.txt" in `BEDPE <https://bedtools.readthedocs.io/en/latest/content/general-usage.html>`_ format::

    $ head K562.neo-loops.txt

    chr13	93270000	93280000	chr13	107860000	107870000	A0,130000,1
    chr13	93270000	93280000	chr13	107870000	107880000	A0,140000,1
    chr13	93270000	93280000	chr13	107980000	107990000	A0,250000,1
    chr13	93280000	93290000	chr13	107860000	107870000	A0,120000,1
    chr13	93280000	93290000	chr13	107870000	107880000	A0,130000,1,C0,130000,1
    chr13	93280000	93290000	chr13	107880000	107890000	A0,140000,1
    chr13	93280000	93290000	chr13	107970000	107980000	A0,230000,1
    chr13	93290000	93300000	chr13	107860000	107870000	A1,110000,1,C0,110000,1
    chr13	93290000	93300000	chr13	107870000	107880000	A1,120000,1,A0,120000,1,C0,120000,1
    chr13	93300000	93310000	chr13	107870000	107880000	C0,110000,1

The last column records the assembly IDs, the genomic distance between two loop anchors on the assembly and whether this
is a neo-loop. For example, for the 1st row above, the loop was detected on the assemblies "A0", the genomic
distance between the two anchors on this assembly is 130K (note that the distance on the reference genome is >14Mb),
and it is a neo-loop as indicated by "1".

Finally, let's reproduce the `figure 1b <https://doi.org/10.1038/s41592-021-01164-w>`_ using the python code below:

    >>> from neoloop.visualize.core import * 
    >>> import cooler
    >>> clr = cooler.Cooler('K562-MboI-allReps-hg38.10K.cool')
    >>> assembly = 'A0      translocation,22,23290555,+,9,130731760,-       translocation,9,131280137,+,13,108009063,+      deletion,13,107848624,-,13,93371155,+   22,22300000     13,93200000'
    >>> vis = Triangle(clr, assembly, n_rows=3, figsize=(7, 4.2), track_partition=[5, 0.4, 0.5], correct='sweight')
    >>> vis.matrix_plot(vmin=0)
    >>> vis.plot_chromosome_bounds(linewidth=2.5)
    >>> vis.plot_loops('K562.neo-loops.txt', face_color='none', marker_size=40, cluster=True)
    >>> vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9) 
    >>> vis.plot_chromosome_bar(name_size=11, coord_size=4.8)
    >>> vis.outfig('K562.A0.pdf')

.. image:: ./images/fig1b.png
        :align: center


Gallery
=======
In addtion to the reconstructed Hi-C maps (.cool), loops (.bedpe), and genes, the visualization module also supports plotting
RNA-Seq/ChIP-Seq/ATAC-Seq signals (.bigwig), peaks (.bed), and motifs (.bed). Below I'm going to share more examples and the
code snippets used to generate the figure.

Code Snippet 1:

    >>> from neoloop.visualize.core import * 
    >>> import cooler
    >>> clr = cooler.Cooler('SCABER-Arima-allReps.10K.cool')
    >>> List = [line.rstrip() for line in open('demo/allOnco-genes.txt')] # please find allOnco-genes.txt in the demo folder of this repository
    >>> assembly = 'A3      deletion,9,38180000,-,9,14660000,+      inversion,9,13870000,-,9,22260000,-     9,38480000      9,24220000'
    >>> vis = Triangle(clr, assembly, n_rows=5, figsize=(7, 5.2), track_partition=[5, 0.8, 0.8, 0.2, 0.5], correct='weight', span=300000, space=0.08)
    >>> vis.matrix_plot(vmin=0, cbr_fontsize=9)
    >>> vis.plot_chromosome_bounds(linewidth=2)
    >>> vis.plot_signal('RNA-Seq', 'enc_SCABER_RNASeq_rep1.bw', label_size=10, data_range_size=9, max_value=0.5, color='#E31A1C')
    >>> vis.plot_signal('H3K27ac', 'SCABER_H3K27ac_pool.bw', label_size=10, data_range_size=9, max_value=20, color='#6A3D9A')
    >>> vis.plot_genes(release=75, filter_=List, fontsize=10)
    >>> vis.plot_chromosome_bar(name_size=13, coord_size=10)
    >>> vis.outfig('SCaBER.NFIB.png', dpi=300)

Figure output 1:

.. image:: ./images/SCaBER.NFIB.png
        :align: center

Note that when you initialize a plotting object, the figure size (**figsize**), the number of tracks (**n_rows**), and the height of each
track (**track_partition**) can all be configured flexibly.

Code Snippet 2:

    >>> from neoloop.visualize.core import * 
    >>> import cooler
    >>> clr = cooler.Cooler('LNCaP-WT-Arima-allReps-filtered.mcool::resolutions/10000')
    >>> assembly = 'C26     translocation,7,14158275,+,14,37516423,+        7,13140000      14,36390000'
    >>> vis = Triangle(clr, assembly, n_rows=6, figsize=(7, 5.3), track_partition=[5, 0.4, 0.8, 0.3, 0.3, 0.5], correct='weight', span=600000, space=0.03)
    >>> vis.matrix_plot(vmin=0, cbr_fontsize=9)
    >>> vis.plot_chromosome_bounds(linewidth=2)
    >>> vis.plot_genes(filter_=['ETV1', 'DGKB', 'MIPOL1'],label_aligns={'DGKB':'right', 'ETV1':'right'}, fontsize=10) 
    >>> vis.plot_signal('DNase-Seq', 'LNCaP.DNase2.hg38.bw', label_size=10, data_range_size=9, max_value=1.8, color='#6A3D9A')
    >>> vis.plot_motif('demo/LNCaP.CTCF-motifs.hg38.txt', subset='+') # an example file LNCaP.CTCF-motifs.hg38.txt can be found at the demo folder of this repository
    >>> vis.plot_motif('demo/LNCaP.CTCF-motifs.hg38.txt', subset='-')
    >>> vis.plot_chromosome_bar(name_size=13, coord_size=10, color_by_order=['#1F78B4','#33A02C'])
    >>> vis.outfig('LNCaP.CTCF-motifs.png', dpi=300)

Figure output 2:

.. image:: ./images/LNCaP.CTCF-motifs.png
        :align: center

Code Snippet 3:

    >>> from neoloop.visualize.core import * 
    >>> import cooler
    >>> clr = cooler.Cooler('LNCaP-WT-Arima-allReps-filtered.mcool::resolutions/10000')
    >>> assembly = 'C26     translocation,7,14158275,+,14,37516423,+        7,13140000      14,36390000'
    >>> vis = Triangle(clr, assembly, n_rows=5, figsize=(7, 5.3), track_partition=[5, 0.4, 0.8, 0.8, 0.5], correct='weight', span=600000, space=0.03)
    >>> vis.matrix_plot(vmin=0, cbr_fontsize=9)
    >>> vis.plot_chromosome_bounds(linewidth=2)
    >>> vis.plot_loops('LNCaP.neoloops.txt', face_color='none', marker_size=40, cluster=True, onlyneo=True) # only show neo-loops
    >>> vis.plot_genes(filter_=['ETV1', 'DGKB', 'MIPOL1'],label_aligns={'DGKB':'right', 'ETV1':'right'}, fontsize=10)
    >>> vis.plot_signal('DNase-Seq', 'LNCaP.DNase2.hg38.bw', label_size=10, data_range_size=9, max_value=1.8, color='#6A3D9A')
    >>> vis.plot_arcs(lw=1.5, cutoff='top', gene_filter=['ETV1'], arc_color='#666666') # ETV1-related neo-loops
    >>> vis.plot_chromosome_bar(name_size=13, coord_size=10, color_by_order=['#1F78B4','#33A02C'])
    >>> vis.outfig('LNCaP.arcs.png', dpi=300)

Figure output 3:

.. image:: ./images/LNCaP.arcs.png
        :align: center

Note that both **plot_loops** and **plot_genes** need to be called before **plot_arcs**.