from intervaltree import Interval, IntervalTree
import collections, logging, sys

LOG_LEVEL = logging.WARNING

# coolbox.utilities.logtools
def get_logger(name, file_=sys.stderr, level=LOG_LEVEL):
    FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
    formatter = logging.Formatter(fmt=FORMAT)
    if isinstance(file_, str):
        handler = logging.FileHandler(file_)
    else:
        handler = logging.StreamHandler(file_)
    handler.setFormatter(formatter)
    log = logging.getLogger(name)
    log.addHandler(handler)
    log.setLevel(level)
    return log

log = get_logger(__name__)

# coolbox.utilities.filetool
def opener(filename):
    """
    Determines if a file is compressed or not

    >>> import gzip
    >>> msg = "hello blablabla"
    >>> tmp_f_raw  = open("/tmp/test_opener.txt", 'w')
    >>> tmp_f_raw.write(msg)
    15
    >>> tmp_f_raw.close()
    >>> tmp_f_gzip = gzip.open('/tmp/test_opener.txt.gz', 'wb')
    >>> tmp_f_gzip.write(to_bytes(msg))
    15
    >>> tmp_f_gzip.close()

    >>> test_raw = opener(tmp_f_raw.name)
    >>> type(test_raw)
    <class '_io.BufferedReader'>
    >>> test_gzip = opener(tmp_f_gzip.name)
    >>> type(test_gzip)
    <class 'gzip.GzipFile'>

    >>> test_raw.close()
    >>> test_gzip.close()

    >>> import os
    >>> os.remove(test_raw.name)
    >>> os.remove(test_gzip.name)

    """
    import gzip
    f = open(filename, 'rb')
    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f

def to_string(s):
    """
    Convert bytes, bytes list to string, string list.

    >>> to_string("hello")
    'hello'
    >>> to_string(b"hello")
    'hello'
    >>> to_string([b'hello', b'world'])
    ['hello', 'world']
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        return s.decode('ascii')
    if isinstance(s, list):
        return [to_string(x) for x in s]
    return s

# coolbox.utilities.genome
class GenomeRange(object):
    """
    Express a range on the genome.

    Attributes
    ----------
    chrom : str
        chromosome

    start : int
        start position

    end : int
        end position

    """

    def __init__(self, *args):
        """
        >>> range1 = GenomeRange("chr1", 1000, 2000)
        >>> str(range1)
        'chr1:1000-2000'
        >>> range2 = GenomeRange("chr2:2000-4000")
        >>> (range2.chrom, range2.start, range2.end)
        ('chr2', 2000, 4000)
        >>> range3 = GenomeRange("chr1", 2000, 1000)
        Traceback (most recent call last):
        ...
        ValueError: Please check that the region end is larger than the region start. Values given: start: 2000, end: 1000
        """
        if len(args) == 1:
            chrom, start, end = GenomeRange.parse_region_string(args[0])
        elif len(args) == 3:
            chrom, start, end = args
        else:
            raise ValueError("inappropriate init arguments. "
                             "correct example: `range1 = GenomeRange(\"chr1:1000-2000\")` or "
                             "`range1 = GenomeRange(\"chr1\", 1000, 2000)`")

        if end < start:
            raise ValueError("Please check that the region end is larger than the region start. "
                             "Values given: start: {}, end: {}".format(start, end))

        self.chrom = chrom
        self.start = start
        self.end = end

    @staticmethod
    def parse_region_string(region_string):
        """
        splits a region string into
        a (chrom, start, end) tuple

        Parameters
        ----------
        region_string : str
            Region string to be parsed, like: "chr:start-end"

        Return
        ------
        result : tuple of str
            Result tuple (chrom, start, end)

        >>> GenomeRange.parse_region_string("chr1:10-20")
        ('chr1', 10, 20)
        >>> GenomeRange.parse_region_string("chr1:0")
        Traceback (innermost last):
         ...
        ValueError: Failure to parse region string, please check that region format should be like "chr:start-end".
        """
        if region_string:
            # separate the chromosome name and the location using the ':' character
            chrom, position = region_string.strip().split(":")

            # clean up the position
            for char in ",.;|!{}()":
                position = position.replace(char, '')

            position_list = position.split("-")
            try:
                region_start = int(position_list[0])
                region_end = int(position_list[1])
            except:
                raise ValueError("Failure to parse region string, please check that region format "
                                 "should be like \"chr:start-end\".")

            return chrom, region_start, region_end

    def change_chrom_names(self):
        """
        >>> range1 = GenomeRange("chr1", 1000, 2000)
        >>> range1.chrom
        'chr1'
        >>> range1.change_chrom_names()
        >>> range1.chrom
        '1'
        >>> range1.change_chrom_names()
        >>> range1.chrom
        'chr1'
        """
        self.chrom = change_chrom_names(self.chrom)

    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)

    @property
    def length(self):
        """
        >>> range1 = GenomeRange("chr1", 0, 1000)
        >>> range1.length
        1000
        """
        return self.end - self.start

    def __eq__(self, other):
        """
        >>> GenomeRange('chr1', 1000, 2000) == GenomeRange("chr1:1000-2000")
        True
        >>> GenomeRange("chr1:1000-2000") == GenomeRange("1:1000-2000")
        False
        """
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    def __contains__(self, another):
        if another.chrom != self.chrom:
            return False
        if another.start < self.start:
            return False
        if another.end > self.end:
            return False
        return True


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.

    >>> change_chrom_names("chr1")
    '1'
    >>> change_chrom_names("1")
    'chr1'
    """
    # TODO: mapping from chromosome names like mithocondria is missing
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom

# coolbox.utilities.bed
class ReadBed(object):
    """
    Reads a bed file. Based on the number of fields
    it tries to guess the type of bed file used. Current options
    are bed3, bed6 and bed12

    Example:
    bed = readBed(open("file.bed", 'r'))
    for interval in bed:
        print(interval['start'])

    """

    def __init__(self, file_handle):

        self.file_type = None
        self.file_handle = file_handle
        self.line_number = 0
        # guess file type
        fields = self.get_no_comment_line()
        fields = to_string(fields)
        fields = fields.split()

        self.guess_file_type(fields)
        self.file_handle.seek(0)
        self.prev_chrom = None
        self.prev_start = -1
        self.prev_line = None

        # list of bed fields
        self.fields = ['chromosome', 'start', 'end',
                       'name', 'score', 'strand',
                       'thick_start', 'thick_end',
                       'rgb', 'block_count',
                       'block_sizes', 'block_starts']

        if self.file_type == 'bed12':
            self.BedInterval = collections.namedtuple('BedInterval', self.fields)
        elif self.file_type == 'bed9':
            self.BedInterval = collections.namedtuple('BedInterval', self.fields[:9])
        else:
            self.BedInterval = collections.namedtuple('BedInterval', self.fields[:6])

    def __iter__(self):
        return self

    def get_no_comment_line(self):
        """
        Skips comment lines starting with '#'
        "track" or "browser" in the bed files
        """
        line = next(self.file_handle)
        line = to_string(line)
        if line.startswith("#") or line.startswith("track") or \
                line.startswith("browser") or line.strip() == '':
            line = self.get_no_comment_line()

        self.line_number += 1
        return line

    def guess_file_type(self, line_values):
        """try to guess type of bed file by counting the fields
        """
        if len(line_values) == 3:
            self.file_type = 'bed3'
        elif len(line_values) == 4:
            self.file_type = 'bedgraph'
        elif len(line_values) == 6:
            self.file_type = 'bed6'
        elif len(line_values) == 12:
            self.file_type = 'bed12'
        elif len(line_values) == 9:
            # this is a case where a specific color is encoded in the 10 field of the bed file
            self.file_type = 'bed9'
        elif len(line_values) > 6:
            # assume bed6
            self.file_type = 'bed6'
            log.warning("Number of fields in BED file is not standard. Assuming bed6\n")
        else:
            # assume bed3
            self.file_type = 'bed3'
            log.warning("Number of fields in BED file is not standard. Assuming bed3\n")
        return self.file_type

    def next(self):
        """
        Return
        ------
        bedInterval object
        """
        line = self.get_no_comment_line()

        bed = self.get_bed_interval(line)
        if self.prev_chrom == bed.chromosome:
            assert self.prev_start <= bed.start, \
                "Bed file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self.file_handle.name, self.prev_line, line)

        self.prev_chrom = bed.chromosome
        self.prev_start = bed.start
        self.prev_line = line

        return bed

    def __next__(self):
        """
        Return
        ------
        bedInterval object
        """
        line = self.get_no_comment_line()

        bed = self.get_bed_interval(line)
        if self.prev_chrom == bed.chromosome:
            assert self.prev_start <= bed.start, \
                "Bed file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self.file_handle.name, self.prev_line, line)

        self.prev_chrom = bed.chromosome
        self.prev_start = bed.start
        self.prev_line = line

        return bed

    def get_bed_interval(self, bed_line):
        r"""
        Processes each bed line from a bed file, casts the values and returns
        a namedtuple object

        >>> bed_line="chr1\t0\t1000\tgene_1\t0.5\t-\t0\t1000\t0\t3\t10,20,100\t20,200,700"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed = bed_f.get_bed_interval(bed_line)
        >>> bed.chromosome
        'chr1'
        >>> bed.block_starts
        [20, 200, 700]

        >>> bed_line="chr2\t0\t1000\tgene_1\t0.5\t-\n"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed_f.get_bed_interval(bed_line)
        BedInterval(chromosome='chr2', start=0, end=1000, name='gene_1', score=0.5, strand='-')
        """

        line_data = bed_line.strip()
        line_data = to_string(line_data)
        line_data = line_data.split("\t")

        if self.file_type == 'bed12':
            assert len(line_data) == 12, "File type detected is bed12 but line {}: {} does " \
                                         "not have 12 fields.".format(self.line_number, bed_line)

        elif self.file_type == 'bed3':
            assert len(line_data) == 3, "File type detected is bed3 but line {}: {} does " \
                                        "not have 3 fields.".format(self.line_number, bed_line)

        elif self.file_type == 'bed6':
            assert len(line_data) == 6, "File type detected is bed6 but line {}: {} does " \
                                        "not have 6 fields.".format(self.line_number, bed_line)
        line_values = []
        for idx, r in enumerate(line_data):
            # first field is always chromosome/contig name
            # and should be cast as a string
            # same for field 3 (name)
            if idx in [0, 3]:
                line_values.append(r)
            # check field strand
            elif idx == 5:
                if r not in ['+', '-', '.']:
                    if r == '1':
                        r = '+'
                    elif r == '-1':
                        r = '-'
                    else:
                        log.warning("*Warning, invalid strand value found {} for line #{}:\n{}\n "
                                    "Setting strand to '.'\n".format(r, bed_line, self.line_number))
                        r = '.'
                line_values.append(r)

            elif idx in [1, 2, 6, 7, 9]:
                # start and end fields must be integers, same for thichStart(6),
                # and thickEnd(7) and blockCount(9) fields
                try:
                    line_values.append(int(r))
                except ValueError:
                    log.warning("Value: {} in field {} at line {} is not an integer\n".format(r, idx + 1,
                                                                                              self.line_number))
                    return dict()
            # check item rgb
            elif idx == 8:
                r = to_string(r)
                rgb = r.split(",")
                if len(rgb) == 3:
                    try:
                        r = list(map(int, rgb))
                    except ValueError as detail:
                        log.warning("Error reading line: #{}. The rgb field {} is not "
                                    "valid.\nError message: {}\n".format(self.line_number, r, detail))
                line_values.append(r)

            elif idx in [10, 11]:
                # this are the block sizes and block start positions
                r = to_string(r)
                r_parts = r.split(',')
                try:
                    r = [int(x) for x in r_parts if x != '']
                except ValueError as detail:
                    log.warning("Error reading line #{}. The block field {} is not "
                                "valid.\nError message: {}\n".format(self.line_number, r, detail))
                line_values.append(r)

            else:
                try:
                    tmp = float(r)
                except ValueError:
                    tmp = r
                except TypeError:
                    tmp = r
                line_values.append(tmp)

        assert line_values[2] > line_values[1], \
            "Start position larger or equal than end for line #{}:\n{}\n".format(self.line_number,
                                                                                 bed_line)

        if self.file_type == 'bed3':
            line_values = line_values[0:3]
            # in case of a bed3, the id, score and strand
            # values are added as ".", 0, "." respectively
            line_values.extend([".", 0, "."])
        elif self.file_type == 'bed6':
            line_values = line_values[0:6]

        return self.BedInterval._make(line_values)

# coolbox.plots.track.base
class TrackPlot(object):
    """
    The TrackPlot object is a holder for all tracks that are to be plotted.
    For example, to plot a bedgraph file a new class that extends TrackPlot
    should be created.

    It is expected that all TrackPlot objects have a plot method.

    """

    def __init__(self, *args, **kwargs):
        if not hasattr(self, 'properties'):
            self.properties = args[0]
        super().__init__()

    def plot_label(self):
        if hasattr(self, 'label_ax'):
            self.label_ax.text(0.15, 0.5, self.properties['title'],
                               horizontalalignment='left', size='large',
                               verticalalignment='center')
