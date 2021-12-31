""" 
Configurable path and file names
--------------------------------
"""

SPADES_EXEC = '/Volumes/Data/bin/SPAdes-3.14.1-Linux/bin/spades.py'

ROOT_DIR = '/Volumes/Data2/'

SPECIMENS_DIR = ROOT_DIR + 'specimens/' # path for specimens directories
SPECIMEN_RESULTS_DIR = 'results'        # results subdirectory for each specimen

RESULTS_DIR = ROOT_DIR + 'results/'     # path for consolidated results
VIRUSES_DIR = RESULTS_DIR + 'viruses/'  # path for per-virus results
SEQUENCES_DIR = VIRUSES_DIR + 'sequences/'  # path for per-virus FASTA files
DIAGRAMS_DIR = VIRUSES_DIR + 'diagrams/'    # path for per-virus diagrams

GB_DIR = ROOT_DIR + 'gb/'
FASTA_DIR = ROOT_DIR + 'fasta/'
GFF3_DIR = ROOT_DIR + 'gff3/'
FEATURE_TREES_PICKLED_FILENAME = ROOT_DIR + 'featuretrees.pickle'

FAMILY_CSV = ROOT_DIR + 'families.csv'

LOGFILE_DIR = ROOT_DIR + 'log/'
LOGFILE_NAME = 'eve_' + time.strftime('%Y-%m-%d_%H:%M:%S') + '.log'
LOGFILE_PATH = Path(LOGFILE_DIR) / LOGFILE_NAME 

if 'logfile' not in globals():
    if not Path(LOGFILE_DIR).exists():
        os.system('mkdir ' + LOGFILE_DIR)
    logFile = open(LOGFILE_PATH , 'w')

""" 
Configurable variables by function
----------------------------------

    *** blastScaffolds ***
    
    EVALUE_VIRUS (floating point number)
       This is the minimum e-value required for a viral hit in a contig.
       
    *** getHits ***
    
    MIN_PIDENT (integer between 0 and 100)
       If the maximum percent identity among all viral hits in a contig is >=  
       this value, then the viral hit with the maximum percent identity is
       deemed to be the "best hit" in the contig.  Otherwise, the longest viral 
       hit in the contig is the "best hit."

    EVALUE_HOST (floating point number)
       This is the minimum e-value required for a host hit in a contig.
       
    MAX_FLANK_DISTANCE (positive integer)
       This is the maximum distance between a pair of host hits in a contig for
       them to be considered a matching pair.

    ALLOWED_OVERLAP (non-negative integer)
       Given a pair of host hits that are within MAX_FLANK_DISTANCE bp of 
       each other in a contig, this is the number of bp that the end of the 
       first/left hit may overlap the beginning of the second/right hit to
       be considered a matching pair.

    OVERLAP_FRACTION (floating point number between 0.0 and 1.0)
       If at least this fraction of a (mosaic) viral hit is overlapped by host  
       hits in a contig, then the viral hit is considered to be present in the 
       reference genome.

    HOST_OVERLAP_FRACTION (floating point number between 0.0 and 1.0)
       If a host hit overlaps at least this fraction of the span of a viral hit,
       then the host hit is considered to overlap the viral hit.
       
    COMPLEXITY_K (positive integer)
       This is the value of k used to compute the complexity of a sequence.
       Complexity is defined to be U_1 * U_2 * ... * U_K, where each U_k
       is the number of unique k-mers that appear in the sequence  divided
       by the maximum possible number in a sequence of that size.
    
    COMPLEXITY_CUTOFF (floating point number between 0.0 and 1.0)
       If a viral hit has complexity less than this value, it is discarded.
    
    MIN_VIRUS_HIT_LENGTH (positive integer)
       A viral hit with length less than this value is discarded.
    
    MAX_HOST_VIRUS_DISTANCE (non-negative integer)
    MIN_FLANKING_HIT_LENGTH (positive integer)
       In a contig, if a host hit is farther than MAX_HOST_VIRUS_DISTANCE bp 
       from the left or right end of the span of the viral hit and its length 
       is less than MIN_FLANKING_HIT_LENGTH, then the host hit is discarded.
       
       MAX_HOST_VIRUS_DISTANCE is also used in getInsertSites to determine if a
       host hit is close enough to a viral hit in their contig to warrant a
       search for annotated features in the host genome.
      
    *** drawXML ***
    
    FLANK_DRAW_LIMIT (positive integer)
       This is the maximum number of host hits to draw in the contig diagrams.
       Pairs of matching flanking regions are applied to the limit first, then
       pairs of inverted flanking regions, then other host hits.
    
    BEST_HITS_ONLY (boolean)
       If True, only draw the contig diagram for the viral hit determined to be 
       the "best hit".  Also used in consolidateAll.
       
    *** consolidateAll ***
    
    OMIT_EVES_IN_REFERENCE (boolean)
       If True, EVEs determined to be present in the host reference genome will
       be omitted from the summary table and diagrams.
    
    BEST_HITS_ONLY (boolean)
       If True, only the viral hit determined to be the "best hit" in each
       contig is included in the summary table and diagrams.
       
    MIN_HITS_TO_SHOW_VIRUS (positive integer)
       If fewer than this number of specimens are found to have a particular
       viral hit, then that virus is not shown in the summary table and diagrams.
    
    DO_CLUSTERING (boolean)
       If True, cluster the aligned viral sequences.
       
    *** kmeans ***
    
    MIN_K, MAX_K (positive integers)
       These are the minimum and maximum values of k to be used in k-means
       clustering.
    
    KMEANS_TRIALS (positive integer)
       For each value of k, this is the number of times the sequences are 
       clustered to find the best set of k clusters.
    
    MAX_KMEANS_ITERATIONS (positive integer)
       This is the maximum number of iterations performed by the k-means 
       clustering algorithm.  If the set of clusters stops changing earlier,
       this number is not reached.
       
    *** getInsertSites ***
    
    FEATURE_SEARCH_DIST (non-negative integer)
       This is the maximum distance in the host genome from a host hit to 
       search for annotated features in the host genome.
       
     MIN_HITS_TO_SHOW_POSITION
       In a cluster of viral hits, this is the minimum number of host hits that 
       must be found at a particular host location to show that location in the 
       insertpositions text file.
       
"""

config = {'MIN_PIDENT': 80,
          'ALLOWED_OVERLAP': 30,
          'BEST_HITS_ONLY': True,
          'COMPLEXITY_CUTOFF': 0.75,
          'COMPLEXITY_K': 3,
          'DO_CLUSTERING': False,
          'EVALUE_HOST': 1e-10,
          'EVALUE_VIRUS': 1e-10,
          'FEATURE_SEARCH_DIST': 10000,
          'FLANK_DRAW_LIMIT': 5,
          'HOST_OVERLAP_FRACTION': 0.8,
          'KMEANS_TRIALS': 15,
          'MAX_FLANK_DISTANCE': 20000,
          'MAX_HOST_VIRUS_DISTANCE': 200,    
          'MAX_K': 12,
          'MAX_KMEANS_ITERATIONS': 100,
          'MIN_FLANKING_HIT_LENGTH': 100,
          'MIN_HITS_TO_SHOW_POSITION': 1,
          'MIN_HITS_TO_SHOW_VIRUS': 2,
          'MIN_K': 2,
          'MIN_VIRUS_HIT_LENGTH': 100,
          'OMIT_EVES_IN_REFERENCE': True,
          'OVERLAP_FRACTION': 0.95}

"""
Specimen parameters
"""

HOST_DB = 'aegyptidb'
VIRUS_DB = 'virusdb'
BASE_FEATURES_GFF3_FILENAME = GFF3_DIR + 'Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3'
REPEAT_FEATURES_GFF3_FILENAME = GFF3_DIR + 'Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3'
CHR_NAMES = {'NC_035107.1': 'Chr1', 'NC_035108.1': 'Chr2', 'NC_035109.1': 'Chr3'}   # for diagrams
CHR_NUMBERS = {'1': 'NC_035107.1', '2': 'NC_035108.1', '3': 'NC_035109.1', 'MT': 'NC_035159.1'} # for gff3 files

""" 
A dictionary used in the getSpecimenLabel function to determine the 
population to which a particular specimen belongs.  If any of the strings 
in the list corresponding to a particular key is contained in the specimen 
name, then that key is used as the population name. 
"""

POPULATIONS = {'Angola': ['Angola'], 
               'Argentina': ['Argentina', 'US_U'], 
               'Australia': ['Australia'], 
               'Brazil': ['Brazil'], 
               'French-Polynesia': ['FrenchPolynesia'],
               'Gabon': ['Gabon'], 
               'Mexico': ['Mexico'], 
               'Philippines': ['Philippines'], 
               'South-Africa': ['South_Africa'], 
               'Thailand': ['Thailand'], 
               'USA': ['USA', 'AZ'], 
               'Vietnam': ['Vietnam']}
               
"""
This dictionary is used to choose a preferred accession number if there are
hits to multiple database sequences belonging to the same species.  For any
species not in this dictionary, the accession number with the longest hit is
chosen.
"""
               
PREFERRED_ACCS = {'Aedes anphevirus':               ['gb|MH037149.1|'], 
                  'Australian Anopheles totivirus': ['ref|NC_035674.1|'], 
                  'Culex ohlsrhavirus':             ['ref|NC_035132.1|'], 
                  'Phasi Charoen-like phasivirus':  ['ref|NC_038263.1|'], 
                  'Wuhan Mosquito Virus 6':         ['gb|MF176251.1|']}
                  
#'MIN_FLANKING_QSTART': 10        # min distance a flanking hit must be from either end of the contig  # REMOVE PROBABLY
#'MAX_FLANKING_HITS': 10
#'FIND_FEATURES': True
