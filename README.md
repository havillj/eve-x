# EVE

## Dependencies

### Python modules

EVE is written in Python 3 (version 3.6 or later) and requires the following nonstandard modules:
- biopython
- dna_features_viewer
- intervaltree
- lxml
- matplotlib
- numpy
- pysam

They should be easy to install using `pip` (or `pip3`).  For example:

```
pip3 install biopython
```

Or, if you don't have root access to your system:

```
pip3 install --user biopython
```

### Samtools

Samtools can be found [here](http://www.htslib.org).

Once installed, assign the path to `samtools` to the variable `SAMTOOLS_EXEC` in [config.py](python/config.py).  For example:

```
SAMTOOLS_EXEC = '/usr/bin/samtools'
```

### SPAdes

The SPAdes assembler is available [here](https://cab.spbu.ru/software/spades/).

Once installed, assign the path to `spades.py` to the variable `SPADES_EXEC` in `config.py`.  For example:

```
SPADES_EXEC = '/usr/local/bin/spades.py'
```

### BLAST

BLAST can be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

Once installed, assign the path to the `blastn` executable to the variable `BLAST_EXEC` in `config.py`.  For example:

```
BLAST_EXEC = '/usr/local/bin/blastn'
```

## Set up

### Directories

The following four directories should be specified in `config.py`:

```
ROOT_DIR = '/Volumes/Data/'              # root directory where all data files related to EVE will be located
SPECIMENS_DIR = ROOT_DIR + 'specimens/'  # path to original specimen BAM files
RESULTS_DIR = ROOT_DIR + 'results/'      # path where all results should be written
BLASTDB_DIR = ROOT_DIR + 'blastdb/'      # path to blast databases
```

If you wish to test your installation, a sample *Aedes aegypti* specimen BAM file can be found [here](https://www.dropbox.com/s/0t7ovuujo0iiroi/La_Lope-Gabon-10.LIN210A1646.sorted.deduped.merged.bam?dl=0).

See the [Directory structure](#directory-structure) section below for details on the contents of the `RESULTS_DIR` directory.

### BLAST databases

You will need two BLAST nucleotide databases: 
- one containing the virus sequences for which you wish to search and
- one containing the host reference genome

Assuming your virus sequences are contained in a FASTA file named `virus.fasta`, you can create the database with

```
makeblastdb -in virus.fasta -dbtype nucl -parse_seqids -title "virusdb" -out virusdb
```

Similarly, if your host reference genome is contained in `host.fasta`, you can create the host database with

```
makeblastdb -in host.fasta -dbtype nucl -parse_seqids -title "hostdb" -out hostdb
```

Once these databases exist, assign their names to the variables `VIRUS_DB` and `HOST_DB` in `config.py`:

```
VIRUS_DB = 'virusdb'
HOST_DB = 'aegyptidb'
```

To make diagrams more readable, also edit the variable `CHR_NAMES` in `config.py` so that the accession number of each host chromosome corresponds to the correct chromosome number.  For example, for *Aedes aegypti*:

```
CHR_NAMES = {'NC_035107.1': 'Chr1', 'NC_035108.1': 'Chr2', 'NC_035109.1': 'Chr3'}
```

### Host genome features

To identify repeat sequences and other genomic features near putative viral insertions, EVE needs files in GFF3 format that list these features.  For example, these files for *Aedes aegpyti* can be obtained from http://vectorbase.org ([here](https://vectorbase.org/common/downloads/release-52/AaegyptiLVP_AGWG/gff/data/VectorBase-52_AaegyptiLVP_AGWG.gff) and [here](https://vectorbase.org/common/downloads/Legacy%20VectorBase%20Files/Aedes-aegypti/Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3.gz)).  Once downloaded and decompressed, add their names to the variables `BASE_FEATURES_GFF3_FILENAME` and `REPEAT_FEATURES_GFF3_FILENAME` in `config.py`.  For example:

```
GFF3_DIR = ROOT_DIR + 'gff3/'
BASE_FEATURES_GFF3_FILENAME = GFF3_DIR + 'Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3'
REPEAT_FEATURES_GFF3_FILENAME = GFF3_DIR + 'Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3'
```

The GFF3 files may identify host chromosomes simply by number.  To translate between these and the chromosomes' accession numbers in the host reference genome, edit the variable `CHR_NUMBERS` in `config.py`.  For example, for *Aedes aegypti*, this would look like:

```
CHR_NUMBERS = {'1': 'NC_035107.1', '2': 'NC_035108.1', '3': 'NC_035109.1', 'MT': 'NC_035159.1'}
```
### Specimen populations

If your specimens come from different populations, you can specify those by editing the `POPULATIONS` dictionary in `config.py`.  The keys are the names of the populations and the values are patterns that appear in the file names of specimens from those populations.  For example, in the following, the region `'USA'` is assigned to any specimen containing either `'USA'` or `'AZ'` in its file name.

```
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
```
### Preferred virus accession numbers

Since there are commonly many genomic sequences in a viral database corresponding to the same viral species, EVE only retains hits from one representative of a species in each contig.  You can optionally specify a preferred accession number for each species by editing the `PREFERRED_ACCS` variable in `config.py`.  For example:

```
PREFERRED_ACCS = {'Aedes anphevirus':               ['gb|MH037149.1|'], 
                  'Australian Anopheles totivirus': ['ref|NC_035674.1|'], 
                  'Culex ohlsrhavirus':             ['ref|NC_035132.1|'], 
                  'Phasi Charoen-like phasivirus':  ['ref|NC_038263.1|'], 
                  'Wuhan Mosquito Virus 6':         ['gb|MF176251.1|']}
```

## Parameters

There are several ways to customize the execution of EVE.  These options are described in detail in the comments of [config.py](python/config.py).

## Usage

To run EVE on the specimen files in the directory `SPECIMENS_DIR`, simply execute
```
python3 eve.py
```

Any of the options in the dictionary `config` in `config.py` may be specified on the command line in the form
```
python3 eve.py --OPTION1=VALUE1 --OPTION2=VALUE2 ...
```

For example:
```
python3 eve.py --DO_CLUSTERING=False --MIN_HITS_TO_SHOW_VIRUS=4
```

EVE will generate many files following the directory structure in the [next section](#directory-structure).

If you wish, you can manually adjust any of the clustered aligned sequence files and then run

```
python3 insertsites.py /path/to/clustered_sequences.fasta
```

to locate putative insertion sites, generate FASTA files for each cluster, and get the geographical distribution of the clusters.  If you wish to also identify genomic features near the insertion sites, specify `--find_features`, as below:

```
python3 insertsites.py --find_features /path/to/clustered_sequences.fasta
```

## Directory structure

In the following directory tree, `N` is used to represent the last in a series, not any particular number.

```
<ROOT_DIR>
├── specimens
│   ├── <SPECIMEN 1>.bam
│   :
│   └── <SPECIMEN N>.bam
│   
└── results
    ├── results.tsv
    ├── specimens
    │   ├── <SPECIMEN 1>
    │   │   ├── <SPECIMEN 1>_unmapped_with_mates.bam
    │   │   ├── diagrams
    │   │   │   ├── <VIRUS_FAMILY 1>
    │   │   │   │   ├── <CONTIG 1>.pdf
    │   │   │   │   :
    │   │   │   │   └── <CONTIG N>.pdf
    │   │   │   :
    │   │   │   └── <VIRUS_FAMILY N>
    │   │   │       ├── 
    │   │   │       :   (as above)
    │   │   │       └── 
    │   │   ├── scaffolds
    │   │   │   ├── blast_scaffolds.csv
    │   │   │   └── scaffolds.fasta
    │   │   ├── sequences
    │   │   │   ├── <SPECIMEN 1>_hits_aligned.fasta
    │   │   │   └── <SPECIMEN 1>_hits_unaligned.fasta
    │   │   ├── xml
    │   │   │   └── <SPECIMEN 1>_hits.xml
    │   │   └── spades
    │   │       ├──
    │   │       :   (SPAdes work files)
    │   │       └── 
    │   :
    │   └── <SPECIMEN N>
    │       ├── 
    │       :   (as above)
    │       └──
    │
    └── viruses
        ├── diagrams
        │   ├── <VIRUS_FAMILY 1>
        │   │   ├── <VIRUS_FAMILY 1 VIRUS 1>_all.pdf
        │   │   :
        │   │   └── <VIRUS_FAMILY 1 VIRUS N>_all.pdf
        │   :
        │   └── <VIRUS_FAMILY N>
        │       ├── 
        │       :   (as above)
        │       └── 
        └── sequences
            ├── <VIRUS_FAMILY 1>
            │   ├── <VIRUS_FAMILY 1>_per_contig_aligned.fasta
            │   ├── <VIRUS_FAMILY 1>_per_contig_unaligned.fasta
            │   ├── <VIRUS_FAMILY 1>_per_specimen_aligned.fasta
            │   ├── <VIRUS_FAMILY 1 VIRUS 1>_per_contig_aligned.fasta
            │   ├── <VIRUS_FAMILY 1 VIRUS 1>_per_contig_unaligned.fasta
            │   ├── <VIRUS_FAMILY 1 VIRUS 1>_per_contig_with_flanks.fasta
            │   ├── <VIRUS_FAMILY 1 VIRUS 1>_per_specimen_aligned.fasta
            │   :
            │   ├── <VIRUS_FAMILY 1 VIRUS N>_per_contig_aligned.fasta
            │   ├── <VIRUS_FAMILY 1 VIRUS N>_per_contig_unaligned.fasta
            │   ├── <VIRUS_FAMILY 1 VIRUS N>_per_contig_with_flanks.fasta
            │   ├── <VIRUS_FAMILY 1 VIRUS N>_per_specimen_aligned.fasta
            │   └── clustered
            │       ├── <VIRUS_FAMILY 1 VIRUS 1>_per_contig_aligned_clustered_k?.fasta
            │       ├── <VIRUS_FAMILY 1 VIRUS 1>_per_contig_aligned_clustered_k?.txt
            │       ├── <VIRUS_FAMILY 1 VIRUS 1>
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_C1.fasta
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_C1_with_flanks.fasta
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_C2.fasta
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_C2_with_flanks.fasta
            │       │   :
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_CN.fasta
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_CN_with_flanks.fasta
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_insertpositions_Chr1.tsv
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_insertpositions_Chr2.tsv
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_insertpositions_Chr3.tsv
            │       │   ├── <VIRUS_FAMILY 1 VIRUS 1>_insertpositions.pdf
            │       │   └── <VIRUS_FAMILY 1 VIRUS 1>_insertpositions.txt
            │       :
            │       ├── <VIRUS_FAMILY 1 VIRUS N>_per_contig_aligned_clustered_k?.fasta
            │       ├── <VIRUS_FAMILY 1 VIRUS N>_per_contig_aligned_clustered_k?.txt
            │       └── <VIRUS_FAMILY 1 VIRUS N>
            │           ├── 
            │           :   (as above)
            │           └── 
            :
            └── <VIRUS_FAMILY N>
                ├── 
                :   (as above)
                └── 
```

## Output file descriptions

For each specimen, two FASTA files are created:
- `<SPECIMEN>_hits_unaligned.fasta` contains the sequences of all viral hits found for that specimen.

- `<SPECIMEN>_hits_aligned.fasta` contains the sequences of all viral hits, aligned to the reference genome of the respective virus.  The FASTA file alternates between specimen and viral reference sequences, each pair representing an alignment. 

For each virus species, four additional FASTA files are created:
- `<VIRUS ACC>_per_specimen_aligned.fasta` contains one entry for each specimen in which this virus was found.  Each entry contains all regions of the virus found in that specimen, aligned to the virus reference genome, which is the first entry in the file.

- `<VIRUS ACC>_per_contig_aligned.fasta` contains one entry per contig per specimen in which the virus was found.  These are also aligned to the virus reference genome and used later to cluster hits.

- `<VIRUS ACC>_per_contig_unaligned.fasta` also contains one entry per contig per specimen, but they are not aligned.  Each of these entries represents a putative NIRVS.

- `<VIRUS ACC>_per_contig_with_flanks.fasta` contains the sequences in the previous file, with flanking regions from the contig in which it was found

Third, for each virus family, three FASTA files are created:
- `<VIRUS FAMILY>_per_specimen_aligned.fasta`
- `<VIRUS FAMILY>_per_contig_aligned.fasta`
- `<VIRUS FAMILY>_per_contig_unaligned.fasta`

In the first two files, sequences are aligned to the reference genome of the virus species with the most hits within a viral family.  The third file contains the raw, unaligned sequences.

Once the hits are k-means clustered, these results are written to new FASTA files with filename ending `_clustered_kx.fasta` (where `x` is the final value of k) and with each sequence identifier prefixed by `Cy_`, where `y` is a cluster number. 

The raw putative insertion sites are written to tab-delimited (TSV) files, one per chromosome, named `<VIRUS ACC>_insertpositions_chrX.tsv` where `X` is a chromosome number.

Putative insertion positions are consolidated by groups of specimens sharing the same set of insertion sites and written to a text file named `<VIRUS ACC>_insertpositions.txt`.
