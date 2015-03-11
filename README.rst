===============================================================================
TREAT - Trypanosome RNA Editing Alignment Tool
===============================================================================

------------------------------------------------------------------------
About
------------------------------------------------------------------------

TREAT is a multiple sequence alignment and visualization tool specifically
designed for sequences containing uridine insertion/deletion RNA editing. This
phenomenon occurs in trypanosomes, a group of unicellular parasitic flagellate
protozoa such as Trypanosoma brucei which causes African sleeping sickness. The
pre-mRNA sequences in trypanosomes are posttranscriptionally edited by the
insertion/deletion of uridylate residues. TREAT aligns sequences (typically
from RNA-Seq experiments) using 3 bases (A,C,G) and assembles editing sites to
detect the extent of uridine (T) editing. TREAT is written in Go and is
released under the GPLv3.0 free software license. 

------------------------------------------------------------------------
Installation
------------------------------------------------------------------------

Fetch from github::

  $ go get install github.com/ubccr/treat/...

Install from source::

  $ git clone https://github.com/ubccr/treat.git treat
  $ cd treat
  $ go build ./...

------------------------------------------------------------------------
Alignment Example
------------------------------------------------------------------------

TREAT alignments require 2 files in FASTA format as input. The fragment of DNA
(read, clone, etc.) to align, for example from an RNA-Seq experiment, and 2
template sequences: Fully Edited and Pre-Edited. The Fully Edited template
represents a mature edited mRNA (completely precisely edited mRNA). The
Pre-Edited template represents the sequence that will be edited in the mature
RNA. The template input file should include the Fully Edited
template sequence first, followed by the Pre-Edited template, and optionally any
alternative edited sequences. For example::

  # simple-templates.fasta
  >Fully Edited
  CTTAATACACTTTTGATTAACAAACTTTAAA
  >Pre-Edited
  CTAATTACACTTTGATAACAAACTAAA

  # simple-sequences.fasta
  >example-1
  CTTAATTACACTTTGATTAACAAACTTTAAA

We save the above sequence files and run the alignment using TREAT::

  $ treat align -t simple-templates.fasta -f simple-sequences.fasta -v
  ================================================================================
  example-1
  ================================================================================
  Editing stops at site: 11
  Junction ends at site: 18
  Junction length: 7
  ================================================================================

  FE: CTTAAT-ACACTTTTGATTAACAAACTTTAAA
  PE: CT-AATTACACTTT-GAT-AACAAACT--AAA
  CL: CTTAATTACACTTT-GATTAACAAACTTTAAA

  IDX: 19  18  17  16  15  14  13  12  11  10  9   8   7   6   5   4   3   2   1   
  CAT: PE  FE  PE  PE  PE  PE  PE  PE  FE+ FE  FE  FE  FE  FE  FE  FE  FE  FE  FE  

  IDX: 0   
  CAT: FE

------------------------------------------------------------------------
Viewing Alignments
------------------------------------------------------------------------

TREAT can optionally store alignments into a database for more complex analysis
and viewing in a web browser. Here we present a slighlty more complex example
where we view sequences from ribosomal protein S12 (RPS12) from Trypanosoma
brucei mitochondria. 

We first create an empty database file to store our data. This command creates
a database file in your current working directory called "rps12.db"::

  $ treat -d `pwd`/rps12.db db create
  **DANGER AREA** Are you sure you want to proceed? [n] y
  You are about to create the ENTIRE database from scratch. Are you sure? [n] y

Next we create a FASTA file (templates.fasta) which contains the Fully Edited,
Pre-Edited and 1 Alternativley Edited template sequences::

  >RPS12-FE Fully Edited
  CTAATACACTTTTGATAACAAACTAAAGTAAAtAtAttttGttttttttGCGtAtGtGAT
  TTTTGtAtGGttGttGtttACGttttGttttAtttGttttAtGttAttAtAtGAGtCCGC
  GAttGCCCAGttCCGGtAACCGACGtGtAttGtAtGCCGtAttttAttTAtAtAAttttG
  tttGGAtGttGCGttGttttttttGttGttttAttGGtttAGttAtGTCAttAtttAttA
  tAGAGGGTGGtGGttttGttGAtttACCCGGtGTAAAGtAttAtACACGTAttGtAAGtt
  AGATTTAGAtATAAGATATGTTTTT
  >RPS12-PE Pre-Edited
  CTAATACACTTTTGATAACAAACTAAAGTAAAAAGGCGAGGATTTTTTGAGTGGGACTGG
  AGAGAAAGAGCCGTTCGAGCCCAGCCGGAACCGACGGAGAGCTTCTTTTGAATAAAAGGG
  AGGCGGGGAGGAGAGTTTCAAAAAGATTTGGGTGGGGGGAACCCTTTGTTTTGGTTAAAG
  AAACATCGTTTAGAAGAGATTTTAGAATAAGATATGTTTTT
  >RPS12-A0 Alternative Editing (Cruz-Reyes 2013) alt_start=27 alt_stop=34
  CTAATACACTTTTGATAACAAACTAAAGTAAAtAtAttttGttttttttGCGtAtGtGAT
  TTTTGtAtGGttGttGtttACGttttGttttAtttGttttAtGttAttAtAtGAGtCCGC
  GAttGCCCAGttCCGGtAACCGACGtGtAttGtAtGCCGtAttttAttTAtAtAAttttG
  tttGGAtGttGCGttGttttttttGttGttttAttGGtttAGttAtGTCAttAtttAttA
  tAGAGGGTGGtGGttttGttGAtttACCtCGttGGttTAtAtAGtAttAtACACGTAttG
  tAAGttAGATTTAGAtATAAGATATGTTTTT

And load the template data into the database::

  $ treat -d `pwd`/rps12.db load_templates -g RPS12 -f templates.fasta -v
  [INFO] Processing template sequences for gene: RPS12
  [INFO] Loaded 3 template sequences for gene RPS12
  [INFO] Done.

Next we create a file with our DNA fragment reads (sample-A.fasta)::

  >cl-1 merge_count=10
  CTAATACACTTTTGATAACAAACTAAAGATATAATATTTTTGTTTTTTTTGCGTATGTGA
  TTTTTGTATGGTTGTTGTTTACGTTTTGTTTTATTTGTTTTATGTTATTATATGAGTCCG
  CGATTGCCCAGTTCCGGTAACCGACGTGTATTGTATGCCGTATTTTATTTATATAATTTT
  GTTTGGATGTTGCGTTGTTTTTTTTGTTGTTTTATTGGTTTAGTTATGTCATTATTTATT
  ATAGAGGGTGGTGGTTTTGTTGATTTACCCGGTGTAAAGTATTATACACGTATTGTAAGT
  TAGATTTAGATATAAGATATGTTTTT
  >cl-2 merge_count=9
  CTAATACACTTTTGATAACAAACTAAAGTAAAAAGGCGAGGATTTTTTGAGTGGGATTCGGT
  ATTTGTTTTATGTTATTATATGAGTCCGCGATTGCCCAGCTCTGGTAACCGACGTGTATTGT
  ATGCCGTATTTTATTTATATAATTTTGTTTGGATGTTGCGTTGTTTTTTTTGTTGTTTTATT
  GGTTTAGTTATGTCATTATTTATTATAGAGGGTGGTGGTTTTGTTGATTTACCCGGTGTAAA
  GTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTT
  >cl-3 merge_count=120
  CTAATACACTTTTGATAACAAACTAAAGTAAAAAGGCGAGGATTTTTTGAGTGGGATTCGGTA
  TTTGTTTTATGTTATTATATGAGTCCGCGATTGCCCAGCTCTGGTAACCGACGTGTATTGTAT
  GCCGTATTTTATTTATATAATTTTGTTTGGATGTTGCGTTGTTTTTTTTGTTGTTTTATTGGT
  TTAGTTATGTCATTATTTATTATAGAGGGTGGTGGTTTTGTTGATTTACCCGGTGTAAAGTAT
  TATACACGTATTGTAAGTTAGATTTAGATATAACATATGTTTTT

And load the sample data using TREAT::

  $ treat -d `pwd`/rps12.db load_samples -g RPS12 -f sample-A.fasta -v
  **DANGER AREA** Delete all alignments for gene RPS12 replicate 0 ? [n] y
  [INFO] Computing average read count...
  [INFO] Total reads accross all files: 139
  [INFO] Normalizing to average read count: 139
  [INFO] Computing total read count for file: sample-A.fasta
  [INFO] Total reads for file: 139
  [INFO] Normalized scaling factor: 1.0
  [INFO] Processing fragments in file: sample-A.fasta
  [INFO] Using sample name: sample-A
  [INFO] Loaded 3 fragment sequences
  [INFO] Done.

We can now start the TREAT server and view the sequences in a web browser::

  $ treat -d `pwd`/rps12.db runserver -v
  [INFO]  * Running on http://127.0.0.1:5000/

.. image:: docs/treat-screen-shot.png

------------------------------------------------------------------------
License
------------------------------------------------------------------------

TREAT is released under the GNU General Public License ("GPL") Version 3.0.
See the LICENSE file. 
