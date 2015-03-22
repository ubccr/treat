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

  $ treat align -t simple-templates.fasta -f simple-sequences.fasta -b t
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

------------------------------------------------------------------------
Viewing Alignments
------------------------------------------------------------------------

TREAT can optionally store alignments into a database for more complex analysis
and viewing in a web browser. Here we present a slighlty more complex example
where we view sequences from ribosomal protein S12 (RPS12) from Trypanosoma
brucei mitochondria. 

Treat accepts sequencing data in FASTA format. An example FASTA file
(templates.fasta) containing the Fully Edited, Pre-Edited and one Alternativley
Edited template sequences is shown below::

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

FASTA file with our DNA fragment reads (sample-A.fasta)::

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

Load the sample data using TREAT::

  $ treat load -g RPS12 -f sample-A.fasta -t templates.fasta
  Total reads across all samples: 350
  Normalizing to average read count:: 350.0000
  Computing total read count for file: sample-A.fa
  Total reads for file: 350
  Normalized scaling factor: 1.0000
  Processing fragments for sample name : sample-A
  Loaded 18 fragment sequences for sample sample-A

We can now start the TREAT server and view the sequences in a web browser::

  $ treat -p 8080

.. image:: docs/treat-screen-shot.png

------------------------------------------------------------------------
License
------------------------------------------------------------------------

TREAT is released under the GNU General Public License ("GPL") Version 3.0.
See the LICENSE file. 
