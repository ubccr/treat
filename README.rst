===============================================================================
TREAT - Trypanosome RNA Editing Alignment Tool
===============================================================================

------------------------------------------------------------------------
About
------------------------------------------------------------------------

TREAT is a multiple sequence alignment and visualization tool specifically
designed to permit the user to analyze variation in sequences caused by Uridine
insertion/deletion RNA editing. This phenomenon occurs in trypanosomes, a group
of unicellular parasitic flagellate protozoa such as Trypanosoma brucei which
causes African sleeping sickness. The pre-mRNA sequences in trypanosomes are
posttranscriptionally edited by the insertion/deletion of uridylate residues.
TREAT aligns sequences using three bases and assembles editing sites to detect
the extent of editing of the fourth base, called the edit base. The edit base
is configurable in TREAT and by default uses 'T'. TREAT is written in Go and is
released under a BSD style free software license. 

------------------------------------------------------------------------
Installation
------------------------------------------------------------------------

Download the latest binary release `here <https://github.com/ubccr/treat/releases>`_::

  $ tar xvzf treat-x.x.x.tar.gz
  $ cd treat-x.x.x

Install from source::

  $ git clone https://github.com/ubccr/treat.git treat
  $ cd treat
  $ go build ./...

Install using go tools::

  $ go get install github.com/ubccr/treat/...

------------------------------------------------------------------------
Simple Alignment
------------------------------------------------------------------------

TREAT can perform a global alignment between two arbitrary sequences detecting
the amount of insertion/deletions (indels) of a single base (called the "edit
base"). The default edit base in TREAT is "T". For example::

  $ treat align -1 ATCTGTATGT -2 ATTCGATTG -b T
  A-TCTGTA-TGT
  ATTC-G-ATTG-

------------------------------------------------------------------------
Alignment with Templates
------------------------------------------------------------------------

TREAT can perform global alignments using template sequences.  TREAT requires
two user provided template sequences: fully edited and pre-edited. The fully
edited template represents a mature edited mRNA transcript (completely
precisely edited mRNA). The pre-edited template represents the sequence that
will be edited in the mature RNA. TREAT accepts in put sequences in FASTA
format. For example::

  # simple-templates.fa
  >Fully Edited
  CTTAATACACTTTTGATTAACAAACTTTAAA
  >Pre-Edited
  CTAATTACACTTTGATAACAAACTAAA

  # simple-sequences.fa
  >example-1
  CTTAATTACACTTTGATTAACAAACTTTAAA

We save the above sequence files and run the alignment using TREAT::

  $ treat align -t simple-templates.fa -f simple-sequences.fa
  ================================================================================
  example-1
  ================================================================================
  EditStop: 11
  JuncEnd: 18
  JuncLen: 7
  ================================================================================

  FE: CTTAA-TACACTTTTGATTAACAAACTTTAAA
  PE: C-TAATTACAC-TTTGA-TAACAAAC--TAAA
  CL: CTTAATTACAC-TTTGATTAACAAACTTTAAA

------------------------------------------------------------------------
Large Scale Alignment Analysis
------------------------------------------------------------------------

TREAT can optionally store alignments into a database for more complex
analysis, searching, and viewing in a web browser. TREAT has been tested on
RNA-Seq data containing millions of sequences reads. Here we present an
example using sequences from ribosomal protein S12 (RPS12) from Trypanosoma
brucei mitochondria. 

TREAT accepts sequencing data in FASTA format. An example FASTA file
(templates.fasta) containing the Fully Edited, Pre-Edited and one alternatively
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

  $ treat --db treat.db load -g RPS12 -f sample-1.fa -t templates.fa
  Total reads across all samples: 139
  Normalizing to average read count:: 139.0000
  Computing total read count for file: sample-1.fa
  Total reads for file: 139
  Normalized scaling factor: 1.0000
  Processing fragments for sample name : sample-1
  Loaded 3 fragment sequences for sample sample-1

A new database file has been created called "treat.db". We can now search the
data using the TREAT command line tool::

  $ treat --db treat.db search -g RPS12 -l 10 --csv
  gene,sample,norm,read_count,alt_editing,has_mutation,edit_stop,junc_end,junc_len,junc_seq
  RPS12,sample-1,10.0000,10,0,0,137,143,6,ATATAATATTTTTG
  RPS12,sample-1,9.0000,9,0,0,95,123,28,TTCGGTATTTGTTTTATGTTATTATATGAGTCCGCGATTGCCCAGCTCTG

Search options are described below::

  $ treat help search
  NAME:
     search - Search database

  USAGE:
     command search [command options] [arguments...]

  OPTIONS:
     --gene, -g                       Gene Name
     --sample, -s                     One or more samples
     --edit-stop "-1"                 Edit stop
     --junc-end "-1"                  Junction end
     --junc-len "-1"                  Junction len
     --alt "0"                        Alt editing region
     --offset, -o "0"                 offset
     --limit, -l "0"                  limit
     --has-mutation                   Has mutation
     --all, -a                        Include all sequences
     --has-alt                        Has Alternative Editing
     --csv                            Output in csv format
     --no-header, -x                  Exclude header from output

Start the TREAT server and view the sequences in a web browser::

  $ treat --db treat.db server -p 8080
  Computing edit stop site cache for gene RPS12...
  Using template dir: /path/to/treat
  Running on http://127.0.0.1:8080

.. image:: docs/treat-screen-shot.png

------------------------------------------------------------------------
License
------------------------------------------------------------------------

TREAT is released under a BSD style license. See the LICENSE file. 
