.. _seqlib-configuration:

SeqLib configuration details
================================

Most parameters are specified within SeqLib objects. Experiment, Condition, and Selection objects have only a name (and output directory if at the root). :ref:`analysis-options`, such as scoring method, are chosen at run time.

Sequencing libraries have :ref:`general-seqlib-parameters`, :ref:`sequence-file-seqlib-parameters`, and other parameter groups depending on the type: 

+----------------------+---------+---------+------------+---------+
| SeqLib type          | Barcode | Variant | Identifier | Overlap |
+======================+=========+=========+============+=========+
| Barcoded Variant     | X       | X       |            |         |
+----------------------+---------+---------+------------+---------+
| Barcoded Identifier  | X       |         | X          |         |
+----------------------+---------+---------+------------+---------+
| Overlap              |         | X       |            | X       |
+----------------------+---------+---------+------------+---------+
| Basic                |         | X       |            |         |
+----------------------+---------+---------+------------+---------+
| Barcodes Only        | X       |         |            |         |
+----------------------+---------+---------+------------+---------+
| Identifiers Only     |         |         | X          |         |
+----------------------+---------+---------+------------+---------+

See :ref:`intro-seqlibs` for descriptions of each type.

.. _general-seqlib-parameters:

General parameters
--------------------------------

* Name

	The object name should be short, descriptive, and not conflict with other object names in the analysis.

* Output Directory
	
	Path to the output directory. This field only appears for the root object.

* Time Point

	The time point must be an integer. All Selections require an input library as time point 0. Time point values may refer to the round of selection or hour of sampling.

* Counts File

	Required for Counts File Mode. Path to an HDF5 file or tab-separated value file that contains counts for this time point. Raw counts from that file will be used for this SeqLib. If an HDF5 file is provided, all tables in the "raw/" group are copied. Sequence file parameters will be ignored. The file must have the suffix ".h5" for HDF5 or one of ".txt" ".tsv" or ".csv" for tab-separated value files.

	.. note:: Tab-separated value files must have exactly two columns separated by a tab. The first line of the file must have the column heading "counts" preceded by a single tab character. The first column contains the barcode, identifier, or HGVS variant string depending on the type of raw counts required by the SeqLib type. The second column contains the count for that element.

.. _sequence-file-seqlib-parameters:

Sequence file parameters
--------------------------------

Enrich2 accepts sequence files in FASTQ_ format. These files may be processed while compressed with gzip or bzip2. The file must have the suffix ".fq" or ".fastq" before compression. 

* Reads

	Required for FASTQ_ File Mode. Path to a FASTQ_ file containing the sequencing reads. For overlap SeqLibs, there are fields for Forward Reads and Reverse Reads.

* Reverse

	Checking this box will reverse-complement reads before analysis. Not present for Overlap SeqLibs.

Read filtering parameters
++++++++++++++++++++++++++++++++++

Filters are applied after read trimming and any read merging.

* Minimum Quality

	Minimum single-base quality. If a single base in the read has a quality score below this value, the read will be discarded.

* Average Quality

	Average read quality. If the mean quality score of all bases in the read is below this value, the read will be discarded.

* Maximum N's

	Maximum number of N nucleotides. If the read contains more than this number of bases called as N, the read will be discarded. This should be set to 0 in most cases.

* Remove Unresolvable Overlaps

	Present for Overlap SeqLibs only. Checking this box discards merged reads with unresolvable discrepant bases (see :ref:`overlap-seqlib-parameters`).

* Maximum Mutations

	Present for SeqLibs with variants only. Maximum number of mutations. If the variant contains more than this number of differences from wild type, the variant is discarded (or aligned if that option is enabled under :ref:`variant-seqlib-parameters`).

.. _barcode-seqlib-parameters:

Barcode parameters
--------------------------------

* Barcode-variant File

	Not present for barcode-only SeqLibs. Path to a tab-separated file in which each line contains a barcode followed by its identifier or linked variant DNA sequence. This file may be processed while compressed with gzip or bzip2. 

* Minimum Count

	Minimum barcode count. If the barcode has fewer counts than this value, it will not be scored and will not contribute to counts of its variant or identifier.

* Trim Start

	Position of the first base to keep when trimming barcodes. All subsequent bases are kept if Trim Length is not specified. Reverse-complementing occurs before trimming. Bases are numbered starting at 1.

* Trim Length

	Number of bases to keep when trimming barcodes. Starts at the first base if Trim Start is not specified. Reverse-complementing occurs before trimming.

.. _variant-seqlib-parameters:

Variant parameters
--------------------------------

* Wild Type Sequence
	
	The wild type DNA sequence. This sequence will be compared to reads or the barcode-variant map when calling variants. All sequences must have the same length and starting position.

* Wild Type Offset

	Integer added to every variant nucleotide position. Used to place variants in the context of a larger sequence.

* Protein Coding

	Checking this box will interpret the wild type sequence as protein coding. The wild type sequence must be in frame.

* Use Aligner

	Checking this box will enable Needleman-Wunsch alignment. Insertion and deletion events will be called.

.. warning:: Using the aligner will dramatically increase run time, and is not recommended for most users.

* Minimum Count

	Minimum variant count. If the variant has fewer counts than this value, it will not be scored and will not contribute to counts of any synonymous elements.

.. _identifier-seqlib-parameters:

Identifier parameters
--------------------------------

* Minimum Count

	Minimum identifier count. If the identifier has fewer counts than this value, it will not be scored.

.. _overlap-seqlib-parameters:

Overlap parameters
--------------------------------

Overlapping read pairs reduce the likelihood of calling sequencing errors as variants. Paired-end Illumina reads are generated such that they overlap in the target region.

When Enrich2 combines forward and reverse reads into merged reads, base quality values in the overlapping region are defined as the higher quality value at each position. Mismatches are resolved by assuming the base with the higher quality value is correct. If mismatched bases have the same quality value, the position is considered unresolvable and replaced by an 'X' base.

* Forward Start

	Position of the first overlapping base in the forward read. Bases are numbered starting at 1. 

* Reverse Start

	Position of the first overlapping base in the reverse read before reverse complementing. Bases are numbered starting at 1. 

* Overlap Length

	Number of bases in the overlapping region.

* Maximum Mismatches

	Maximum number of mismatches in the overlapping region. If a merged read has more than this number of mismatches, the read pair will be discarded.

* Overlap Only

	Checking this box will trim the merged reads to the overlapping region.

