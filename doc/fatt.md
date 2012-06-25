fatt
====

This tool manipulates FASTA/FASTQ files.

extract
-------
You can extract sequences with the specified names::

    fatt extract --seq chr1 foo.fastq > chr1.fastq

You can give more than one name::

    fatt extract --seq chr3 --seq chr2 foo.fastq > chr2and3.fastq

If you have millions of names, use the file input::

    generate-names-by-yourprogram > names.txt
    fatt extract --file names.txt foo.fastq > picked.fastq

If you do not like temporary files, than use stdin input mode::

    generate-names-by-yourprogram | fatt extract --stdin foo.fastq > picked.fastq

You can output sequences OTHER THAN the specified::

    fatt extract --reverse --seq chr1 foo.fastq > withoutchr1.fastq

If a file contains some duplicated sequences and you want to eliminate them::

    fatt extract --unique sequences_with_lots_of_duplication.fastq > uniqueseqs.fastq

If both --start=x and --end=y options are specified, you can obtain sequences ranging from x-th (inclusive) to y-th (exclusive)::

    fatt extract --start=200 --end=300 foo.fastq > foo_200_to_300.fastq

This feature might be useful especially for n-fold cross-validation. Alternatively, you can do the same thing by::

    fatt extract --start=200 --num=100 foo.fastq > foo_200_to_300.fastq

Note that --start and --end take 0-origin numbers.

count
-----
You can count the number of the sequences in each given file.

    fatt count foo.fastq

name
----
You can output the name of the sequences in each given file.

    fatt name foo.fastq > foo.names.txt

chksamename
-----------
You can get the name of the sequences that appear more than
once in the given files. Note that it looks only names,
so differnt sequences with the same name will also be reported.

    fatt chksamename duplication.fastq

When you find some duplicated sequences, you probably want to use 'fatt extract --unique'
to filter out the duplicated sequences.

len
---
You can calculate the length of the sequences in given files.

    fatt len foo.fasta

index
-----
It creates an index on the name of the sequences in each given file.
Subsequent access may get faster if the file is very large and you
retrieve only a few sequences.

    fatt index foo.fasta

Index files are SQLite3 database that contains the name, the position, and the rank
of the sequences in each given file. The file name of the index is the original
FASTA/FASTQ file name plus '.index'. For example, the above example creates foo.fasta.index.
Therefore, the directory must be writable. It does not overwrite if there is any existing file.
This command accesses storage quite randomly, so avoid using remote file systems 
for performance where possible.

guessqvtype
----------
There are several types of FASTQ formats. They differ in the base of Quality Value.
This command takes FASTQ files and guesses the base of QV.

    fatt guessqvtype foo.fastq

tocsv
-----
It converts input FASTA/FASTQ files into CSV files.

    fatt tocsv foo.fastq > foo.csv

If you like TSV instead of CSV, give --tsv.

    fatt tocsv --tsv foo.fastq > foo.tsv

In the both cases, you can specify --noheader to suppress the header output.


fold
----
If the sequences or the QVs are too long in a single line, you can fold at
the specified length (70 chars by default).

    fatt fold foo_long_lines.fastq > foo_folded.fastq

Give --len=n to fold at n characters.

    fatt fold --len=50 foo_long_lines.fastq > foo_folded.fastq

unfold
------
It collects nucleotide characters into a single line. Most Illumina reads are already in this format.

    fatt unfold foo.fastq > foo_unfolded.fastq


help
----
You can see the description of a subcommand. For example, if you do not remember
the options for 'fatt extract', than you would type

    fatt help extract
