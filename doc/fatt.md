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

Index files are SQLite3 database that contains the name and position of the sequences
in each given file. The file name of the index is the original FASTA/FASTQ file name
plus '.index'. For example, the above example creates foo.fasta.index. Therefore,
the directory must be writable. It does not overwrite if there is any existing file.
This command accesses storage quite randomly, so avoid using remote file systems 
for performance where possible.

help
----
You can see the description of a subcommand. For example, if you do not remember
the options for 'fatt extract', than you would type

    fatt help extract
