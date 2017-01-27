fatt
=====

This tool manipulates FASTA/FASTQ files.

extract
--------
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
------
You can count the number of the sequences in each given file.

    fatt count foo.fastq

name
-----
You can output the name of the sequences in each given file.

    fatt name foo.fastq > foo.names.txt

chksamename
------------
You can get the name of the sequences that appear more than
once in the given files. Note that it looks only names,
so differnt sequences with the same name will also be reported.

    fatt chksamename duplication.fastq

When you find some duplicated sequences, you probably want to use 'fatt extract --unique'
to filter out the duplicated sequences.

len
----
You can calculate the length of the sequences in given files.

    fatt len foo.fasta

stat
-----
You can show the statistics of input sequences by the following command.

    fatt stat foo.fasta foo2.fastq

If you give multiple input files, the input files are considered as a single (large) file.
When you give --html option, the output is formatted in HTML.

    fatt stat --html foo.fasta

In addition, you may use --json to output in JSON format.

    fatt stat --json foo.fasta

index
------
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
------------
There are several types of FASTQ formats. They differ in the base of Quality Value.
This command takes FASTQ files and guesses the base of QV.

    fatt guessqvtype foo.fastq

[A Wikipedia article about FASTQ
(en)](http://en.wikipedia.org/wiki/FASTQ_format) explains variants of
FASTQ formats.

convertqv
----------
Convert FASTQ files into anothers with different QV base/range.
This would help you convert old Solexa/Illumina FASTQ into a Sanger
FASTQ file, or vice versa.

    fatt convertqv --fromillumina --tosanger foo.fastq > foo.sanger.fastq

tocsv
------
It converts input FASTA/FASTQ files into CSV files.

    fatt tocsv foo.fastq > foo.csv

If you like TSV instead of CSV, give --tsv.

    fatt tocsv --tsv foo.fastq > foo.tsv

In the both cases, you can specify --noheader to suppress the header output.

tofasta
--------
It converts input FASTQ files into FASTA files.

    fatt tofastq foo.fastq > foo.fasta

If the input is not in FASTQ format, it aborts with an error message.

fold
-----
If the sequences or the QVs are too long in a single line, you can fold at
the specified length (70 chars by default).

    fatt fold foo_long_lines.fastq > foo_folded.fastq

Give --len=n to fold at n characters.

    fatt fold --len=50 foo_long_lines.fastq > foo_folded.fastq

unfold
-------
It collects nucleotide characters into a single line. Most Illumina reads are already in this format.

    fatt unfold foo.fastq > foo_unfolded.fastq

composition
------------
It calculates all 1- to 3-mer frequencies in the given files.

    fatt composition foo.fastq

It recognizes 'A' and 'a' as different characters by default. To ignore
cases, add '--ignorecase'.

split
------
Splits (possibly) huge files into smaller chunks of files.

    fatt split --num=10 huge.fastq

huge.fastq will be split into 10 files of similar sizes.
Alternatively you can tell the maxinum number of bases for a single
output file. The following example splits huge.fastq into files of
equal to or slightly larger than 10 Mbp (except for the last file).

    fatt split --max=10000000000 huge.fastq

fatt counts all nucleotides by default, but you can tell it to ignore
N's (ignorecase) when you give --excn option. This might help you when
the file contains lots of N's.

It usually outputs the exact number of files, but sometimes it cannot;
such an extreme example is giving --num=100 for human chromosomes, for
which the number of the sequences is less than 100.

To get the exact number of the output files, you can do one of the
followings.

Output files are named prefix.1, prefix.2, ..., and so on. You can check
whether they exist from prefix.1 and you will find the number if you
failed to find a file with a particular suffix.

Another way to do that is giving --retstat. With this, fatt will return
the number by the exit code. The following loop in bash will process the split
files. The exit code is only 7 bits in width, so you cannot use this
method when the number of partitions exceeds 100.

    fatt split --num=10 --retstat huge.fastq
    for i in 1..$?; do
        do something with huge.fastq.$i
    done

Last but not the least, you can use --retfile option to return the
number of partitions by file. The following example might tell you in a
second.

    fatt split --num=1000 --filestat splitnum.txt huge.fastq
    for i in 1..`cat splitnum.txt`; do
        do something with huge.fastq.$i
    done

edit
-----
It allows us to edit FASTA/FASTQ files. It takes a Genome Edit Script
(GES) and optional input files.

    fatt edit edit.ges input1.fastq input2.fastq input3.fastq

The file type of inputs is automatically determined.

    fatt edit edit.ges input1.fa input2.fa input3.fa

Note that fatt looks the content of input files, but not the file
extensions. In other words, if a file 'input.fastq' contains FASTA
sequences, it is considered as a FASTA file.

Next we explain the format of GES.
GES is similar to shell script, but it differs in commands.
It ignores blank lines and lines starting with '#':

    # You can write comment
    
    # Blank line is ignored

The other lines contain commands.
Each line starts with a command. Arguments for that command (if any) follow it:

    command arg1 arg2 arg3 ...

If you give an argument containing space characters, quote it with '"'.
Here is an example:

    setdesc "Chromosome I" "E. coli MG1655 Chromosome I (circular)"

Available commands are the following.

### loadall
You can load an entire file (FASTA/FASTQ).

    loadall input.fastq

Note that you cannot use both FASTA and FASTQ files. It means that if
the first file loaded is FASTQ, then the next file must be also FASTQ.
This restriction applies to other loading commands, too.

### loadone
It loads a single sequence from a file (FASTA/FASTQ).
This operation requires that the target input file has an index.
See 'fatt index' for details.

    loadone input.fastq seq1

It utilizes the index of the input file, achieving very fast access.

### saveall
You can save the entire sequences in memory to a file.
The type of the format (FASTA or FASTQ) is automatically determined.

    saveall output.fastq

Again, note that the output format is determined by the content, not by
the file name.

### saveone
It saves a single sequence into a file.

    saveone output.fastq seq1

### rename
It renames a sequence into another.

    rename seq1 new_seq1

### trim3 and trim5
It trims the 5'- or 3'-end of the specified sequence by the specified
amount (in bp).

    # Trims the 3'-end of seq1 by 3 bp.
    trim3 seq1 3
    # Trims the 5'-end of seq1 by 10 bp.
    trim5 seq2 10

### dupseq
It duplicates the specified sequence.

    dupseq seq1 newseq1

### split
It splits the specified sequence by the specified position.

    # seq1 is split into two sequences, seq1_left and seq1_right.
    split seq1 10 seq1_left seq1_right
    # seq1_left is 10 bp in length, and the rest goes to seq1_right
    # seq1 is removed.

### join
It joins two sequences (the first and the second arguments) into one
(the third argument).

    # seq1 and seq2 are joined
    join seq1 seq2 seq12
    # ex)
    #    seq1: AAA
    #    seq2: CCC
    #      =>
    #         seq12: AAACCC

The two sequences are removed after the new sequence is created.

### setdesc
It sets the description of the given sequence.

    setdesc seq1 "Chromosome I"

### print
It prints the content of the given sequence.

    print seq1

You can specify the range.

    # This will print nuleotides from 10 bp (0-origin, inclusive) to 20 bp (0-origin, exclusive)
    print seq1 10 20

help
----
You can see the description of a subcommand. For example, if you do not remember
the options for 'fatt extract', then you would type

    fatt help extract
