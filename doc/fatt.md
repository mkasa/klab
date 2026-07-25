fatt
=====

This tool manipulates FASTA/FASTQ files.

Compressed input
-----------------
Every command accepts gzip-compressed files in place of plain FASTA/FASTQ files.
Compression is detected automatically from the file content (not the extension),
so you may freely mix plain and gzip-compressed files::

    fatt count foo.fastq bar.fastq.gz

Index-based access (see ``index`` below) also works on gzip-compressed files.
Note, however, that seeking into a plain gzip file is inherently slow, because
the data up to the requested sequence has to be decompressed from the start.

For random access, prefer **BGZF** (the blocked-gzip format produced by
``bgzip``/htslib, used by samtools/tabix). A BGZF file is still an ordinary gzip
file, so every command reads it transparently, but it is split into independent
~64 KiB blocks. When you ``fatt index`` a BGZF file, fatt stores each sequence's
BGZF *virtual offset*, so ``fatt extract`` jumps straight to the right block and
decompresses only that block instead of the whole prefix. If a ``bgzip`` ``.gzi``
index sits next to the file, fatt reads it to locate the blocks; otherwise it
derives the block boundaries directly from the file. (Indexes created by older
versions of fatt on a BGZF file still work, but fall back to the slow seek with a
note suggesting you recreate them with ``fatt index --force``.)

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

Note that --start and --end take 0-origin numbers. You may also give --start
alone (from the x-th sequence to the end of the file) or --end alone (from the
beginning of the file to the y-th sequence, exclusive).

If an index file (see 'index' below) sits next to the input, fatt uses it
automatically to jump straight to the requested sequences.

    --index      Create the index (<file>.index) if there is none, and use it.
    --noindex    Do not use an index even if <file>.index exists.
    --force      With --index, replace an existing index.

The result does not depend on whether an index is used; the index only makes the
access faster. An index cannot be used together with --reverse (and therefore
not with --unique either), because those have to read the whole file anyway.

    fatt extract --index --seq chr1 foo.fastq > chr1.fastq
    fatt extract --noindex --seq chr1 foo.fastq > chr1.fastq

fatt refuses to use an index whose recorded offset does not point at the
requested sequence (i.e. a stale index), and it tells you to recreate it with
'fatt index --force'.

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

The duplicated names are written to the standard output, so you can save them:

    fatt chksamename duplication.fastq > dups.txt

Give --read to output only the read names (the file names are omitted).

When you find some duplicated sequences, you probably want to use 'fatt extract --unique'
to filter out the duplicated sequences.

len
----
You can calculate the length of the sequences in given files.

    fatt len foo.fasta

Give --name to put the name of each sequence in front of its length
(i.e. '<name><TAB><length>').

    fatt len --name foo.fasta

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
Give --force to remove an existing index first.

    fatt index --force foo.fasta

Sequence names must be unique, because the name is the primary key of the index.
If a name occurs twice, fatt tells you which name it was and creates no index
(use 'fatt chksamename' to list all of them).

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

    fatt tofasta foo.fastq > foo.fasta

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

clean
------
It replaces every character that is not a nucleotide with a nucleotide, which
is handy when a downstream tool only understands A/C/G/T.

    fatt clean foo.fasta > foo_clean.fasta

By default, [^ACGTNacgtn] is replaced and N/n is left alone; give --processn to
replace N/n as well. The replacement character is chosen by exactly one of the
following options (--random is the default):

    --a         Replace with 'A'
    --c         Replace with 'C'
    --g         Replace with 'G'
    --t         Replace with 'T'
    --n         Replace with 'N'
    --random    Replace with A/C/G/T at random
    --processn  Also replace N/n

    fatt clean --n --processn foo.fasta > foo_all_n.fasta

Quality strings of a FASTQ input are copied through untouched.

composition
------------
It calculates all 1- to 3-mer frequencies in the given files.

    fatt composition foo.fastq

It recognizes 'A' and 'a' as different characters by default. To ignore
cases, add '--ignorecase'.

The following options are available.

    --ignorecase  Treat 'A' and 'a' as the same character
    --monomer     Show only the 1-mer statistics
    --bimer       Show only the 2-mer statistics
    --trimer      Show only the 3-mer statistics
    --dapicheck   Show DAPI-staining related statistics instead
    --countends   Also count the n-mers that run over the ends of the sequences;
                  such an end is shown as '*'

split
------
Splits (possibly) huge files into smaller chunks of files.

    fatt split --num=10 huge.fastq

huge.fastq will be split into 10 files of similar sizes.
Alternatively you can tell the maxinum number of bases for a single
output file. The following example splits huge.fastq into files of
equal to or slightly larger than 10 Mbp (except for the last file).

    fatt split --max=10000000 huge.fastq

Note that --max counts BASES, not bytes.

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

Last but not the least, you can use --filestat option to return the
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

### delete
It removes the specified sequence from memory.

    delete seq1

### complement
It constructs the reverse complement of the specified sequence (the first
argument) and names it as the second argument.

    complement seq1 seq1_rc

Note that the source sequence is REMOVED, just like 'split' and 'join' remove
their inputs. If you want to keep it, duplicate it first:

    dupseq seq1 seq1_copy
    complement seq1_copy seq1_rc

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

Exit status
------------
fatt returns 0 when everything went fine, and a non-zero status when anything
failed (an unknown command, an input file that could not be opened or that is
neither FASTA nor FASTQ, a broken or stale index, a requested read name that is
not in the file, an invalid quality value, and so on). Diagnostics always go to
the standard error, so a failure never gets mixed into the output.

    fatt extract --seq chr1 foo.fasta > chr1.fasta || echo "extraction failed"

Two subcommands use the exit status for their own purpose: 'fatt split
--retstat' returns the number of the output files, and 'fatt edit' returns 2
when the Genome Edit Script fails.
