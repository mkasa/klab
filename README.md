klab
====

This public repository contains a suite of programs that are
used in research in Kasahara lab.

The programs include various codes such as small utilities
or tools for bioinformatics.

Installation
------------
This package follows the standard installation process of waf::

	$ ./waf configure
	$ ./waf build
	$ ./waf install

If you like GNU autotools-style configure, you can instead do like this::

	$ ./configure
	$ make
	$ make install

although these commands are just a wrapper for the former commands.

When you wish to use scripts that use Perl modules, you need to add
`--enable-perl` for `waf configure`, namely::

	$ ./waf configure --enable-perl

or::

	$ ./configure --enable-perl

You may need to be root for installing the Perl libraries.

List
----

|name|type|description|note|
|---|---|---|---|
|fatt|bio|FASTA/FASTQ manipulation tool.|see also doc/fatt.md|
|sieve|general|Random sampling of text files.|maybe useful for cross-validation|
|convertsequence|bio|Format conversion of sequence files.|see also perldoc|
|fixshebang|general|Fix shebang lines.|see also perldoc|
|gcc-color|general|Color the output of gcc/g++.|not extensively used.|
|icc-color|general|Color the output of Intel C++.|
|mydaemon|general|Automatically set ups crontab to ensure your daemon is running. (Works without root)|
|rep|general|Simple wrapper for different VCS such as svn/git/hg.|see also perldoc|
|sq|general|Execute SQL queries over CSV files.|see also perldoc. requires a bunch of Perl modules. DBD::CSV is required.|
|mddoc|general|Simple wrapper to view formatted Markdown (and restructured) texts via text browser|Requires Markdown.pl or pandoc.|
|gmddoc|general|Simple wrapper to view formatted GitHub-flavored Markdown via (graphic) web browser|Requires grip (python module).|
|sha_scan|general|Find duplicated files by collecting SHA1 hashes of files in a given directory.|
|rep|general|Repository utilitiy (one command, manage svn/cvs/git/hg at once!)|
|gfwhere|general|Find inconsistently replicated files in GlusterFS|
|json2csv|general|Convert JSON into CSV|
|csv2html|general|Convert CSV into HTML|
|csv2md|general|Convert CSV into a table in Markdown extra|
|split_paf|bio|Split alignments(PAF format) by looking CIGAR string|

Licenses
--------
The programs are licensed under the modified BSD Licenses
unless otherwise stated in source code.

The princple is that 3rd party libraries and their derivatives
are basically licensed under their original licenses, while
what we developed from scratch are licensed under the modified BSD.
The 3rd party libraries include SQLite3 (http://www.sqlite.org/), sqdbcpp
(http://code.google.com/p/sqdbcpp/), imgcat2 (modified from imgcat, https://iterm2.com/3.2/documentation-images.html).

