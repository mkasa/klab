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

<table style="border=1 solid">
	<tr><th>name</th><th>type</th><th>description</th><th>note</th></tr>
	<tr><td>fatt</td><td>bio</td><td>FASTA/FASTQ manipulation tool.</td><td>see also doc/fatt.md</td></tr>
	<tr><td>sieve</td><td>general</td><td>Random sampling of text files.</t><td>maybe useful for cross-validation</td></tr>
    <tr><td>convertsequence</td><td>bio</td><td>Format conversion of sequence files.</td><td>see also perldoc</td></tr>
    <tr><td>fixshebang</td><td>general</td><td>Fix shebang lines.</td><td>see also perldoc</td></tr>
    <tr><td>gcc-color</td><td>general</td><td>Color the output of gcc/g++.</td><td>not extensively used.</td></tr>
    <tr><td>icc-color</td><td>general</td><td>Color the output of Intel C++.</td><td></td></tr>
    <tr><td>mydaemon</td><td>general</td><td>Automatically set ups crontab to ensure your daemon is running. (Works without root)</td><td></td></tr>
    <tr><td>rep</td><td>general</td><td>Simple wrapper for different VCS such as svn/git/hg.</td><td>see also perldoc</td></tr>
    <tr><td>sq</td><td>general</td><td>Execute SQL queries over CSV files.</td><td>Obsolete. Use q instead.</td></tr>
    <tr><td>mddoc</td><td>general</td><td>Simple wrapper to view formatted Markdown (and restructured) texts via text browser</td><td>Obsolete. Use glow.</td></tr>
    <tr><td>gmddoc</t><td>general</td><td>Simple wrapper to view formatted GitHub-flavored Markdown via (graphic) web browser</td><td>Obsolete. Standalone grip can do this now.</td></tr>
    <tr><td>sha_scan</td><td>general</td><td>Find duplicated files by collecting SHA1 hashes of files in a given directory.</td><td></td></tr>
    <tr><td>rep</td><td>general</td><td>Repository utilitiy (one command, manage svn/cvs/git/hg at once!)</td><td></td></tr>
    <tr><rd>taw</td><td>general</td><td>Tiny Amazon EC2 Wrapper.</td><td>Moved to an independent repository</td></tr>
    <tr><td>gfwhere</td><td>general</td><td>Find inconsistently replicated files in GlusterFS</td><td></td></tr>
    <tr><td>json2csv</td><td>general</td><td>Convert JSON into CSV</td><td></td></tr>
    <tr><td>csv2html</td><td>general</td><td>Convert CSV into HTML</td><td></td></tr>
    <tr><td>csv2md</td><td>general</td><td>Convert CSV into a table in Markdown extra</td><td></td></tr>
    <tr><td></td><td></td><td></td><td></td></tr>
</table>

Licenses
--------
The programs are licensed under the modified BSD Licenses
unless otherwise stated in source code.

The princple is that 3rd party libraries and their derivatives
are basically licensed under their original licenses, while
what we developed from scratch are licensed under the modified BSD.
The 3rd party libraries include SQLite3 (http://www.sqlite.org/), sqdbcpp
(http://code.google.com/p/sqdbcpp/), imgcat2 (modified from imgcat, https://iterm2.com/3.2/documentation-images.html).

