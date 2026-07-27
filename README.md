klab
====

This public repository contains a suite of programs that are
used in research in Kasahara lab.

The programs include various codes such as small utilities
or tools for bioinformatics.

Installation
------------
This package follows the standard installation process of waf:

```console
$ ./waf configure
$ ./waf build
$ ./waf install
```

If you like GNU autotools-style configure, you can instead do like this:

```console
$ ./configure
$ make
$ make install
```

although these commands are just a wrapper for the former commands.
Note that out-of-tree builds are not supported; run `./configure` from the
top of the source tree.

When you wish to use scripts that use Perl modules, you need to add
`--enable-perl` for `waf configure`, namely:

```console
$ ./waf configure --enable-perl
```

or:

```console
$ ./configure --enable-perl
```

You may need to be root for installing the Perl libraries.

Dependencies
------------
The Python scripts are self-contained: each one starts with
`#!/usr/bin/env -S uv run --script` and carries PEP 723 inline metadata, so
[uv](https://docs.astral.sh/uv/) fetches what they need on first run and you
have nothing to install.

The Perl scripts do have external dependencies. There are two ways to get them:

* `cpanfile` lists the non-core CPAN modules. Use it with
  [Carton](https://metacpan.org/dist/Carton) (`carton install`, then
  `carton exec -- ./script/...`) or with
  [cpm](https://metacpan.org/dist/App-cpm) (`cpm install -L local`).
* `environment.yml` is a conda/mamba environment (`conda env create -f
  environment.yml`). This is usually the better choice, because in addition to
  the Perl modules it also provides the **external binaries** the scripts shell
  out to -- BLAST, blat, primer3, seqkit, jq, R, xz -- and a compiler toolchain
  for building `fatt` and `sieve`. Carton and cpm cannot install those.

Both files say in comments which script needs each entry, and which few
dependencies cannot be satisfied from either source.

List
----

|name|type|description|note|
|---|---|---|---|
|fatt|bio|FASTA/FASTQ manipulation tool.|see also doc/fatt.md|
|sieve|general|Random sampling of text files.|maybe useful for cross-validation|
|convertsequence|bio|Format conversion of sequence files.|see also perldoc|
|fixshebang|general|Fix shebang lines.|see also perldoc|
|gcc-color|general|Color the output of gcc/g++.|not extensively used.|
|icc-color|general|Color the output of Intel C++.|see also perldoc|
|mydaemon|general|Automatically set ups crontab to ensure your daemon is running. (Works without root)||
|rep|general|Repository utility; a single wrapper command to manage svn/cvs/git/hg at once.|see also perldoc|
|sha_scanp|general|Find duplicated files by collecting SHA-256 hashes of files in a given directory.||
|gfwhich|general|Show which GlusterFS hosts hold the files in the current directory, and find inconsistently replicated ("split-brain") files.|see also perldoc|
|json2csv|general|Convert JSON into CSV||
|csv2html|general|Convert CSV into HTML||
|split_paf|bio|Split alignments (PAF format) by looking CIGAR string||
|cco|general|Generate a conventional-commit message from the staged changes; can install/remove a `prepare-commit-msg` hook.|see also perldoc|
|headtail|general|Show the first and the last few lines of each input (head and tail in one command).|see also perldoc|
|imgcat2|general|Display images inline on terminals that support iTerm2's inline image protocol, with control over the display size.|see also doc/imgcat2.md|
|ispcr|bio|In-silico PCR; finds the fragments a primer pair would amplify from target sequences.|needs BLAST and BioPerl; see also perldoc|
|plotr|general|Draw a plot with R and display it on the terminal (iTerm2).|see also perldoc|
|recompressbyxz|general|Recompress gzip/bzip2 (or plain) files with xz, verifying the round trip before the original is replaced.||
|reduce_genome_feature|bio|Reduce genome features (GC fraction, N fraction, BED coverage, gene counts) into fixed-size bins and emit them as TSV.|`--help` for usage|
|pbjellysummary2json|bio|Convert a PBJelly summary into JSON.||
|quastreport2json|bio|Convert a QUAST contig report into JSON.||
|BLASTM8Parse.pm|bio|Perl module that parses BLAST -m8 (tabular) output.|installed only with `--enable-perl`|

The `undocumentedscript/` directory holds 21 further Perl files that are kept
in the repository for reference only: they are neither built, installed, nor
supported, and they are not covered by the list above.

Licenses
--------
The programs are licensed under the modified BSD Licenses
unless otherwise stated in source code.

The princple is that 3rd party libraries and their derivatives
are basically licensed under their original licenses, while
what we developed from scratch are licensed under the modified BSD.
The 3rd party libraries include SQLite3 (http://www.sqlite.org/), sqdbcpp
(http://code.google.com/p/sqdbcpp/), imgcat2 (modified from imgcat, https://iterm2.com/3.2/documentation-images.html).

