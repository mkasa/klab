# cpanfile -- CPAN dependencies of the Perl scripts in this repository.
#
# The Python scripts in script/ are self-contained: they carry PEP 723 inline
# metadata and run under `#!/usr/bin/env -S uv run --script`, so they need no
# dependency management at all.  The Perl scripts had none until this file.
#
# Only NON-CORE modules are listed.  Modules that ship with perl itself
# (JSON::PP, Getopt::Long, File::Temp, File::Spec, File::stat, Pod::Usage,
# Term::ANSIColor, Carp, Cwd, FileHandle, Fcntl, FindBin, ...) are deliberately
# omitted; core-ness was checked with Module::CoreList.
#
# The in-repo modules -- BLASTM8Parse, BLASTCSVParse, MiscLib, Primer3,
# ParallelExec -- are LOCAL and are resolved through `use lib $FindBin::Bin`.
# They are not CPAN distributions and must never be listed here.
#
#
# USAGE
# -----
# With Carton (installs into ./local, records exact versions in cpanfile.snapshot):
#
#     carton install
#     carton exec -- ./script/convertsequence --help
#
# With cpm (much faster, no snapshot):
#
#     cpm install -L local
#     PERL5LIB=$PWD/local/lib/perl5 ./script/convertsequence --help
#
# Optional feature groups (see the `feature` blocks below) are off by default:
#
#     carton install --with-feature=daemon
#     cpm install -L local --with-feature=daemon --with-feature=undocumented
#
# Neither Carton nor cpm can install the EXTERNAL BINARIES these scripts shell
# out to (BLAST, blat, primer3, seqkit, jq, R, ...).  If you want those too, use
# environment.yml with conda/mamba instead -- see the Dependencies section of
# README.md.
#
#
# A NOTE ON BioPerl
# -----------------
# BioPerl was split into many small distributions during the 1.7 series, so the
# modules these scripts use no longer all come from one place:
#
#     Bio::Seq, Bio::SeqIO                     -> BioPerl            (1.7.8)
#     Bio::DB::GenBank, Bio::DB::Query::GenBank-> Bio-DB-NCBIHelper  (1.7.8)
#     Bio::Perl                                -> Bio-Procedural     (1.7.4)
#     Bio::Seq::SeqWithQuality                 -> gone from CPAN; no longer used
#     Bio::Tools::RestrictionEnzyme            -> gone from CPAN, see below
#
# Version numbers in the Bio* world are three-part ("1.7.8"), which version.pm
# parses as a v-string (1.007008), so the minimum-version comparisons below
# behave as one would expect.

requires 'perl', '5.010';

# --- Required by installed scripts (script/, i.e. what `waf install` puts on PATH) ---

# script/convertsequence, script/ispcr
# Provides Bio::Seq and Bio::SeqIO.  1.7.2 is the oldest release with the
# current module layout; there is no reason to demand 1.7.8 specifically.
requires 'Bio::SeqIO', '1.7.2';
requires 'Bio::Seq',   '1.7.2';

# script/convertsequence (`use Bio::DB::GenBank`, and `require
# Bio::DB::Query::GenBank` on the --query path).  These moved OUT of the BioPerl
# distribution into Bio-DB-NCBIHelper, so BioPerl alone is not enough.
requires 'Bio::DB::GenBank', '1.7.2';

# NOTE: Bio::Seq::SeqWithQuality is NOT required.
#
# script/ispcr, undocumentedscript/blastm8alignments2alignmentsummary.pl and
# undocumentedscript/blastreport2dotplot.pl each used to carry a
# `use Bio::Seq::SeqWithQuality;` line.  That module was dropped from BioPerl
# after the 1.6 series and is not on CPAN in any distribution today (404 on
# MetaCPAN; the successor is Bio::Seq::Quality).  It would have pinned the whole
# repository to the legacy monolithic BioPerl 1.6.924.
#
# All three imports turned out to be DEAD: the module was named in the `use`
# line and never referenced anywhere else in those files.  The lines have been
# removed, so the scripts now run against a current BioPerl and nothing needs to
# be listed here.

# --- Optional: script/mydaemon's daemon mode ---
#
# `feature` rather than `recommends` on purpose.  Proc::PID::File is loaded
# lazily -- script/mydaemon does `unless(eval { require Proc::PID::File; 1 })`
# and degrades gracefully -- and it is only reached on the daemon code path, so
# it is genuinely optional rather than merely "nice to have".  A feature block
# says exactly which capability it buys and can be opted into by name
# (`--with-feature=daemon`), whereas `recommends` is an unnamed all-or-nothing
# flag (`--with-recommends`) that would drag in every other recommendation too.
feature 'daemon', 'script/mydaemon: PID-file locking for --daemon mode' => sub {
    requires 'Proc::PID::File', '1.27';
};

# --- Optional: undocumentedscript/ ---
#
# README.md states that undocumentedscript/ is kept "for reference only: they
# are neither built, installed, nor supported".  Its dependencies are therefore
# behind a feature flag so that a normal `carton install` does not pay for them.
feature 'undocumented', 'undocumentedscript/: unsupported reference scripts' => sub {

    # undocumentedscript/xmlex.pl
    # 2.44 is the version that has been stable for a decade and is what both
    # bioconda and conda-forge ship; 2.59 is current but is not required.
    requires 'XML::Parser', '2.44';

    # undocumentedscript/blastreport2dotplot.pl
    # One distribution supplies both PostScript::Simple (0.09) and
    # PostScript::Simple::EPS (0.02); requiring the former pulls in the latter.
    requires 'PostScript::Simple', '0.09';

    # undocumentedscript/blastreport2dotplot.pl
    # `use Bio::Perl` -- the procedural BioPerl facade, which now lives in its
    # own Bio-Procedural distribution rather than in BioPerl.
    requires 'Bio::Perl', '1.7.4';

    # NOTE: Bio::Tools::RestrictionEnzyme is NOT required.  Its only consumer was
    # undocumentedscript/blastcsv2rflpmap.pl, which has been deleted.  The module
    # was removed from BioPerl before 1.6.924 and is 404 on MetaCPAN; its
    # successors (Bio::Restriction::*) live only in the bioperl/Bio-Restriction
    # GitHub repository and have never been released to CPAN, so the script could
    # not have been made installable without being ported.

    # NOTE: undocumentedscript/ParallelExec.pm and parablast.pl / parablat.pl do
    # `eval { require TGEW; 1 }`.  TGEW is a private grid-engine wrapper that
    # exists neither in this repository nor on CPAN, so it is intentionally NOT
    # a dependency.  All three call sites already treat it as optional and print
    # an explanatory message when it is missing.
};
