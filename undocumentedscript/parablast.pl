#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile tempdir);
use File::Copy;
use File::Spec;

my $flag_man        = 0;
my $flag_help       = 0;
my $flag_forcesplit = 0;
my $flag_force      = 0;
my $param_blastexecname = 'blastall';
my $debug           = 0;
my $num_cpus        = 8;

GetOptions( 'help|?'     => \$flag_help,
	        'man'        => \$flag_man,
	        'debug'      => \$debug,
	        'numcpus=i'  => \$num_cpus,
	        'forcesplit' => \$flag_forcesplit,
	        'force'      => \$flag_force,
	        'blastname=s' => \$param_blastexecname,
	        ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man;
$flag_forcesplit = 1 if($flag_force);

# TGEW.pm provides the Sun Grid Engine support. It is not distributed with
# this repository, so it is loaded lazily: without it we simply run in the
# local (SMP) mode.
my $tgew_is_available = eval { require TGEW; 1 } ? 1 : 0;
if(!$tgew_is_available) {
    print STDERR "TGEW.pm is not available; assuming a non-SGE (SMP) environment.\n" if($debug);
    if($ENV{SGE_ROOT} || $ENV{SGE_CELL}) {
        print STDERR "WARNING: SGE_ROOT/SGE_CELL is set, so this looks like a Sun Grid Engine\n";
        print STDERR "         environment, but TGEW.pm could not be loaded:\n";
        print STDERR "         $@";
        print STDERR "         Jobs will be run locally instead of being submitted to the grid.\n";
    }
}
sub tgewIsInstalled {
    return 0 unless($tgew_is_available);
    return TGEW::is_installed();
}

# ParallelExec.pm may pull in TGEW.pm, so loading it can fail on a plain
# checkout. Load it lazily so that we can say what went wrong.
#
# NOTE: 'require' (unlike 'use') never calls import(), so nothing is imported
#       from ParallelExec here; its subs are called fully qualified below.
unless(eval { require ParallelExec; 1 }) {
    my $loaderror = $@;
    print STDERR "ERROR: cannot load ParallelExec.pm\n";
    print STDERR "       $loaderror";
    print STDERR "       (ParallelExec.pm itself requires TGEW.pm, which is not part of this repository.)\n"
        unless($tgew_is_available);
    exit 1;
}

# Quote a string so that it can safely be used as a single shell word.
sub shellquote {
    my $s = shift;
    $s = '' unless(defined $s);
    $s =~ s/'/'\\''/g;
    return "'$s'";
}

my $outputFormatType = '0';

if(tgewIsInstalled()) {
	if($num_cpus > 288) {
	    print STDERR "You specified -numcpus=$num_cpus option, but the number of CPUs cannot exceed 288 due to a technical limitation.\n";
	    exit 1;
	}
} else {
	if($num_cpus > 40) {
    	print STDERR "You specifed -numcpus=$num_cpus option, but number of CPUs cannot exceed 40 due to a limitation of zsh.\n";
    	exit 1;
    }
}

# parse BLAST options
my $blastoptions = '';
my $databasename = '';
my ($blastqueryfile, $blastoutputfile);
{
    # NOTE: 'while(my $arg = shift)' would stop at the first false argument,
    #       and BLAST options routinely take '0' as their value
    #       (e.g. '-F 0'), which silently discarded everything after it.
    my $takeNextArgument = sub {
        my $optname = shift;
        unless(@ARGV) {
            print STDERR "BLAST option '$optname' requires a value.\n";
            exit 1;
        }
        return shift @ARGV;
    };
    while(@ARGV) {
	    my $arg = shift @ARGV;
	    if($arg =~ /^-/) {
	        if ($arg =~ /^-i(.*)/) {
	            $blastqueryfile = ($1 ne '') ? $1 : $takeNextArgument->('-i');
	        } elsif ($arg =~ /^-o(.*)/) {
	            $blastoutputfile = ($1 ne '') ? $1 : $takeNextArgument->('-o');
	        } else {
	            $blastoptions .= ' ' . shellquote($arg);
	            if($arg =~ /^-m(\d*)$/) {
	                $outputFormatType = $1;
	                if($outputFormatType eq '') {
	                    $outputFormatType = $takeNextArgument->('-m');
	                    $blastoptions .= ' ' . shellquote($outputFormatType);
	                }
	            } elsif ($arg =~ /^-d(.*)/) {
	                if($1 ne '') {
	                    $databasename = $1;
	                } else {
	                    $databasename = $takeNextArgument->('-d');
	                    $blastoptions .= ' ' . shellquote($databasename);
	                }
	            }
	        }
	    } else {
	        $blastoptions .= ' ' . shellquote($arg);
	    }
    }
    if($debug) {
	    if(defined $blastoutputfile) {
	        print STDERR "Output file        = $blastoutputfile\n";
	    } else {
	        print STDERR "Output file        = STDOUT\n";
	    }
	    print STDERR "Query file         = ", (defined $blastqueryfile ? $blastqueryfile : '(not specified)'), "\n";
	    print STDERR "Format type        = $outputFormatType\n";
	    print STDERR "DATABASE           = ", ($databasename ne '' ? $databasename : '(not specified)'), "\n";
    }
}

# check if query file exists
if(!defined $blastqueryfile) {
    print STDERR "Query file is not specified. \n";
    exit 1;
}
if($databasename eq '') {
    print STDERR "No BLAST database was specified. Give '-d <database>' in the BLAST arguments.\n";
    exit 1;
}
unless(blastDatabaseExists($databasename)) {
    print STDERR "Database '$databasename' does not exist.\n";
    print STDERR "(looked for $databasename.nsq / .psq / .nal / .pal and for volume files\n";
    print STDERR " such as $databasename.00.nsq)\n";
    exit 1;
}

# A BLAST database may be a single volume ('nt.nsq'), an alias ('nt.nal' /
# 'nt.pal') or a set of volumes ('nt.00.nsq', 'nt.01.nsq', ...).
sub blastDatabaseExists {
    my $db = shift;
    for my $suffix (qw(nsq psq nal pal)) {
        return 1 if(-e "$db.$suffix");
    }
    my ($dirpath, $base);
    if($db =~ m{^(.*)[/\\]([^/\\]*)$}) {
        $dirpath = ($1 eq '') ? '/' : $1;
        $base    = $2;
    } else {
        $dirpath = '.';
        $base    = $db;
    }
    my $dh;
    return 0 unless(opendir($dh, $dirpath));
    my $found = 0;
    while(defined(my $entry = readdir($dh))) {
        if($entry =~ /^\Q$base\E\.\d+\.[np]sq$/) {
            $found = 1;
            last;
        }
    }
    closedir($dh);
    return $found;
}

# check if the specified output format is supported.
if($outputFormatType ne '0' && $outputFormatType ne '8' && $outputFormatType ne '9') {
    print STDERR "-m $outputFormatType option is not supported by paraBLAST.\n";
    print STDERR "Please contact the developers to support -m $outputFormatType if you feel it appropriate.\n";
    exit 1;
}

# check blast
# --blastname is user input, so it must not be interpolated into a shell
# command unquoted.
my $blastpath;
{
    my $whichcmd = 'which ' . shellquote($param_blastexecname);
    $blastpath = `$whichcmd`;
    $blastpath = '' unless(defined $blastpath);
    chomp $blastpath;
}
unless($blastpath ne '' && -x $blastpath) {
    print STDERR "Could not locate the path of '$param_blastexecname'\n";
    print STDERR "Make sure that BLAST is properly installed on your system and it is placed on PATH\n";
    exit 1;
} else {
    if($debug) {
		print STDERR "BLAST path         = $blastpath\n";
    }
}

# splitfasta.pl lives next to this script; fall back to PATH if it does not.
my @splitfastacmd;
{
    my $sibling = File::Spec->catfile($FindBin::Bin, 'splitfasta.pl');
    @splitfastacmd = (-e $sibling) ? ($^X, $sibling) : ('splitfasta.pl');
}

# split query
my @splitqueryfiles;
my $haveToDeleteSplitFiles = 0;
my $splitworkdir = undef; # private temporary directory, removed on success
{
    my $symlink_exists = eval { symlink("",""); 1 };
    if($symlink_exists) { # UNIX-like OS
	    # tempfile(OPEN => 0) only reserves a *name* (a race), and with no
	    # directory in the template both the name and every chunk would be
	    # created in the current working directory. Use a private
	    # directory under the system temporary directory instead.
	    $splitworkdir = tempdir('parablastXXXXXX', TMPDIR => 1, CLEANUP => 0);
	    my $splitqueryfilename = File::Spec->catfile($splitworkdir, 'query');
	    # the symlink lives in another directory, so it needs an absolute target
	    my $absolutequeryfile = File::Spec->rel2abs($blastqueryfile);
	    symlink $absolutequeryfile, $splitqueryfilename
	        or die "Could not symlink '$absolutequeryfile' to '$splitqueryfilename': $!";
	    my @cmdline = (@splitfastacmd, '--equalbase', $splitqueryfilename, $num_cpus);
	    print STDERR "  % ", join(' ', map { shellquote($_) } @cmdline), "\n";
	    if(system(@cmdline)) {
	        die "Could not execute splitfasta.pl";
	    }
	    unlink $splitqueryfilename; #remove symbolic link
	    for(0..$num_cpus-1) {
	        my $p = $_ + 1;
	        my $splitpartfastafilename = "${splitqueryfilename}.$p";
	        if(-e $splitpartfastafilename) {
		        push(@splitqueryfiles, $splitpartfastafilename);
	        } else {
	            if($_ == 0) {
	                print STDERR "ERROR: splitfasta.pl did not produce any split query file.\n";
	                print STDERR "       Check that '$blastqueryfile' is a non-empty FASTA file.\n";
	                exit 1;
	            }
	            # for small FASTA
	            print STDERR "The given FASTA file has small number of reads.\n";
	            print STDERR "You gave $num_cpus CPUs, but the number of blocks is $_\n";
	            print STDERR "The number of CPUs decreased to $_, and then proceeds.\n";
	            print STDERR "This is not to warn you but just to inform you about the processing.\n";
	            $num_cpus = $_;
	            last;
	        }
	    }
	    $haveToDeleteSplitFiles = 1;
    } else { # OS without symbolic link (i.e. Windows)
	    if(-e "${blastqueryfile}.1" && !$flag_forcesplit) {
	        print STDERR "There seems to be split query file already.\n";
	        print STDERR "Abort for safety.\n";
	        print STDERR "If you intend to overwrite these split query files, add --forcesplit option.\n";
	        exit 1;
	    } else {
	        $haveToDeleteSplitFiles = 1;
	    }
	    my @cmdline = (@splitfastacmd, '--equalbase', $blastqueryfile, $num_cpus);
	    print STDERR "  % ", join(' ', map { shellquote($_) } @cmdline), "\n";
	    if(system(@cmdline)) {
	        die "Could not execute splitfasta.pl";
	    }
	    for(0..$num_cpus-1) {
	        my $p = $_ + 1;
	        push(@splitqueryfiles, "$blastqueryfile.$p");
	    }
    }
    while(@splitqueryfiles > 0 && !-e $splitqueryfiles[-1]) {
	    pop(@splitqueryfiles);
	    print STDERR "Decreased the number of parallelism\n";
	    $num_cpus--;
    }
    # With no chunk at all the combine step below would run 'cat' with no
    # operand, which reads standard input and hangs forever.
    if(@splitqueryfiles == 0) {
        print STDERR "ERROR: no split query file was produced from '$blastqueryfile'.\n";
        print STDERR "       Nothing to do; aborting.\n";
        rmdir $splitworkdir if(defined $splitworkdir);
        exit 1;
    }
}

my $errorHasOccured = undef;

# execute BLAST in parallel
EXECUTEBLASTINPARALLEL: {
    my @commandlines;
    for(0..$num_cpus-1) {
	    my $outputfile = "$splitqueryfiles[$_].out";
	    if(-e $outputfile && !$flag_force) {
	        print STDERR "'$outputfile' already exists! Abort for safety.\n";
	        print STDERR "Try -force option if you want to overwrite them\n";
	        exit 1;
	    }
	    # ParallelExec writes these strings into a shell script / makefile,
	    # so every word has to be quoted. '$(O)' is a marker that
	    # ParallelExec strips, so it must stay outside the quotes.
	    push(@commandlines,
	         shellquote($blastpath) .
	         ' -i ' . shellquote($splitqueryfiles[$_]) .
	         ' -o ' . "\$(O)" . shellquote($outputfile) .
	         $blastoptions);
    }
    if($debug) {
	    print STDERR "Command lines : \n";
	    for(@commandlines) {
	        print STDERR "  $_\n";
	    }
	}
    my %retobj = ParallelExec::parallelExecute(@commandlines);
    if($retobj{error}) {
	    $errorHasOccured = "Parallel execution of BLAST failed";
	    last;
    }
    for(0..$num_cpus-1) {
	    if(${$retobj{errorlevels}}[$_]) {
	        $errorHasOccured = "Job $_ (execution of BLAST) failed";
	        last EXECUTEBLASTINPARALLEL;
        }
    }
}

# clean divided query
# When a job failed the chunk files are kept so that the run can be inspected
# and/or restarted; deleting them would throw away the partial results.
if($haveToDeleteSplitFiles) {
    unless($debug || $errorHasOccured) {
	    for(@splitqueryfiles) {
	        print STDERR "  rm $_\n";
	        unlink $_;
	    }
    } else {
	    print STDERR "  Do not delete temporary files shown below. \n";
	    for(@splitqueryfiles) {
	        print STDERR "    $_\n";
	    }
    }
}

# put the results together
unless($errorHasOccured) {
    my @splitresultfiles = map {"$_.out"} @splitqueryfiles;
    if($outputFormatType eq '8' || $outputFormatType eq '9') {
	    print STDERR "Combine strategy : blast -m $outputFormatType\n";
	    my $cmdline = "cat " . join(' ', map { shellquote($_) } @splitresultfiles);
	    $cmdline .= " > " . shellquote($blastoutputfile) if(defined $blastoutputfile && $blastoutputfile ne '');
	    print STDERR "Put all things together\n";
	    print STDERR "% $cmdline\n";
	    if(system($cmdline)) {
	        $errorHasOccured = "Error while executing cat command";
	    }
    } elsif($outputFormatType eq '0') {
	    print STDERR "Combine strategy : blast -m 0 (pairwise)\n";
	    # A plain 'cat' would repeat the program/reference header once per
	    # chunk. Keep the header of the first chunk only. The per-chunk
	    # 'Database:' / statistics trailers are kept as they are: they hold
	    # each chunk's own statistics and cannot be merged meaningfully.
	    print STDERR "NOTE: the combined -m 0 report keeps one 'Database:'/statistics block\n";
	    print STDERR "      per chunk; use -m 8 or -m 9 if you need a single tabular file.\n";
	    my $destination = (defined $blastoutputfile && $blastoutputfile ne '') ? $blastoutputfile : undef;
	    unless(defined $destination) {
	        # no -o was given, so BLAST wrote to files and we dump to STDOUT
	        for my $index (0..@splitresultfiles-1) {
	            printFileFromFirstBodyLine(sub { /^Query=/ }, $splitresultfiles[$index], $index == 0);
	        }
	    } else {
	        my $firstone = $splitresultfiles[0];
	        unless(copy($firstone, $destination)) {
	            $errorHasOccured = "Could not copy '$firstone' to '$destination': $!";
	        } else {
	            for(1..@splitresultfiles-1) {
	                appendFileFromFirstBodyLine(sub { /^Query=/ }, $splitresultfiles[$_], $destination);
	            }
	        }
	    }
    } else {
        print STDERR "Unknown output type '$outputFormatType'\n";
        $errorHasOccured = "Unknown output type '$outputFormatType'";
    }
}

# clean devided results
{
    unless($debug || $errorHasOccured) {
	    for(@splitqueryfiles) {
	        my $outputfile = "$_.out";
	        print STDERR "  rm $outputfile\n";
	        unlink $outputfile;
	    }
	    rmdir $splitworkdir if(defined $splitworkdir);
    } else {
	    print STDERR "  Do not delete temporary files shown below. \n";
	    for(@splitqueryfiles) {
	        print STDERR "    $_.out\n";
	    }
    }
}

if($errorHasOccured) {
    print STDERR "ERROR: $errorHasOccured\n";
    print STDERR "ERROR: parablast.pl did not complete successfully.\n";
    print STDERR "       The per-chunk files listed above were kept for inspection.\n";
    exit 1;
}

# Append $sourcefile to $destfile, dropping everything before the first line
# for which $isbody_function is true (the repeated BLAST program header).
sub appendFileFromFirstBodyLine {
    my $isbody_function = shift;
    my $sourcefile = shift;
    my $destfile   = shift;
    open(APPENDFILE, '<', $sourcefile) or do {
	    print STDERR "Cannot read '$sourcefile'\n";
	    return;
    };
    unless(open(APPENDEDFILE, '>>', $destfile)) {
	    print STDERR "Cannot append to '$destfile'\n";
	    close APPENDFILE;
	    return;
    }
    my $isinheader = 1;
    while(<APPENDFILE>) {
	    if($isinheader) {
	        next unless(&{$isbody_function});
	        $isinheader = 0;
	    }
	    print APPENDEDFILE;
    }
    close APPENDEDFILE;
    close APPENDFILE;
    return 0;
}

# Same as above but writes to STDOUT; $keepheader selects the first chunk,
# whose header has to be preserved.
sub printFileFromFirstBodyLine {
    my $isbody_function = shift;
    my $sourcefile = shift;
    my $keepheader  = shift;
    open(APPENDFILE, '<', $sourcefile) or do {
	    print STDERR "Cannot read '$sourcefile'\n";
	    return;
    };
    my $isinheader = $keepheader ? 0 : 1;
    while(<APPENDFILE>) {
	    if($isinheader) {
	        next unless(&{$isbody_function});
	        $isinheader = 0;
	    }
	    print;
    }
    close APPENDFILE;
    return 0;
}

=pod

=head1 NAME

parablast - Parallel BLAST

=head1 SYNOPSIS

parablast [options for parablast] -- [arguments for BLAST]

Options:
   -numcpus=nn          specify the number of parallelism
   -help            brief help message
   -man             full documentation
   -forcesplit      overwrite when splitting query file. (On windows)

=head1 OPTIONS

=over 8

=item B<-numcpus>

Specify the number of parallelism. If you give -numcpus=10, Query sequences are separated into 10, and 10 BLAST instances are executed. Due to a limitation of zsh, the number of parallelism cannot exceed 40 in the current implementation.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-forcesplit>

Overwrite files when splitting query file.  On UNIX systems, it first symbolic links to temporary, then it splits temporary files into many. However, on Windows system, there is no mechanism of symbolic link. Split files are named "originalname.x" where x is index number of splitting.  These files are overwritten if B<-forcesplit> option is specified.

=back

=head1 DESCRIPTION

B<This program> will divide queries into as many parts as specifed, execute BLAST in parallell,
put the results together into one file.

=cut
