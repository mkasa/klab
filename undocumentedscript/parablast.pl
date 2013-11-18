#!/usr/bin/env perl

use strict;

use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile tempdir);
use ParallelExec;
use File::Copy;
use TGEW;

my $flag_man        = 0;
my $flag_help       = 0;
my $flag_forcesplit = 0;
my $flag_force      = 0;
my $param_blastexecname = 'blastall';
my $debug           = 0;
my $num_cpus        = 8;

GetOptions( 'help|?'     => \$flag_help,
	        'debug'      => \$debug,
	        'numcpus=i'  => \$num_cpus,
	        'forcesplit' => \$flag_forcesplit,
	        'force'      => \$flag_force,
	        'blastname=s' => \$param_blastexecname,
	        ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man;
$flag_forcesplit = 1 if($flag_force);

my $outputFormatType = '0';

if(TGEW::is_installed()) {
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
    while(my $arg = shift) {
	    if($arg =~ /^-/) {
	        if ($arg =~ /^-i(.*)/) {
	            if(defined $1 && $1 ne '') {
	                $blastqueryfile = $1;
	            } else {
	                $blastqueryfile = shift;
	            }
	        } elsif ($arg =~ /^-o(.*)/) {
	            if(defined $1 && $1 ne '') {
	                $blastoutputfile = $1;
	            } else {
	                $blastoutputfile = shift;
	            }
	        } else {
	            $blastoptions .= " $arg";
	            if($arg =~ /^-m(\d*)/) {
	                $outputFormatType = $1;
	                unless(defined $outputFormatType && $outputFormatType ne '') {
	                    $outputFormatType = shift;
	                    $blastoptions .= " $outputFormatType";
	                }
	            } elsif ($arg =~ /^-d(.*)/) {
	                if(defined $1 && $1 ne '') {
	                    $databasename = $1;
	                } else {
	                    $databasename = shift;
	                    $blastoptions .= " $databasename";
	                }
	            }
	        }
	    } else {
	        $blastoptions .= " $arg";
	    }
    }
    if($debug) {
	    if(defined $blastoutputfile) {
	        print STDERR "Output file        = $blastoutputfile\n";
	    } else {
	        print STDERR "Output file        = STDOUT\n";
	    }
	    print STDERR "Query file         = $blastqueryfile\n";
	    print STDERR "Format type        = $outputFormatType\n";
	    print STDERR "DATABASE           = $databasename\n";
    }
}

# check if query file exists
if(!defined $blastqueryfile) {
    print STDERR "Query file is not specified. \n";
    exit 1;
}
if(!-e "$databasename.nsq" && !-e "$databasename.psq") {
    print STDERR "Database file does not exist.\n";
    exit 1;
}

# check if the specified output format is supported.
if($outputFormatType ne '0' && $outputFormatType ne '8' && $outputFormatType ne '9') {
    print STDERR "-m $outputFormatType option is not supported by paraBLAST.\n";
    print STDERR "Please contact the developers to support -m $outputFormatType if you feel it appropriate.\n";
    exit 1;
}

# check blast
my $blastpath = `which $param_blastexecname`; chomp $blastpath;
unless(-x $blastpath) {
    print STDERR "Could not locate the path of '$param_blastexecname'\n";
    print STDERR "Make sure that BLAST is properly installed on your system and it is placed on PATH\n";
    exit 1;
} else {
    if($debug) {
		print STDERR "BLAST path         = $blastpath\n";
    }
}

# split query
my @splitqueryfiles;
my $haveToDeleteSplitFiles = 0;
{
    my $symlink_exists = eval { symlink("",""); 1 };
    if($symlink_exists) { # UNIX-like OS
	    my $splitqueryfilenamebase = 'parablastqueryXXXXX';
	    my (undef, $splitqueryfilename) = tempfile($splitqueryfilenamebase, OPEN => 0);
	    symlink $blastqueryfile, $splitqueryfilename or die "Could not simlink to temporary";
	    my $cmdline = "splitfasta.pl --equalbase $splitqueryfilename $num_cpus";
	    print STDERR "  %$cmdline\n";
	    if(system $cmdline) {
	        die "Could not execute splitfasta.pl";
	    }
	    unlink $splitqueryfilename; #remove symbolic link
	    for(0..$num_cpus-1) {
	        my $p = $_ + 1;
	        my $splitpartfastafilename = "${splitqueryfilename}.$p";
	        if(-e $splitpartfastafilename) {
		        push(@splitqueryfiles, $splitpartfastafilename);
	        } else {
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
	    my $cmdline = "perl splitfasta.pl --equalbase $blastqueryfile $num_cpus";
	    print STDERR "  %$cmdline\n";
	    if(system $cmdline) {
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
	    push(@commandlines, "$blastpath -i $splitqueryfiles[$_] -o \$(O)$outputfile $blastoptions");
    }
    if($debug) {
	    print STDERR "Command lines : \n";
	    for(@commandlines) {
	        print STDERR "  $_\n";
	    }
	}
    my %retobj = parallelExecute(@commandlines);
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
if($haveToDeleteSplitFiles) {
    unless($debug) {
	    for(@splitqueryfiles) {
	        print STDERR "  rm $_\n";
	        unlink $_;
	    }
    } else {
	    print STDERR "  Debugging mode. Do not delete temporary files shown below. \n";
	    for(@splitqueryfiles) {
	        print STDERR "    $_\n";
	    }
    }
}

# put the results together
unless($errorHasOccured) {
    my @splitresultfiles = map {"$_.out"} @splitqueryfiles;
    if($outputFormatType eq '0' || $outputFormatType eq '8' || $outputFormatType eq '9') {
	    print STDERR "Combine strategy : blast -m $outputFormatType\n";
	    my $cmdline;
	    $cmdline = "cat " . join(' ',  @splitresultfiles);
	    $cmdline .= " > $blastoutputfile" if(defined $blastoutputfile && $blastoutputfile ne '');
	    print STDERR "Put all things together\n";
	    print STDERR "% $cmdline\n";
	    system $cmdline;
	    if($?) {
	        print STDERR "Error while executing cat command\n";
	    }
    } else {
        print STDERR "Unknown output type\n";
    }
}

# clean devided results
{
    unless($debug) {
	    for(@splitqueryfiles) {
	        my $outputfile = "$_.out";
	        print STDERR "  rm $outputfile\n";
	        unlink $outputfile;
	    }
    } else {
	    print STDERR "  Debugging mode. Do not delete temporary files shown below. \n";
	    for(@splitqueryfiles) {
	        print STDERR "    $_.out\n";
	    }
    }
}

sub appendFileWithHeaderSkip($$$) {
    my $isheader_function = shift;
    my $sourcefile = shift;
    my $destfile   = shift;
    open APPENDFILE, "< $sourcefile" or return;
    unless(open APPENDEDFILE, ">> $destfile") {
	    print STDERR "Cannot append to '$destfile'\n";
	    close APPENDFILE;
	    return;
    }
    while(<APPENDFILE>) {
		last if(&{$isheader_function});
    }
    while(<APPENDFILE>) {
		print APPENDEDFILE;
    }
    close APPENDEDFILE;
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
