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
my $debug           = 0;
my $num_cpus        = 8;

GetOptions( 'help|?'     => \$flag_help,
	        'debug'      => \$debug,
	        'numcpus=i'  => \$num_cpus,
	        'forcesplit' => \$flag_forcesplit,
	        'force'      => \$flag_force
	        ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man;
$flag_forcesplit = 1 if($flag_force);

my $outputFormatType = 'psl';

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

# parse BLAT options
my $blatoptions = '';
my ($blatsubjectfile, $blatqueryfile, $blatoutputfile);
{
    my @blatargs;
    while(my $arg = shift) {
	    if($arg =~ /^-/) {
	        $blatoptions .= " $arg";
	        if($arg =~ /^-out=([\w\d]+)/) {
	            $outputFormatType = $1;
	        }
	    } else {
	        push(@blatargs, $arg);
	    }
    }
    if(@blatargs != 3) {
	    print STDERR "BLAT requires 3 arguments.\n";
	    print STDERR "  usage: blat database query output\n";
	    print STDERR "  usage: parablat [parablat options] -- [blat arguments]\n";
	    print STDERR "just type blat to see arguments of BLAT\n";
	    print STDERR "parablat.pl -man to see manual\n";
	    exit 1;
    }
    ($blatsubjectfile, $blatqueryfile, $blatoutputfile) = @blatargs;
    if($debug) {
	    print STDERR "BLAT subject file = $blatsubjectfile\n";
	    print STDERR "BLAT query   file = $blatqueryfile\n";
	    print STDERR "BLAT output  file = $blatoutputfile\n";
    }
}

# check blat
my $blatpath = `which blat`; chomp $blatpath;
unless(-x $blatpath) {
    print STDERR "Could not locate the path of blat\n";
    print STDERR "Make sure that blat is properly installed on your system and it is placed on PATH\n";
    exit 1;
} else {
    if($debug) {
	    print STDERR "BLAT path         = $blatpath\n";
    }
}

# split query
my @splitqueryfiles;
my $haveToDeleteSplitFiles = 0;
{
    my $symlink_exists = eval { symlink("",""); 1 };
    if($symlink_exists) { # UNIX-like OS
	    my $splitqueryfilenamebase = 'parablatqueryXXXXX';
	    my (undef, $splitqueryfilename) = tempfile($splitqueryfilenamebase, OPEN => 0);
	    symlink $blatqueryfile, $splitqueryfilename or die "Could not simlink to temporary";
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
	    if(-e "${blatqueryfile}.1" && !$flag_forcesplit) {
	        print STDERR "There seems to be split query file already.\n";
	        print STDERR "Abort for safety.\n";
	        print STDERR "If you intend to overwrite these split query files, add --forcesplit option.\n";
	        exit 1;
	    } else {
	        $haveToDeleteSplitFiles = 1;
	    }
	    my $cmdline = "perl splitfasta.pl --equalbase $blatqueryfile $num_cpus";
	    print STDERR "  %$cmdline\n";
	    if(system $cmdline) {
	        die "Could not execute splitfasta.pl";
	    }
	    for(0..$num_cpus-1) {
	        my $p = $_ + 1;
	        push(@splitqueryfiles, "$blatqueryfile.$p");
	    }
    }
    while(@splitqueryfiles > 0 && !-e $splitqueryfiles[-1]) {
	    pop(@splitqueryfiles);
	    print STDERR "Decreased the number of parallelism\n";
	    $num_cpus--;
    }
}

my $errorHasOccured = undef;

# execute BLAT in parallel
EXECUTEBLATINPARALLEL: {
    my @commandlines;
    for(0..$num_cpus-1) {
	    my $outputfile = "$splitqueryfiles[$_].out";
	    if(-e $outputfile && !$flag_force) {
	        print STDERR "'$outputfile' already exists! Abort for safety.\n";
	        print STDERR "Try -force option if you want to overwrite them\n";
	        exit 1;
	    }
	    push(@commandlines, "$blatpath $blatsubjectfile $splitqueryfiles[$_] \$(O)$outputfile $blatoptions");
    }
    if($debug) {
	    print "Command lines : \n";
	    for(@commandlines) {
	        print "  $_\n";
	    }
	}
    my %retobj = parallelExecute(@commandlines);
    if($retobj{error}) {
	    $errorHasOccured = "Parallel execution of BLAT failed";
	    last;
    }
    for(0..$num_cpus-1) {
	    if(${$retobj{errorlevels}}[$_]) {
	        $errorHasOccured = "Job $_ (execution of BLAT) failed";
	        last EXECUTEBLATINPARALLEL;
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
    if($outputFormatType eq 'blast8' || $outputFormatType eq 'blast9' || $outputFormatType eq 'maf' || $outputFormatType eq 'wublast' || $outputFormatType eq 'sim4') {
        print STDERR "Combining strategy : $outputFormatType\n";
	    # just concatinating all
	    my $cmdline = "cat " . join(' ',  @splitresultfiles) . " > $blatoutputfile";
	    print STDERR "Put all things together\n";
	    print STDERR "% $cmdline\n";
	    system $cmdline;
	    if($?) {
	        print STDERR "Error while executing cat command\n";
	    }
    } elsif($outputFormatType eq 'blast') {
		print STDERR "Combining strategy : blast\n";
	    my $firstone = $splitresultfiles[0];
	    copy($firstone, $blatoutputfile);
	    for(1..@splitresultfiles-1) {
	        appendFileWithHeaderSkip(sub { /^Reference/ }, $splitresultfiles[$_], $blatoutputfile);
	    }
    } elsif($outputFormatType eq 'psl' || $outputFormatType eq 'pslx') {
		print STDERR "Combining strategy : psl/pslx\n";
	    my $firstone = $splitresultfiles[0];
	    copy($firstone, $blatoutputfile);
	    for(1..@splitresultfiles-1) {
	        appendFileWithHeaderSkip(sub { /^---------/ }, $splitresultfiles[$_], $blatoutputfile);
	    }
    } elsif($outputFormatType eq 'axt') {
		print STDERR "Combining strategy : axt\n";
		my $f = { nextid => 0 };
	    for(0..@splitresultfiles-1) {
	    	appendAxt($f, $splitresultfiles[$_], $blatoutputfile);
	    }
    } else {
        print STDERR "Unknown output type '$outputFormatType'\n";
    }
}

sub appendAxt($$$)
{
	my $idobj          = shift;
	my $axtfile        = shift;
	my $parablatoutput = shift;
	my $fh = new FileHandle "< $axtfile";
	my $outfh = new FileHandle ">> $parablatoutput";
	if($fh) {
		unless($outfh) {
			$fh->close();
			print STDERR "ERROR: could not open $parablatoutput\n";
			exit 1;
		}
		while(<$fh>) {
			if(/^(\d+)\s+(.*)$/) {
				my $idnum = $1;
				my $rest  = $2;
				my $newid = $idobj->{currentid}++;
				print $outfh "$newid $rest\n";
				my $s1 = <$fh>;
				print $outfh $s1;
				my $s2 = <$fh>;
				print $outfh $s2;
				my $blank = <$fh>;
				print $outfh "\n";
			} else {
				print STDERR "WARNING: illegal header\n";
				print STDERR;
			}
		}
		$fh->close();
		$outfh->close();
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

parablat - Parallel BLAT

=head1 SYNOPSIS

parablat [options for parablat] -- [arguments for BLAT]

Options:
   -numcpus=nn          specify the number of parallelism
   -help            brief help message
   -man             full documentation
   -forcesplit      overwrite when splitting query file. (On windows)

=head1 OPTIONS

=over 8

=item B<-numcpus>

Specify the number of parallelism. If you give -numcpus=10, Query sequences are separated into 10, and 10 BLAT instances are executed.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-forcesplit>

Overwrite files when splitting query file.  On UNIX systems, it first symbolic links to temporary, then it splits temporary files into many. However, on Windows system, there is no mechanism of symbolic link. Split files are named "originalname.x" where x is index number of splitting.  These files are overwritten if B<-forcesplit> option is specified.

=back

=head1 DESCRIPTION

B<This program> will divide queries into as many parts as specifed, execute BLAT in parallell,
put the results together into one file.

=cut
