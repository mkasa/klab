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
use FileHandle;

my $flag_man        = 0;
my $flag_help       = 0;
my $flag_forcesplit = 0;
my $flag_force      = 0;
my $debug           = 0;
my $num_cpus        = 8;

GetOptions( 'help|?'     => \$flag_help,
	        'man'        => \$flag_man,
	        'debug'      => \$debug,
	        'numcpus=i'  => \$num_cpus,
	        'forcesplit' => \$flag_forcesplit,
	        'force'      => \$flag_force
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

# ParallelExec.pm currently pulls in TGEW.pm unconditionally, so loading it may
# fail on a plain checkout. Load it lazily so that we can say what went wrong.
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

my $outputFormatType = 'psl';

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

# parse BLAT options
my $blatoptions = '';
my ($blatsubjectfile, $blatqueryfile, $blatoutputfile);
{
    my @blatargs;
    # NOTE: 'while(my $arg = shift)' would stop at a false argument such as a
    #       file literally named '0'.
    while(@ARGV) {
	    my $arg = shift @ARGV;
	    if($arg =~ /^-/) {
	        $blatoptions .= ' ' . shellquote($arg);
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
my $blatpath = `which blat`;
$blatpath = '' unless(defined $blatpath);
chomp $blatpath;
unless($blatpath ne '' && -x $blatpath) {
    print STDERR "Could not locate the path of blat\n";
    print STDERR "Make sure that blat is properly installed on your system and it is placed on PATH\n";
    exit 1;
} else {
    if($debug) {
	    print STDERR "BLAT path         = $blatpath\n";
    }
}

# splitfasta.pl lives next to this script; fall back to PATH if it does not.
# Both the UNIX and the Windows branch below must invoke it the same way.
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
	    $splitworkdir = tempdir('parablatXXXXXX', TMPDIR => 1, CLEANUP => 0);
	    my $splitqueryfilename = File::Spec->catfile($splitworkdir, 'query');
	    # the symlink lives in another directory, so it needs an absolute target
	    my $absolutequeryfile = File::Spec->rel2abs($blatqueryfile);
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
	                print STDERR "       Check that '$blatqueryfile' is a non-empty FASTA file.\n";
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
	    if(-e "${blatqueryfile}.1" && !$flag_forcesplit) {
	        print STDERR "There seems to be split query file already.\n";
	        print STDERR "Abort for safety.\n";
	        print STDERR "If you intend to overwrite these split query files, add --forcesplit option.\n";
	        exit 1;
	    } else {
	        $haveToDeleteSplitFiles = 1;
	    }
	    my @cmdline = (@splitfastacmd, '--equalbase', $blatqueryfile, $num_cpus);
	    print STDERR "  % ", join(' ', map { shellquote($_) } @cmdline), "\n";
	    if(system(@cmdline)) {
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
    # With no chunk at all the combine step below would run 'cat' with no
    # operand, which reads standard input and hangs forever.
    if(@splitqueryfiles == 0) {
        print STDERR "ERROR: no split query file was produced from '$blatqueryfile'.\n";
        print STDERR "       Nothing to do; aborting.\n";
        rmdir $splitworkdir if(defined $splitworkdir);
        exit 1;
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
	    # ParallelExec writes these strings into a shell script / makefile,
	    # so every word has to be quoted. '$(O)' is a marker that
	    # ParallelExec strips, so it must stay outside the quotes.
	    push(@commandlines,
	         shellquote($blatpath) . ' ' .
	         shellquote($blatsubjectfile) . ' ' .
	         shellquote($splitqueryfiles[$_]) . ' ' .
	         "\$(O)" . shellquote($outputfile) .
	         $blatoptions);
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
    if($outputFormatType eq 'blast8' || $outputFormatType eq 'blast9' || $outputFormatType eq 'wublast' || $outputFormatType eq 'sim4') {
        print STDERR "Combining strategy : $outputFormatType\n";
	    # just concatinating all
	    my $cmdline = "cat " . join(' ', map { shellquote($_) } @splitresultfiles)
	                . " > " . shellquote($blatoutputfile);
	    print STDERR "Put all things together\n";
	    print STDERR "% $cmdline\n";
	    if(system($cmdline)) {
	        $errorHasOccured = "Error while executing cat command";
	    }
    } elsif($outputFormatType eq 'maf') {
        print STDERR "Combining strategy : maf\n";
	    # A plain 'cat' would repeat the '##maf version=' header once per
	    # chunk in the middle of the file, which most MAF parsers reject.
	    my $firstone = $splitresultfiles[0];
	    unless(copy($firstone, $blatoutputfile)) {
	        $errorHasOccured = "Could not copy '$firstone' to '$blatoutputfile': $!";
	    } else {
	        for(1..@splitresultfiles-1) {
	            appendMafWithoutHeader($splitresultfiles[$_], $blatoutputfile);
	        }
	    }
    } elsif($outputFormatType eq 'blast') {
		print STDERR "Combining strategy : blast\n";
	    my $firstone = $splitresultfiles[0];
	    unless(copy($firstone, $blatoutputfile)) {
	        $errorHasOccured = "Could not copy '$firstone' to '$blatoutputfile': $!";
	    } else {
	        for(1..@splitresultfiles-1) {
	            appendFileWithHeaderSkip(sub { /^Reference/ }, $splitresultfiles[$_], $blatoutputfile);
	        }
	    }
    } elsif($outputFormatType eq 'psl' || $outputFormatType eq 'pslx') {
		print STDERR "Combining strategy : psl/pslx\n";
	    my $firstone = $splitresultfiles[0];
	    unless(copy($firstone, $blatoutputfile)) {
	        $errorHasOccured = "Could not copy '$firstone' to '$blatoutputfile': $!";
	    } else {
	        for(1..@splitresultfiles-1) {
	            appendFileWithHeaderSkip(sub { /^---------/ }, $splitresultfiles[$_], $blatoutputfile);
	        }
	    }
    } elsif($outputFormatType eq 'axt') {
		print STDERR "Combining strategy : axt\n";
		# appendAxt() opens the destination with '>>', so it has to be
		# truncated first; otherwise re-running would double the output.
		my $truncated = 0;
		if(open(my $truncatefh, '>', $blatoutputfile)) {
		    close($truncatefh);
		    $truncated = 1;
		}
		unless($truncated) {
		    $errorHasOccured = "Could not create '$blatoutputfile': $!";
		} else {
		    my $f = { currentid => 0 };
		    for(0..@splitresultfiles-1) {
	    	    appendAxt($f, $splitresultfiles[$_], $blatoutputfile);
		    }
		}
    } else {
        print STDERR "Unknown output type '$outputFormatType'\n";
        $errorHasOccured = "Unknown output type '$outputFormatType'";
    }
}

sub appendAxt
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
			next if(/^\s*$/); # blank line between two blocks
			if(/^(\d+)\s+(.*)$/) {
				my $rest  = $2;
				my $newid = $idobj->{currentid}++;
				my $s1 = <$fh>;
				my $s2 = <$fh>;
				unless(defined $s1 && defined $s2) {
					print STDERR "WARNING: truncated axt block at the end of '$axtfile'\n";
					last;
				}
				my $blank = <$fh>; # the (optional) separating blank line
				print $outfh "$newid $rest\n";
				print $outfh $s1;
				print $outfh $s2;
				print $outfh "\n";
			} else {
				print STDERR "WARNING: illegal header in '$axtfile'\n";
				print STDERR $_;
			}
		}
		$fh->close();
		$outfh->close();
	} else {
		print STDERR "ERROR: could not open $axtfile\n";
		$outfh->close() if($outfh);
		exit 1;
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
    print STDERR "ERROR: parablat.pl did not complete successfully.\n";
    print STDERR "       The per-chunk files listed above were kept for inspection.\n";
    exit 1;
}

sub appendFileWithHeaderSkip {
    my $isheader_function = shift;
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

# Append a MAF file, dropping its leading '##maf'/'#' header block so that the
# combined file has exactly one header.
sub appendMafWithoutHeader {
    my $sourcefile = shift;
    my $destfile   = shift;
    open(my $in, '<', $sourcefile) or do {
	    print STDERR "Cannot read '$sourcefile'\n";
	    return;
    };
    my $out;
    unless(open($out, '>>', $destfile)) {
	    print STDERR "Cannot append to '$destfile'\n";
	    close $in;
	    return;
    }
    my $isinheader = 1;
    while(<$in>) {
	    if($isinheader) {
	        next if(/^#/ || /^\s*$/);
	        $isinheader = 0;
	    }
	    print $out $_;
    }
    close $out;
    close $in;
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
