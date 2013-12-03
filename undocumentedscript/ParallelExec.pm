use strict;

use File::Temp qw(tempfile tempdir);
use FileHandle;
use TGEW;

sub parallelExecute(@) {
	if(TGEW::is_installed()) {
		parallelExecute_SGE(@_);
	} else {
		parallelExecute_SMP(@_);
	}
}

# parallel execution of given commands (for Sun Grid Engine environment)
# returns hash.
# $ret = parallelExecute(('command1', 'command2', 'command3'));
sub parallelExecute_SGE(@)
{
    my @commands = @_;
    my %rethash = ();

    my $makefiletemplate = 'tempXXXXX';
    my (undef, $makefile) = tempfile($makefiletemplate, OPEN => 0);

    my $stdouttemplate = 'tempsXXXXX';
    my @stdoutfile;

    for(@commands) {
	    my (undef, $sofile) = tempfile($stdouttemplate, OPEN => 0);
	    push(@stdoutfile, $sofile);
    }

    {
	    my $fh = new FileHandle "> $makefile";
	    $rethash{error} = 1;
	    return unless(defined $fh);
	    print $fh "all : ";
	    print $fh join(' ', @stdoutfile);
	    print $fh "\n\n";
	    for(my $i = 0; $i < @commands; $i++) {
	        print $fh "$stdoutfile[$i] :\n";
	        print $fh "\t$commands[$i] > $stdoutfile[$i]\n\n";
	    }
	    $fh->close();
    }

    my $tgepath = `which tge_make`; chomp $tgepath;
    my $commandline = "$tgepath -f $makefile -tgelock";
    if(0) {
	    print STDERR "\% $commandline\n";
	} else {
		system "$commandline -n";
	}
    $rethash{shellstdout} = `$commandline`;
    $rethash{shellerrorlevel} = $?;

    # print STDERR "ERRORLEVEL = $rethash{shellerrorlevel}\n";

    SETRETVAL: {
	    if($?) {
			print STDERR "ERROR CODE: $?\n";
			last;
		}
	    my @errorlevels;
	    for(@commands) { push(@errorlevels, $rethash{shellerrorlevel}); } # tenuki
	    $rethash{errorlevels} = \@errorlevels;

	    my @outs;
	    for(my $i = 0; $i < @commands; $i++) {
	        my $pfd2 = new FileHandle "< $stdoutfile[$i]";
	        # print STDERR "OP:$stdoutfile[$i]\n";
	        unless(defined $pfd2) {
				push(@outs, []);
			} else {
	            print STDERR "FLOP:$stdoutfile[$i]\n";
	        	my @outputstrings = <$pfd2>;
	        	push(@outs, \@outputstrings);
	        	$pfd2->close();
			}
	    }
	    $rethash{stdouts} = \@outs;

	    $rethash{error} = 0;
        # print STDERR "NOERR\n";
    }
    unlink $makefile;
    for(my $i = 0; $i < @commands; $i++) {
	    unlink $stdoutfile[$i];
    }
    return %rethash;
}

sub removeIOTags($)
{
	my $cmdstr = shift;
	$cmdstr =~ s/\$\([IO]\)//gi;
	return $cmdstr;
}

# parallel execution of given commands (for SMP/zsh environment)
# returns hash.
# $ret = parallelExecute_SMP(('command1', 'command2', 'command3'));
sub parallelExecute_SMP(@)
{
    my @commands = @_;
    my %rethash = ();

    my $zshtemplate = 'tempXXXXX';
    my (undef, $zshfile) = tempfile($zshtemplate, OPEN => 0);

    my $zshpath = `which zsh`; chomp $zshpath;
    unless(-x $zshpath) {
	    print STDERR "This script requires 'zsh', which is not found on system.\n";
	    print STDERR "Please check if zsh is properly installed.\n";
	    print STDERR "zsh must be on PATH environmental variable.\n";
	    exit 1;
    }

    my $stdouttemplate = 'tempsXXXXX';
    my @stdoutfile;
    my $statustemplate = 'tempeXXXXX';
    my @statusfile;

    for(@commands) {
	    my (undef, $sofile) = tempfile($stdouttemplate, OPEN => 0);
	    push(@stdoutfile, $sofile);
	    my (undef, $stfile) = tempfile($statustemplate, OPEN => 0);
	    push(@statusfile, $stfile);
    }

    open PARALLELEXEFH, ">$zshfile";
    print PARALLELEXEFH "#!$zshpath\n";
    for(my $i = 0; $i < @commands; $i++) {
	    my $cmd  = removeIOTags($commands[$i]);
	    my $sout = $stdoutfile[$i];
	    my $stat = $statusfile[$i];
	    print PARALLELEXEFH "($cmd; print \$? > $stat) > $sout &\n";
    }
    print PARALLELEXEFH "wait\n";
    close PARALLELEXEFH;

    $rethash{error} = 1;
    $rethash{shellstdout} = `zsh $zshfile`;
    $rethash{shellerrorlevel} = $?;
    SETRETVAL: {
	last if($?);
	my @errorlevels;
	for(my $i = 0; $i < @commands; $i++) {
	    open PARALLELEXEFH2, "<$statusfile[$i]" or last SETRETVAL;
	    my $line = <PARALLELEXEFH2>; chomp $line;
	    push(@errorlevels, $line);
	    close PARALLELEXEFH2;
	}
	$rethash{errorlevels} = \@errorlevels;

	my @outs;
	for(my $i = 0; $i < @commands; $i++) {
	    open PARALLELEXEFH2, "<$stdoutfile[$i]" or last SETRETVAL;
	    my @outputstrings = <PARALLELEXEFH2>;
	    push(@outs, \@outputstrings);
	    close PARALLELEXEFH2;
	}
	$rethash{stdouts} = \@outs;

	$rethash{error} = 0;
    }
    unlink $zshfile;
    for(my $i = 0; $i < @commands; $i++) {
	    unlink $stdoutfile[$i];
	    unlink $statusfile[$i];
    }
    return %rethash;
}

1;
