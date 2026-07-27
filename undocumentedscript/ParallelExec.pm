# ParallelExec.pm - run a list of shell command lines in parallel, either
# through Sun Grid Engine (tge_make) or locally through zsh.
#
#   use ParallelExec qw(parallelExecute);
#   my %ret = parallelExecute('cmd1', 'cmd2', 'cmd3');
#
# Nothing is exported by default; import exactly the subs you use, or call
# them fully qualified as ParallelExec::parallelExecute(...). The latter is
# what a caller that loads this module with 'require' (rather than 'use')
# has to do, because 'require' never runs import().
#
# NOTE: TGEW.pm (Sun Grid Engine wrapper) is loaded lazily, only when we are
#       about to decide whether the SGE backend can be used. It used to be a
#       hard 'use TGEW;', which made this module impossible to compile - and
#       therefore the SMP backend impossible to use - on any machine without
#       TGEW installed.

package ParallelExec;

use strict;
use warnings;

use File::Temp;
use File::Spec;
use FileHandle;
use Exporter 'import';

our @EXPORT_OK = qw(
    parallelExecute
    parallelExecute_SGE
    parallelExecute_SMP
    removeIOTags
);
our %EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

# Quote a string for /bin/sh.
sub _shell_quote($)
{
    my $s = shift;
    return "''" unless(defined $s && $s ne '');
    $s =~ s/'/'\\''/g;
    return "'" . $s . "'";
}

# Escape a shell command line so that it can be put in a makefile recipe.
# '$' is make's variable sigil, so it has to be doubled or make eats it.
sub _make_escape($)
{
    my $s = shift;
    return '' unless(defined $s);
    $s =~ s/\$/\$\$/g;
    return $s;
}

# Number of concurrent jobs to run in the SMP backend.
sub _cpu_count()
{
    if(defined $ENV{PARALLELEXEC_MAXJOBS} && $ENV{PARALLELEXEC_MAXJOBS} =~ /^\d+$/ && $ENV{PARALLELEXEC_MAXJOBS} > 0) {
        return $ENV{PARALLELEXEC_MAXJOBS} + 0;
    }
    my $n;
    if($^O =~ /^darwin/i || $^O =~ /bsd/i) {
        $n = `sysctl -n hw.ncpu 2>/dev/null`;
    } elsif(-r '/proc/cpuinfo') {
        my $count = 0;
        if(open(my $cfh, '<', '/proc/cpuinfo')) {
            while(<$cfh>) { $count++ if(/^processor\s*:/); }
            close $cfh;
        }
        $n = $count;
    } else {
        $n = `getconf _NPROCESSORS_ONLN 2>/dev/null`;
    }
    $n = '' unless(defined $n);
    chomp $n;
    return ($n =~ /^\d+$/ && $n > 0) ? $n + 0 : 1;
}

sub _is_sge_available()
{
    return 0 unless(eval { require TGEW; 1 });
    my $r = eval { TGEW::is_installed() };
    return $r ? 1 : 0;
}

sub parallelExecute(@) {
	if(_is_sge_available()) {
		return parallelExecute_SGE(@_);
	} else {
		return parallelExecute_SMP(@_);
	}
}

# parallel execution of given commands (for Sun Grid Engine environment)
# returns hash.
# $ret = parallelExecute(('command1', 'command2', 'command3'));
sub parallelExecute_SGE(@)
{
    my @commands = @_;
    my %rethash = ();
    $rethash{error} = 1;

    # All the temporary files live in a private (mode 0700) directory that is
    # removed when $tmpdirobj goes out of scope. The previous code passed a
    # RELATIVE template to tempfile(..., OPEN => 0), which neither created nor
    # locked the file and dropped temp*/temps* droppings into the caller's cwd.
    my $tmpdirobj = File::Temp->newdir();
    my $tmpdir    = "$tmpdirobj";

    my $makefile = File::Spec->catfile($tmpdir, 'parallelexec.mk');
    my @stdoutfile;
    my @stderrfile;
    my @statusfile;
    for(my $i = 0; $i < @commands; $i++) {
	    push(@stdoutfile, File::Spec->catfile($tmpdir, "out.$i"));
	    push(@stderrfile, File::Spec->catfile($tmpdir, "err.$i"));
	    push(@statusfile, File::Spec->catfile($tmpdir, "status.$i"));
    }

    {
	    my $fh = new FileHandle "> $makefile";
	    unless(defined $fh) {
		    print STDERR "Cannot create a temporary makefile '$makefile': $!\n";
		    return %rethash;
	    }
	    print $fh "all : ";
	    print $fh join(' ', @stdoutfile);
	    print $fh "\n\n";
	    for(my $i = 0; $i < @commands; $i++) {
	        # removeIOTags() must run BEFORE the '$' escaping, otherwise the
	        # $(I)/$(O) tags would turn into shell command substitutions.
	        my $cmd = _make_escape(removeIOTags($commands[$i]));
	        print $fh "$stdoutfile[$i] :\n";
	        # '-' so that make keeps going and every per-command exit status
	        # really gets recorded in its own status file.
	        print $fh "\t-$cmd > $stdoutfile[$i] 2> $stderrfile[$i]; echo \$\$? > $statusfile[$i]\n\n";
	    }
	    $fh->close();
    }

    my $tgepath = `which tge_make`;
    $tgepath = '' unless(defined $tgepath);
    chomp $tgepath;
    unless($tgepath ne '' && -x $tgepath) {
	    print STDERR "This function requires 'tge_make', which is not found on system.\n";
	    print STDERR "Please check if TGE is properly installed and on PATH.\n";
	    return %rethash;
    }
    my $commandline = _shell_quote($tgepath) . ' -f ' . _shell_quote($makefile) . ' -tgelock';
    $rethash{shellstdout} = `$commandline`;
    $rethash{shellerrorlevel} = $?;

    # print STDERR "ERRORLEVEL = $rethash{shellerrorlevel}\n";

    SETRETVAL: {
	    if($rethash{shellerrorlevel}) {
			print STDERR "ERROR CODE: $rethash{shellerrorlevel}\n";
			last;
		}
	    my @errorlevels;
	    for(my $i = 0; $i < @commands; $i++) {
	        push(@errorlevels, _read_status($statusfile[$i]));
	    }
	    $rethash{errorlevels} = \@errorlevels;

	    $rethash{stdouts} = _read_all_files(\@stdoutfile);
	    $rethash{stderrs} = _read_all_files(\@stderrfile);

	    $rethash{error} = 0;
        # print STDERR "NOERR\n";
    }
    return %rethash;
}

# Read a per-command exit status file. A killed command may leave the file
# empty; report that as a failure rather than as undef (which numeric
# comparisons would silently treat as success).
sub _read_status($)
{
    my $statusfile = shift;
    open(my $fh, '<', $statusfile) or return -1;
    my $line = <$fh>;
    close $fh;
    return -1 unless(defined $line);
    chomp $line;
    return -1 unless($line =~ /^-?\d+$/);
    return $line + 0;
}

sub _read_all_files($)
{
    my $filesref = shift;
    my @outs;
    for my $f (@$filesref) {
        if(open(my $fh, '<', $f)) {
            my @lines = <$fh>;
            close $fh;
            push(@outs, \@lines);
        } else {
            push(@outs, []);
        }
    }
    return \@outs;
}

sub removeIOTags($)
{
	my $cmdstr = shift;
	return '' unless(defined $cmdstr);
	$cmdstr =~ s/\$\([IO]\)//gi;
	return $cmdstr;
}

# parallel execution of given commands (for SMP/zsh environment)
# returns hash.
# $ret = parallelExecute_SMP(('command1', 'command2', 'command3'));
#
# At most _cpu_count() commands run at the same time; set the environment
# variable PARALLELEXEC_MAXJOBS to override.
sub parallelExecute_SMP(@)
{
    my @commands = @_;
    my %rethash = ();
    $rethash{error} = 1;

    my $zshpath = `which zsh`;
    $zshpath = '' unless(defined $zshpath);
    chomp $zshpath;
    unless($zshpath ne '' && -x $zshpath) {
	    print STDERR "This script requires 'zsh', which is not found on system.\n";
	    print STDERR "Please check if zsh is properly installed.\n";
	    print STDERR "zsh must be on PATH environmental variable.\n";
	    exit 1;
    }

    # See the comment in parallelExecute_SGE() about the temporary directory.
    my $tmpdirobj = File::Temp->newdir();
    my $tmpdir    = "$tmpdirobj";

    my $zshfile = File::Spec->catfile($tmpdir, 'parallelexec.zsh');
    my @stdoutfile;
    my @stderrfile;
    my @statusfile;
    for(my $i = 0; $i < @commands; $i++) {
	    push(@stdoutfile, File::Spec->catfile($tmpdir, "out.$i"));
	    push(@stderrfile, File::Spec->catfile($tmpdir, "err.$i"));
	    push(@statusfile, File::Spec->catfile($tmpdir, "status.$i"));
    }

    my $maxjobs = _cpu_count();
    {
        open(my $wfh, '>', $zshfile) or do {
            print STDERR "Cannot create a temporary script '$zshfile': $!\n";
            return %rethash;
        };
        print $wfh "#!$zshpath\n";
        for(my $i = 0; $i < @commands; $i++) {
	        my $cmd  = removeIOTags($commands[$i]);
	        my $sout = $stdoutfile[$i];
	        my $serr = $stderrfile[$i];
	        my $stat = $statusfile[$i];
	        print $wfh "($cmd; print \$? > $stat) > $sout 2> $serr &\n";
	        # Cap the number of simultaneous jobs; without this,
	        # parallelExecute(@thousand_commands) forked 1000 processes at once.
	        print $wfh "wait\n" if(($i + 1) % $maxjobs == 0 && $i + 1 < @commands);
        }
        print $wfh "wait\n";
        unless(close $wfh) {
            print STDERR "Cannot write a temporary script '$zshfile': $!\n";
            return %rethash;
        }
    }

    my $commandline = _shell_quote($zshpath) . ' ' . _shell_quote($zshfile);
    $rethash{shellstdout} = `$commandline`;
    $rethash{shellerrorlevel} = $?;
    SETRETVAL: {
	last if($rethash{shellerrorlevel});
	my @errorlevels;
	for(my $i = 0; $i < @commands; $i++) {
	    push(@errorlevels, _read_status($statusfile[$i]));
	}
	$rethash{errorlevels} = \@errorlevels;

	$rethash{stdouts} = _read_all_files(\@stdoutfile);
	$rethash{stderrs} = _read_all_files(\@stderrfile);

	$rethash{error} = 0;
    }
    return %rethash;
}

1;
