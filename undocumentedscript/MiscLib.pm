=pod

=head1 NAME

MiscLib.pm - Miscellaneous library

=head1 SYNOPSIS

	use MiscLib;

	my $a = 3;
	assert($a == 3, "a must be 3"); # croak unless $a == 3

	my $b = abs(-4); # 4

	my $maxvalue = max(1, 2, 3, 4, 5); # 5

	my $minvalue = min(1, 2, 3, 4, 5); # 1

	my @uniqueval = unique(4, 4, 1, 5, 3, 2, 5, 3, 2, 3);
	print join(',', sort(@uniqueval)),"\n"; # 1, 2, 3, 4, 5

	my @arr = ("ab", "foo", "bar", "baz");
	print join(',', reverse(@arr)),"\n"; # baz,bar,foo,ab

	my $blastpath = getExecutablePath("blastall");
	die "'blastall' is not found" unless(defined $blastpath);
	system "$blastpath -d nr -i input.fa -p blastn";

    if(isWindows()) {
    	print "This script is executed on Windows\n";
    } else {
    	print "This script is executed on non-Windows(maybe on UNIX system)\n";
    }

    my @dirs = printAndExecute("ls -al");   # array context. croak on error
    my $path = printAndExecute("which ls"); # scalar context. croak on error

    print tolower("AbCd1"), "\n"; # prints "abcd1\n";
    print toupper("AbCd1"), "\n"; # prints "ABCD1\n";
    print "Yes" if(!islower("a2s")); # prints "Yes";
    print "Yes" if(islower("as")); # prints "Yes";
    print "Yes" if(isupper("YZ")); # prints "Yes";
    print "Yes" if(!isupper("Yz")); # prints "Yes";

    my $cpath = get_canonical_path_name("./abc//defgag////./efasdf/f/");
    print "Yes" if($cpath eq './abc/defgag/efasdf/f'); # prints "Yes";
    my $cpath2 = get_canonical_path_name("/usr/local/bin/../lib/libz.so");
    print "Yes" if($cpath2 eq '/usr/local/lib/libz.so'); # prints "Yes";

	my $rpath = get_trunk_file_name("/usr/local/bin/hogehoge");
	print "Yes" if($rpath eq 'hogehoge'); # prints "Yes";

	my $dir  = get_dir_name("/usr/local/bin/hogehoge");
	print "Yes" if($dir eq '/usr/local/bin'); # prints "Yes";

	my $wpath = get_windows_path("c:/usr/local/bin/a.exe");
	print "Yes" if($dir eq "c:\\use\\local\\bin\\a.exe"); # prints "Yes";

   	my $upath = get_unix_path("c:\\use\\local\\bin\\a.exe");
	print "Yes" if($dir eq "c:/usr/local/bin/a.exe"); # prints "Yes";

	my $abspath = get_absolute_path("hogehoge\abc.exe");
	print "Yes" if($dir eq "c:\\currentdirectory\\hogehoge\\abc.exe"); # On Windows

   	my $abspath = get_absolute_path("bin/abc.pl");
	print "Yes" if($dir eq "/currentdirectory/bin/abc.pl"); # On UNIX

=head1 DESCRIPTION

B<MiscLib.pm> have various useful functions. Please see the synopsis to
understand the usage.

=cut

use Carp;
use Cwd;

sub assert($$)
{
    my $value = shift;
    my $message = shift;
    if($value == 0) {
    	unless($message eq '') {
        	croak "Assertion failed '$message'";
        } else {
        	croak "Assertion failed";
        }
    }
}

sub abs($)
{
    my $val = shift;
    return $val if($val >= 0);
    return -$val;
}

sub max(@)
{
	my $result = shift;
	for(@_) {
    	$result = $_ if($_ > $result);
	}
	return $result;
}

sub min(@)
{
	my $result = shift;
	for(@_) {
		$result = $_ if($_ < $result);
	}
	return $result;
}

sub unique(@)
{
	my %count;
	return grep(!$count{$_}++, @_);
}

sub reverse(@)
{
	my @retval;
	my $len = @_;
	my $lenminusone = $len - 1;
	for(my $i = 0; $i < $len; $i++) {
		$retval[$lenminusone - $i] = $_[$i];
	}
	return @retval;
}

sub ensureExecutable($)
{
	my $exefilename = shift;
	my $path = `which $exefilename`;
	chomp $path;
	chop $path if($path =~ /\r$/);
	return $path if($^O =~ /win/i);
	return undef unless(-x $path);
	return $path;
}

sub getExecutablePath($)
{
	my $v = shift;
	return ensureExecutable($v);
}

sub isWindows()
{
	return 1 if($ENV{'OS'} =~ /Win/);
	return 0;
}

sub printAndExecute($)
{
	my $commandLine = shift;
	print STDERR "% $commandLine\n";
	if(wantarray) {
		my @result = `$commandLine`;
		croak "Error occured while executing the command line above.\nError code($?)\nError message:$!" if($?);
		return @result;
	} else {
		my $result = `$commandLine`;
		croak "Error occured while executing the command line above.\nError code($?)\nError message:$!" if($?);
		return $result;
	}
}

sub islower($)
{
	my $char = shift;
	return $char =~ /^[a-z]+$/;
}

sub isupper($)
{
	my $char = shift;
	return $char =~ /^[A-Z]+$/;
}

sub tolower($)
{
	my $str = shift;
	$str =~ tr/A-Z/a-z/;
	return $str;
}

sub toupper($)
{
	my $str = shift;
	$str =~ tr/a-z/A-Z/;
	return $str;
}

sub get_canonical_path_name($)
{
	my $pathName = shift;
	$pathName .= '/';
	$pathName = './' . $pathName unless($pathName =~ m|^/|);
	$pathName =~ s|\\|/|g;
	$pathName =~ s|//|/|g;
	$pathName =~ s|/\./|/|g;
	$pathName =~ s!/([^/\.]+|\.[^/\.]+|\.\.[^/\.]+)/\.\./!/!g;
	$pathName =~ s|/$||;
	$pathName =~ s|^\./||;
	return $pathName;
}

sub get_trunk_file_name($)
{
	my $pathName = shift;
	$pathName =~ s|\\|/|g;
	$pathName =~ s|^.*/||;
	return $pathName;
}

sub get_dir_name($)
{
	my $pathName = shift;
	$pathName =~ s|\\|/|g;
	if($pathName =~ m|/|) {
		$pathName =~ s|/[^/]+$||;
		return $pathName;
	} else {
		return './';
	}
}

sub get_windows_path($)
{
	my $pathName = shift;
	$pathName =~ s|/|\\|g;
	return $pathName;
}

sub get_unix_path($)
{
	my $pathName = shift;
	$pathName =~ s|\\|/|g;
	return $pathName;
}

sub get_absolute_path($)
{
	my $pathName = shift;
	if(isWindows()){
		return $pathName if($pathName =~ m|^[A-Z]\:[\\/]|);
		my $cwd = getcwd();
		if($cwd =~ /[\\\/]$/) {
			return $cwd . $pathName;
		} else {
			return $cwd . "\\" . $pathName;
		}
	} else {
		return $pathName if($pathName =~ m|^/|);
		my $cwd = getcwd();
		if($cwd =~ /[\\\/]$/) {
			return $cwd . $pathName;
		} else {
			return $cwd . "/" . $pathName;
		}
	}
}

1;

