#!/usr/bin/env perl

use strict;

use Getopt::Long;
use Pod::Usage;
use DBI;

my $flag_man      = 0;
my $flag_help     = 0;
my $flag_matepair = 1;
my $debug     = 0;

my $param_url          = 'http://default/';
my $param_species      = 'human';
my $param_revision     = '200401';
my $param_date         = 'default';
my $param_color        = '0,0,0';
my $param_trackname    = 'AlignedSequence';
my $param_trackcomment = '';
my $param_dburl        = 'http://default/query.cgi?name=';
my $param_leastlength  = 100; # bp
my $param_coverage     = 0.0;

GetOptions( 'help|?'     => \$flag_help,
	    'man'        => \$flag_man,
	    'debug'      => \$debug,
	    'url=s'      => \$param_url,
	    'species=s'  => \$param_species,
	    'revision=s' => \$param_revision,
	    'color=s'    => \$param_color,
	    'date=s'     => \$param_date,
	    'name=s'     => \$param_trackname,
	    'comment=s'  => \$param_trackcomment,
	    'dburl=s'    => \$param_dburl,
	    'leastlen=i' => \$param_leastlength,
	    'matepair'   => \$flag_matepair,
	    'leastcov=f' => \$param_coverage
	    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man; 

my $blastcsvfilename                = shift;
my $outputeasygenetrackfilename     = shift;
my $outputalignmentdatabasefilename = shift;

if($blastcsvfilename eq '' || $outputeasygenetrackfilename eq '' || $outputalignmentdatabasefilename eq '') {
    print STDERR "usage: blastcsv2easygenetrack.pl [options] <BLAST CSV> <output easy gene track> <alignment database for web>\n";
    print STDERR "   --help      show help\n";
    print STDERR "   --man       show manual\n";
    print STDERR "   --url=      set description URL\n";
    print STDERR "   --species=  set species. (ex: --species=medaka)\n";
    print STDERR "   --revision= set revision. (ex: --revision=200412)\n";
    print STDERR "   --color=    set background color. (ex: --color=255,255,255) \n";
    print STDERR "   --date=     set date. (ex: --date=20041108) \n";
    print STDERR "   --name=     set track name. (ex: --name=alignedSequence) \n";
    print STDERR "   --comment=  set track comment.\n";
    print STDERR "   --dburl=    set database URL of your track.\n";
    print STDERR "   --leastlen= set least length of alignments. (in bp) \n";
    exit 1;
}

if($param_date eq 'default') {
    $param_date = `date +%Y%m%d`;
    chomp $param_date;
}

if(-e $outputalignmentdatabasefilename) {
    print STDERR "Unlink existing $outputalignmentdatabasefilename (and its journal file)\n";
    unlink $outputalignmentdatabasefilename;
    unlink "${outputalignmentdatabasefilename}-journal";
}

print STDERR "Parameters\n";
print STDERR "  Least length of alignments : $param_leastlength bp\n";
print STDERR "  Description URL            : $param_url\n";
print STDERR "  Species                    : $param_species\n";
print STDERR "  Revision                   : $param_revision\n";
print STDERR "  Date                       : $param_date\n";
print STDERR "  Track Name                 : $param_trackname\n";
print STDERR "  Track Comment              : $param_trackcomment\n";
print STDERR "  Database URL               : $param_dburl\n";
print STDERR "\n";

my $databasehandle = DBI->connect("dbi:SQLite:dbname=$outputalignmentdatabasefilename","","");
{
    my $result  = $databasehandle->do('create table alignmentinfo (
                                       id integer primary key,
                                       targetname  varchar(64),
                                       queryname   varchar(64),
                                       targetstart integer,
                                       targetend   integer,
                                       querystart  integer,
                                       queryend    integer,
                                       targetalignstr text,
                                       queryalignstr  text
                                       );')	|| die $databasehandle->errstr;

    $result = $databasehandle->do('PRAGMA default_synchronous = OFF;') || die $databasehandle->errstr;
    $result = $databasehandle->do('PRAGMA synchronous = OFF;')         || die $databasehandle->errstr;
}

open OFH, "> $outputeasygenetrackfilename" or die "Cannot open '$outputeasygenetrackfilename'";

open FH, "< $blastcsvfilename" or die "Cannot open '$blastcsvfilename'";
print OFH "geneTrack name=$param_trackname comment=\"$param_trackcomment\" description_url=$param_url ";
print OFH "color=$param_color species=$param_species revision=$param_revision date=$param_date\n";

print STDERR "Processing ...\n";

my $insert_sql_statement = "insert into alignmentinfo values (?, ?, ?, ?, ?, ?, ?, ?, ?);";
my $insert_sql_handle    = $databasehandle->prepare($insert_sql_statement);

unless($flag_matepair) {
    print STDERR "Scaffold/Read mode \n";
    my $unique_id = 1;
    while(<FH>) {
	chomp;
	chop if(/\r$/);
	next if(/^#/);
	my ($qname, $qlen, $qstart, $qend, $sname, $slen, $sstart, $send, $mcount, $mlength, $qstring, $sstring)
	    = split(/,/);
	unless($qname eq '' || $sname eq '') {
	    if($mlength >= $param_leastlength && $mlength >= $param_coverage * $qlen) {
		my $disp_name = $qname;
		my $id = $unique_id;
		print OFH "gene target=$sname name=$disp_name range=$sstart,$send url=${param_dburl}?id=$id\n";
		my $result = $insert_sql_handle->execute($unique_id, $sname, $qname, $sstart, $send, $qstart, $qend, $sstring, $qstring);
		unless(defined $result) {
		    print STDERR "Warning : SQL error at ID $unique_id\n";
		}
		$unique_id++;
	    }
	}
    }
} else {
    print STDERR "Matepair mode \n";
    my %read2matcharray;
    while(<FH>) {
	chomp;
	chop if(/\r$/);
	next if(/^#/);
	my ($qname, $qlen, $qstart, $qend, $sname, $slen, $sstart, $send, $mcount, $mlength, $qstring, $sstring)
	    = split(/,/);
	unless($qname eq '' || $sname eq '') {
	    if($mlength >= $param_leastlength && $mlength >= $param_coverage * $qlen) {
		push(@{$read2matcharray{$qname}}, $_);
	    }
        }        
    }   
    # map<string(clonename), pair<forward matchstring, backward matchstring> >
    my %clonename2matchpair; 
    while(my ($readname, $matcharrayref) = each %read2matcharray) {
	my $bestmatchstring = undef;
	my $bestmatchmcount = 0;
	for(@{$matcharrayref}) {
	    my ($qname, $qlen, $qstart, $qend, $sname, $slen, $sstart, $send, $mcount, $mlength, $qstring, $sstring)
		= split(/,/);
	    if($bestmatchmcount < $mcount) {
		$bestmatchmcount = $mcount;
		$bestmatchstring = $_;
	    }
	}
	if(defined $bestmatchstring) {
	    if($readname =~ /^(.*)\.([bg])\d+/i) {
		my $clonename     = $1;
		my $directionchar = $2;
		if($directionchar =~ /^b$/i ) {
		    $clonename2matchpair{$clonename}->{forward}  = $bestmatchstring;
		} else {
		    $clonename2matchpair{$clonename}->{backward} = $bestmatchstring;
		}
	    }
	}
    }
    my $unique_id = 1;
    while(my ($clonename, $matchpairhashref) = each %clonename2matchpair) {
	my $is_f_available = defined $matchpairhashref->{forward};
	my $is_b_available = defined $matchpairhashref->{backward};
	my ($f_qname, $f_qlen, $f_qstart, $f_qend, $f_sname, $f_slen, $f_sstart, $f_send, $f_mcount, $f_mlength, $f_qstring, $f_sstring)
	    = split(/,/, $matchpairhashref->{forward});
	my ($b_qname, $b_qlen, $b_qstart, $b_qend, $b_sname, $b_slen, $b_sstart, $b_send, $b_mcount, $b_mlength, $b_qstring, $b_sstring)
	    = split(/,/, $matchpairhashref->{backward});
	my $is_f_plusstrand = ($f_sstart < $f_send) ^ ($f_qstart < $f_qend);
	my $is_b_plusstrand = ($b_sstart < $b_send) ^ ($b_qstart < $b_qend);

	my $are_bothends_available = $is_f_available && $is_b_available;
	my $are_bothends_on_the_same_scaffold = $f_sname eq $b_sname;
	my $are_bothends_facing_each_other    = ($is_f_plusstrand && !$is_b_plusstrand) || (!$is_f_plusstrand && $is_b_plusstrand);
	if($are_bothends_available && $are_bothends_on_the_same_scaffold && $are_bothends_facing_each_other) {
	    my $disp_name = $clonename;
	    my $id = $unique_id;
	    my $range_start = min($f_sstart, $f_send, $b_sstart, $b_send);
	    my $range_end   = max($f_sstart, $f_send, $b_sstart, $b_send);
	    print OFH "gene target=$f_sname name=$disp_name range=$range_start,$range_end exon=$f_sstart,$f_send exon=$b_sstart,$b_send\n";
	} else {
	    if($is_f_available) {
		my $strandchar = $is_f_plusstrand ? '+' : '-';
		print OFH "gene target=$f_sname name=$f_qname range=$f_sstart,$f_send strand=$strandchar\n";
	    }
	    if($is_b_available){
		my $strandchar = $is_b_plusstrand ? '+' : '-';
		print OFH "gene target=$b_sname name=$b_qname range=$b_sstart,$b_send strand=$strandchar\n";
	    }
	}
    }
}

print STDERR "Creating Index\n";
{
    my $result  = $databasehandle->do('create index alignmentinfo_idx on alignmentinfo (id);')
	|| die $databasehandle->errstr;
}
print STDERR "Done.\n";
$databasehandle->disconnect();

sub min()
{
    my $result = shift;
    while(my $v = shift) {
	$result = $v if($v < $result);
    }
    return $result;
}

sub max()
{
    my $result = shift;
    while(my $v = shift) {
	$result = $v if($v > $result);
    }
    return $result;
}


=pod

=head1 NAME

blastcsv2easygenetrack - Making 'easy gene track' from BLAST CSV

=head1 SYNOPSIS

blastcsv2easygenetrack.pl [options] <BLAST CSV> <output easy gene track> <alignment database for web>

Options:
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-url>

Specify URL that describes your track. Example is shown below.
	-url=http://www.somegene.com/query.cgi?gene=ABC1

=item B<-species>

Specify species of your track. Example is shown below.
	-species=human

=item B<-revision>

Specify revision of genome data of your track. Example is shown below.
	-revision=200412

=item B<-color>

Specify background color of your track. The color is specified in 'R,G,B' style. Example of 'background color black' is shown below.
	-color=0,0,0

=item B<-date>

Specify date of your track. Current date is used for default. Example is shown below.
	-date=20041105

=item B<-name>

Specify name of your track. Example is shown below.
	-name=AlignedSomeSeq

=item B<-comment>

Specify comment of your track.

=item B<-dburl>

Specify database URL of your track.  Example is shown below.
	-dburl=http://www.somegene.com/query.cgi?id=

=item B<-leastlen>

Specify least length of alignments. -leastlen=200 will filter out alignments of less than 200bp.
Default is 100.

=item B<-leastcov>

Specify least coverage of alignments. Coverages are calculated by dividing match length by query length.

=item B<-matepair>

Enables matepair mode. If you specify this option, forward and backward clones are paired to see
clone coverage easily.

=back

=head1 DESCRIPTION

B<blastcsv2easygenetrack.pl> will read the BLAST CSV input file and
output corresponding easy gene track.

=cut
