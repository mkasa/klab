#!/usr/bin/env perl

use strict;

use Getopt::Long;
use Pod::Usage;
use XML::Parser;
use File::stat;
use DBI;
# also DBD::SQLite is nececcary.

my $flag_man            = 0;
my $flag_help           = 0;
my $flag_eliminatespace = 0;
my $flag_xml            = 0;
my $flag_newline        = 0;
my $flag_index          = 0;
my $debug               = 0;
my @param_attributes    = ();

GetOptions( 'help|?'       => \$flag_help,
	    'man'          => \$flag_man,
	    'debug'        => \$debug,
	    'delspace'     => \$flag_eliminatespace,
	    'xml'          => \$flag_xml,
	    'newline'      => \$flag_newline,
	    'index'        => \$flag_index,
	    'attributes=s' => \@param_attributes,
	    ) or pod2usage(2);
pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man; 

my $search_pattern = shift;
my $input_filename = shift;

if(($search_pattern eq '' && !$flag_index) || $input_filename eq '') {
    print STDERR "usage: xmlex.pl [options] <pattern> <input XML>\n";
    print STDERR "See perldoc for details\n";
    exit 1;
}
$flag_xml = 1 if($flag_newline);

my $pattern_nodename;
my @pattern_attributeconstraint;
my $hasMatchedSomething = 0;
if($search_pattern ne '') {
    if($search_pattern =~ /^([^@]+)(@(.*))?$/) {
	$pattern_nodename = $1;
	print STDERR "node name = '$pattern_nodename'\n" if($debug);
	my $attrconststr = $3;
	my @attrconstarr = split(/,/, $attrconststr);
	for(@attrconstarr) {
	    if(/^([^=]+)(=(.*))?$/) {
		my $attrname = $1;
		my $value = $3;
		print STDERR "attr '$attrname', value '$value'\n" if($debug);
		push(@pattern_attributeconstraint, {key => $attrname, val => $value});
	    } else {
		print STDERR "Invalid attribute constraint '$_'\n";
		exit 1;
	    }
	}
    } else {
	print STDERR "Invalid pattern '$search_pattern'\n";
	exit 1;
    }
}

my $index_filename = "${input_filename}.index";

if($flag_index && !-e $index_filename) {
    print STDERR "Creating index file ...\n";
    print STDERR "Done.\n";
    if($debug) {
	for(@param_attributes) {
	    print "Specified attr. : $_\n";
	}
    }
    die "Not implemented yet";
    my $databasehandle = DBI->connect("dbi:SQLite:dbname=${index_filename}","","");
    my $result  = $databasehandle->do('create table xmlexidxattr (
                                               tagname       text,
                                               attributename text
                                       );');
    #	|| die $databasehandle->errstr; # error is OK. (for already existing database)
    $result = $databasehandle->do('PRAGMA default_synchronous = OFF;') || die $databasehandle->errstr;
    $result = $databasehandle->do('PRAGMA synchronous = OFF;')         || die $databasehandle->errstr;
    my $insert_sql_statement = "insert into xmlexidxattr values (?, ?);";
    my $insert_sql_handle    = $databasehandle->prepare($insert_sql_statement);
    for(@param_attributes) {
	if(/^([^@]+)(@(.*))?$/) {
	    my $nodename = $1;
	    print STDERR "node name = '$nodename'\n" if($debug);
	    my $attrconststr = $3;
	    my @attrconstarr = split(/,/, $attrconststr);
	    for(@attrconstarr) {
		if(/^([^=]+)(=(.*))?$/) {
		    my $attrname = $1;
		    my $value = $3;
		    print STDERR "attr '$attrname', value '$value'\n" if($debug);
		    push(@pattern_attributeconstraint, {key => $attrname, val => $value});
		} else {
		    print STDERR "Invalid attribute constraint '$_'\n";
		    exit 1;
		}
	    }
	} else {
	    print STDERR "Invalid pattern '$search_pattern'\n";
	}
    }
    # TODO: write next;
}
exit if($flag_index && $search_pattern eq '');

my $index_mode = 0;
if(-e $index_filename) {
    my $inputfilestat = stat($input_filename);
    my $indexfilestat = stat($index_filename);
    if($inputfilestat->mtime <= $indexfilestat->mtime) {
	$index_mode = 1;
    } else {
	print STDERR "Index file '$index_filename' is created at ", scalar localtime $indexfilestat->mtime, "\n";
	print STDERR "Input file '$input_filename' is created at ", scalar localtime $inputfilestat->mtime, "\n";
	print STDERR "Index file is older. Ignores index file for safety.\n";
    }
}

if($index_mode) {
    print STDERR "Index mode is not yet implemented\n";
} else {
    my @start_tag_pile;
    my @start_tag_pile_tobeoutput;
    my @end_tag_pile;
    my @end_tag_pile_tobeoutput;
    my $isShowing = 0;
    my $hasShownHeader = 0;

    sub elemTextHandler ($$) {
      my $self = shift;
      my $text = shift;
      if($isShowing) {
	  $text =~ s/\s+//g if($flag_eliminatespace);
	  print $text;
      }
    }
    sub elemEndHandler ($$) {
      my $self = shift;
      my $name = shift;
      pop(@start_tag_pile_tobeoutput);
      pop(@start_tag_pile);
      pop(@end_tag_pile);
      my $tooutput_close = pop(@end_tag_pile_tobeoutput);
      my $tooutput = $tooutput_close;
      if($isShowing) {
	  $tooutput = 1;
	  --$isShowing;
      }
      print "</$name>" if($tooutput);
      print "\n" if($tooutput_close && $flag_newline);
    }
    sub elemStartHandler ($$%) {
      my ($self, $name, %attrs) = @_;
      my $tagstr = $name;
      while(my ($k, $v) = each %attrs) {
	  $tagstr .= " $k=\"$v\"";
      }
      push(@start_tag_pile, $tagstr);
      push(@start_tag_pile_tobeoutput, 1);
      push(@end_tag_pile, $name);
      push(@end_tag_pile_tobeoutput, 0);
      
      if($isShowing) {
	  print "<$tagstr>";
	  $isShowing++;
      } else {
	  if($pattern_nodename eq $name) {
	      my $has_matched_the_criteria = 1;
	      for my $attrconst (@pattern_attributeconstraint) {
		  my $name = $attrconst->{key};
		  my $val  = $attrconst->{val};
		  if($val eq '') {
		      unless(defined $attrs{$name}) {
			  $has_matched_the_criteria = 0; last;
		      }
		  } else {
		      if($attrs{$name} !~ $val) {
			  $has_matched_the_criteria = 0; last;
		      }
		  }
	      }
	      if($has_matched_the_criteria) {
		  if($flag_xml) {
		      unless($hasShownHeader) {
			  extract_header_from_xml($input_filename);
			  $hasShownHeader = 1;
		      }
		      for(my $i = 0; $i < @start_tag_pile_tobeoutput; $i++) {
			  if($start_tag_pile_tobeoutput[$i]) {
			      my $tstr = $start_tag_pile[$i];
			      print "<$tstr>";
			      print "\n" if($flag_newline && $i + 1 != @start_tag_pile_tobeoutput);
			      $start_tag_pile_tobeoutput[$i] = 0;
			  }
		      }
		      for(@end_tag_pile_tobeoutput) {
			  $_ = 1;
		      }
		  } else {
		      print "<$tagstr>";
		  }
		  $isShowing = 1;
		  $hasMatchedSomething = 1;
	      }
	  }
      }
    }

    my $parser = new XML::Parser(ErrorContext=>3);
    $parser->setHandlers(Start=>\&elemStartHandler,
                         Char=>\&elemTextHandler,
                         End=>\&elemEndHandler);
    $parser->parsefile($input_filename);
}
print "\n" if($hasMatchedSomething);

sub extract_header_from_xml($) {
    my $filename = shift;
    open EXTRACTHEADER, "< $filename" or die "Cannot open '$filename'";
    my $restline = 10;
    while(<EXTRACTHEADER>) {
	print if(/<\?(.*)\?>/);
	$restline--;
	last if($restline <= 0);
    }
    close EXTRACTHEADER;
}

=pod

=head1 NAME

xmlex.pl - XML extractor

=head1 SYNOPSIS

xmlex.pl [options] <pattern> <input XML file>

Options:
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-delspace>

Remove all spaces from XML text.

=item B<-xml>

Outputs header and all parent tags.

=item B<-newline>

Outputs parent tags with newline.

=item B<-index>

Use index mode.

=item B<-attribute>

Use attributes specifed by B<-attribute> when creating index.
For exmaple, giving -attribute=country@name will index <country> tags with
attribute 'name'. You can give multiple attributes by separating them by ','.
For example, givine -attribute=country@name,place will index tags like
<country name="USA" place="America">. Attributes that are not specified by
the option are just not indexed.
If you want to index various tags, just give -attributes multiple times.
Giving '-attribute=country@place,name -attribute=state@name' will
index <country place="sth" name="sth"> and <state name="sth">.

=back

=head1 DESCRIPTION

B<xmlex.pl> will read search the given pattern in XML.

=head1 EXAMPLES

Example XML file (abc.xml)
    <hello>
        <world>
             <country name="USA">
                  <state name="Verginia" />
             </country>
             <country name="Japan">
                  <prefecture name="Nagano" />
             </country>
        </world>
    </hello>

Exmaple1:
    xmlex.pl country@name=3 abc.xml
    <hello><world><country name="Japan">
    <prefecture name="Nagano" />
    </country></world></hello>

=cut
