#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use XML::Parser;

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

sub isempty { my $s = shift; return (!defined $s || $s eq ''); }

if((isempty($search_pattern) && !$flag_index) || isempty($input_filename)) {
    print STDERR "usage: xmlex.pl [options] <pattern> <input XML>\n";
    print STDERR "See perldoc for details\n";
    exit 1;
}
$search_pattern = '' unless(defined $search_pattern);
$flag_xml = 1 if($flag_newline);

my $pattern_nodename;
my @pattern_attributeconstraint;
my $hasMatchedSomething = 0;
if($search_pattern ne '') {
    # <tagname>[@<constraint>[,<constraint>...]]
    # A second '@' is a syntax error, not part of the attribute name.
    if($search_pattern =~ /^([^@]+)(?:@([^@]*))?$/) {
	$pattern_nodename = $1;
	print STDERR "node name = '$pattern_nodename'\n" if($debug);
	my $attrconststr = $2;
	my @attrconstarr = defined $attrconststr ? split(/,/, $attrconststr) : ();
	for(@attrconstarr) {
	    # <attr>            attribute must be present
	    # <attr>=<value>    attribute must be exactly <value>
	    # <attr>~=<regexp>  attribute must match the regular expression
	    if(/^([^=~]+)(?:(~?=)(.*))?$/) {
		my $attrname = $1;
		my $operator = $2;
		my $value    = $3;
		print STDERR "attr '$attrname', op '", (defined $operator ? $operator : '(exists)'),
		             "', value '", (defined $value ? $value : ''), "'\n" if($debug);
		if(!defined $operator || !defined $value || $value eq '') {
		    # bare attribute name (or a trailing '='): existence test
		    push(@pattern_attributeconstraint, {key => $attrname, op => undef});
		} elsif($operator eq '=') {
		    push(@pattern_attributeconstraint, {key => $attrname, op => '=', val => $value});
		} else {
		    my $compiled = eval { qr/$value/ };
		    unless(defined $compiled) {
			print STDERR "Invalid regular expression '$value' in attribute constraint '$_':\n";
			print STDERR "  $@";
			exit 1;
		    }
		    push(@pattern_attributeconstraint, {key => $attrname, op => '~=', re => $compiled});
		}
	    } else {
		print STDERR "Invalid attribute constraint '$_'\n";
		exit 1;
	    }
	}
    } else {
	print STDERR "Invalid pattern '$search_pattern'\n";
	print STDERR "Expected '<tagname>' or '<tagname>\@<attr>[=<value>][,...]'\n";
	exit 1;
    }
}

my $index_filename = "${input_filename}.index";

# Index mode has never been implemented: the creation code used to print
# "Done." and then die, and the lookup code just printed a message and exited
# successfully, producing no output at all. Fail loudly when it is requested,
# and simply ignore a stale index file otherwise.
if($flag_index) {
    print STDERR "ERROR: index mode (-index) is not implemented.\n";
    print STDERR "       (the index file would be '$index_filename')\n";
    if($debug) {
	for(@param_attributes) {
	    print STDERR "Specified attr. : $_\n";
	}
    }
    exit 1;
}

unless(-r $input_filename) {
    print STDERR "ERROR: cannot open the input XML file '$input_filename'",
                 (-e $input_filename ? " (not readable)" : " (no such file)"), "\n";
    exit 1;
}

if(-e $index_filename) {
    print STDERR "NOTE: '$index_filename' exists, but index mode is not implemented.\n";
    print STDERR "      The index file is ignored and '$input_filename' is scanned normally.\n";
}

{
    my @start_tag_pile;
    my @start_tag_pile_tobeoutput;
    my @end_tag_pile;
    my @end_tag_pile_tobeoutput;
    my $isShowing = 0;
    my $hasShownHeader = 0;

    sub elemTextHandler {
      my $self = shift;
      my $text = shift;
      if($isShowing) {
	  $text =~ s/\s+//g if($flag_eliminatespace);
	  # XML::Parser hands us the *unescaped* text, so it has to be
	  # escaped again to stay well-formed XML.
	  print xmlEscapeText($text);
      }
    }
    sub elemEndHandler {
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
    sub elemStartHandler {
      my ($self, $name, %attrs) = @_;
      my $tagstr = $name;
      for my $k (sort keys %attrs) {
	  $tagstr .= ' ' . $k . '="' . xmlEscapeAttributeValue($attrs{$k}) . '"';
      }
      push(@start_tag_pile, $tagstr);
      push(@start_tag_pile_tobeoutput, 1);
      push(@end_tag_pile, $name);
      push(@end_tag_pile_tobeoutput, 0);
      
      if($isShowing) {
	  print "<$tagstr>";
	  $isShowing++;
      } else {
	  if(defined $pattern_nodename && $pattern_nodename eq $name) {
	      my $has_matched_the_criteria = 1;
	      for my $attrconst (@pattern_attributeconstraint) {
		  my $attrname = $attrconst->{key};
		  my $operator = $attrconst->{op};
		  my $actual   = $attrs{$attrname};
		  if(!defined $actual) {
		      $has_matched_the_criteria = 0; last;
		  }
		  if(!defined $operator) {
		      # existence test only; already satisfied
		  } elsif($operator eq '=') {
		      unless($actual eq $attrconst->{val}) {
			  $has_matched_the_criteria = 0; last;
		      }
		  } else {
		      unless($actual =~ $attrconst->{re}) {
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

# Escape a text node. XML::Parser unescapes the input, so '&' and '<' have to
# be put back before the text is written out again.
sub xmlEscapeText {
    my $s = shift;
    return '' unless(defined $s);
    $s =~ s/&/&amp;/g;
    $s =~ s/</&lt;/g;
    $s =~ s/>/&gt;/g;
    return $s;
}

# Same, for a value that goes inside a double quoted attribute.
sub xmlEscapeAttributeValue {
    my $s = shift;
    return '' unless(defined $s);
    $s =~ s/&/&amp;/g;
    $s =~ s/</&lt;/g;
    $s =~ s/>/&gt;/g;
    $s =~ s/"/&quot;/g;
    return $s;
}

sub extract_header_from_xml {
    my $filename = shift;
    open(my $fh, '<', $filename) or die "Cannot open '$filename': $!";
    my $restline = 10;
    while(<$fh>) {
	# print only the XML declaration itself: printing the whole line would
	# emit the root element too when the two share a line.
	if(/(<\?.*?\?>)/) {
	    print "$1\n";
	    last;
	}
	$restline--;
	last if($restline <= 0);
    }
    close $fh;
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

The pattern is C<< <tagname> >>, optionally followed by C<@> and a comma
separated list of attribute constraints:

=over 4

=item C<attr>

the attribute must be present, whatever its value

=item C<< attr=value >>

the attribute must be exactly C<value> (a plain string comparison, so
C<name=US> does B<not> match C<USA>)

=item C<< attr~=regexp >>

the attribute must match the Perl regular expression C<regexp>. The
expression is unanchored, so C<< name~=US >> matches C<USA> as well; use
C<< name~=^US$ >> for an exact match. An invalid expression is reported
instead of aborting with a regular expression compilation error.

=back

Only one C<@> is allowed in a pattern.

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
    xmlex.pl -xml -delspace country@name=Japan abc.xml
    <hello><world><country name="Japan"><prefecture name="Nagano"></prefecture></country></world></hello>

Exmaple2 (regular expression; matches both USA and Japan):
    xmlex.pl -delspace 'country@name~=[Aa]' abc.xml
    <country name="USA"><state name="Verginia"></state></country><country name="Japan"><prefecture name="Nagano"></prefecture></country>

=cut
