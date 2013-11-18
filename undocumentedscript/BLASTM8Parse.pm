# my @matches = parseBLASTM8ResultFile("blastm8.txt");
# for(@matches) {
#   my $query_id   = $_->{qid};
#   my $subject_id = $_->{sid};
#   ... so on.
# }

sub parseBLASTM8ResultFile($) {
    my $m8file = shift;
    open  BLASTM8, "< $m8file" or return undef;
    my @retval;
    while(<BLASTM8>) {
	my $br = parseBLASTM8ResultLine($_);
	push(@retval, $br);
    }
    close BLASTM8;
    return @retval;
} 

sub parseBLASTM8ResultLine($) {
    my $line = shift;
    my @args = split(/\t/);
    my %rval;
    $rval{'qid'}      = $args[ 0]; # Query ID
    chop $rval{'qid'} if($rval{'qid'} =~ /\|$/);
    $rval{'sid'}      = $args[ 1]; # Subject ID
    chop $rval{'sid'} if($rval{'sid'} =~ /\|$/);
    $rval{'identity'} = $args[ 2]; # Identity
    $rval{'alignlen'} = $args[ 3]; # Alignment length
    $rval{'mismatch'} = $args[ 4]; # Number of mismatch
    $rval{'gapopen'}  = $args[ 5]; # Number of gap opening
    $rval{'qstart'}   = $args[ 6]; # Query start position   (1-origin)
    $rval{'qend'}     = $args[ 7]; # Query end position     (1-origin)
    $rval{'sstart'}   = $args[ 8]; # Subject start position (1-origin)
    $rval{'send'}     = $args[ 9]; # Subject end position   (1-origin)
    $rval{'eval'}     = $args[10]; # e-value
    $rval{'bitscore'} = $args[11]; # bit-score
    return \%rval;
}

sub createBLASTM8Line($) {
    my $parsedobj = shift;
    my @line = ( $parsedobj->{'qid'},
		 $parsedobj->{'sid'},
		 $parsedobj->{'identity'},
		 $parsedobj->{'alignlen'},
		 $parsedobj->{'mismatch'},
		 $parsedobj->{'gapopen'},
		 $parsedobj->{'qstart'},
		 $parsedobj->{'qend'},
		 $parsedobj->{'sstart'},
		 $parsedobj->{'send'},
		 $parsedobj->{'eval'},
		 $parsedobj->{'bitscore'}
		);
    return join("\t", @line);
}

1;
