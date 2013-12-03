#!/usr/bin/env perl

use strict;

use PostScript::Simple;
use PostScript::Simple::EPS;
use Bio::Perl;
use Bio::Seq;
use Bio::Seq::SeqWithQuality;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

my $ratio = 0.00001;
my $xbias = 2.0;
my $ybias = 2.0;
my $dx = 0.0; # image is rotated, so dx means vertical difference
my $dy = 2.0; # dy means horizontal difference
my $max_number_of_matched_alignment = 2147483646;
my $gaugestep    =  50000; #  50kbp
my $biggaugestep = 100000; # 100kbp
my $scaffoldcutoffratio  = 0.01;
my $scaffoldcutofflength = 3000;
my $maxpileratio         = 0;
my $debug = 0;
my $flag_doblast = 0;
my $flag_man     = 0;
my $flag_help    = 0;
my $flag_nolandscape = 0;
my $flag_useblat = 0;
my $numcpu = 1;
my $fastMap = 0;
my $tileSize =0;
my $fontname = "Times-Roman";


# matchthreshold
# blast matches whose length is less than $matchthreshold are discarded
my $matchthreshold = 70;

GetOptions( 'ratio=f' => \$ratio,
	    'xbias=f' => \$xbias,
	    'ybias=f' => \$ybias,
	    'dx=f' => \$dy,
	    'dy=f' => \$dx,
	    'gstep=i' => \$gaugestep,
	    'bgstep=i' => \$biggaugestep,
	    'matchthres=i' =>\$matchthreshold,
            'scaffoldcutoff=f' => \$scaffoldcutoffratio,
	    'leastscaffold=i'  => \$scaffoldcutofflength,
	    'maxpileratio=f'   => \$maxpileratio,
	    'maxlines=i' => \$max_number_of_matched_alignment,
	    'blast'   => \$flag_doblast,
	    'man'     => \$flag_man,
	    'help|?'  => \$flag_help,
	    'debug'   => \$debug,
	    'nolandscape' => \$flag_nolandscape,
	    'blat'    => \$flag_useblat,
	    'numcpu=i'  => \$numcpu,
	    'tileSize=i' => \$tileSize,
	    'fastMap' => \$fastMap,
	    'fontname=s' => \$fontname,
	    ) or pod2usage(2);

my $genome            = shift;
my $assembly          = shift;
my $blastreport       = shift;
my $outputepsprefix   = shift;

$flag_help = 1 if($genome eq '' or $assembly eq '' or $blastreport eq '' or $outputepsprefix eq '');

pod2usage(1) if $flag_help;
pod2usage(-verbose => 2) if $flag_man; 

unless(-e "$blastreport"){
    print STDERR "Can not found BLAST report file '$blastreport' \n";
    if($flag_doblast) {
	print STDERR "'$blastreport' is not found, but -blast option is specified\n";
	print STDERR "Will try to create MegaBLAST report\n\n";
	unless(-e "${genome}.nhr" && -e "${genome}.nsq" && -e "${genome}.nin"){
	    print STDERR "Creating Database\n";
	    my $createdbcmd = "makeblastdb -dbtype nucl -in $genome -out $genome";
	    print STDERR "% $createdbcmd\n";
	    `$createdbcmd`;
	    if($? > 0) {
		print STDERR "Error occured while formatdb. Abort.\n";
		exit 1;
	    }
	} else {
	    print STDERR "Already createdbed\n";
	}
	print STDERR "Mega BLASTing...\n";
	my $blastcmdline = "blastn -db $genome -query $assembly -evalue 1e-4 -perc_identity 93 -outfmt 7 -out $blastreport -num_threads 16";
	print STDERR "% $blastcmdline\n";
	`$blastcmdline`;
	if($? > 0) {
	    print STDERR "Error occured while BLASTing. Abort.\n";
	    exit 1;
	}
	print STDERR "Done\n";
    } elsif($flag_useblat) {
	print STDERR "'$blastreport' is not found, but -blat option is specified\n";
	print STDERR "Will try to create MegaBLAST report by BLAT\n\n";
	my $cmdline;
	if($numcpu == 1) {
	    $cmdline = "blat -out=blast9 -extendThroughN";
	} else {
	    $cmdline = "parablat.pl -numcpu=$numcpu -- -out=blast9 -extendThroughN";
	}
	if($fastMap) {
	    $cmdline .= " -fastMap";
	}
	if($tileSize > 0) {
	    $cmdline .= " -tileSize=$tileSize";
	}
	$cmdline .= " $genome $assembly $blastreport";
	print STDERR "Execute : $cmdline\n" if($debug);
	if(system "$cmdline") {
	    die "Could not execute blat($cmdline)\n";
	}
    } else {
	exit 2;
    }
} else {
    if($flag_doblast) {
	print STDERR "You specified -blast option, but there exists '$blastreport' already\n";
	print STDERR "and it seems there is no need to run BLAST. Will use the existing file.\n";
    } elsif($flag_useblat) {
	print STDERR "You specified -blat option, but there exists '$blastreport' already\n";
	print STDERR "and it seems there is no need to run BLAT. Will use the existing file.\n";
    }
}

print "Genome FASTA   = $genome\n";
print "Assembly FASTA = $assembly\n\n";

my $genomesize = 0;
my %chromosome_length;
{
    my $seqin  = Bio::SeqIO->new('-file'   => "<$genome",
				 '-format' => 'fasta'     );
    while(my $seq = $seqin->next_seq) {
	my $seqname = $seq->display_id;
	print "Loaded '$seqname'\n" if ($debug);
	my $nseq = $seq->seq;
	$chromosome_length{$seqname} = length $nseq;
	print "leght = $chromosome_length{$seqname}\n" if($debug);
    }
}

my %scaffold_length;
{
    my $seqin  = Bio::SeqIO->new('-file'   => "<$assembly",
				 '-format' => 'fasta'     );
    while(my $seq = $seqin->next_seq) {
	my $seqname = $seq->display_id;
	my $nseq = $seq->seq;
	$scaffold_length{$seqname} = length $nseq;
    }
}

# for use in figure mode 2
my $flinestep = 10000;

# color table
#          0    1    2    3    4    5    6    7    8    9   10   11   12    13    14   ext,  oth
my @r = (  0,   0, 255, 255,   0, 255, 128, 128,   0, 128, 255, 128,   0,  255,  128,  255,    0);
my @g = (  0, 255, 255,   0, 255, 128, 128,   0, 128, 255, 128,   0, 255,  128,  255,    0,    0);
my @b = (255,   0,   0, 255, 255, 128,   0, 128, 255, 128,   0, 255, 128,  255,    0,    0,    0);
my $maxnormalcolor = 14;
my $extcolor = 15;
my $othcolor = 16;

# map<string, vector<BlastMatch* >* > matches;
my %matches;
# my @blastmatches = @$matches{'chromosome1'}; # vector<BlastMatch*>
# my $blastmatch   = $blastmatches[0];         # BlastMatch*
# my $sid          = $blastmatch->{sid}        # subject ID

{
        print STDERR "Parsing BLAST report...\n";
	open FH, "$blastreport" or die "Cannot open megaBLAST report '$blastreport'";
	my $num = 0;
	my $pushednum = 0;
	my @array;
	while(<FH>){
	    chomp;
	    next if(/^#/ or /^$/);
	    my %br = parseblastresult($_);
	    my $subjectID = $br{sid};
	    if ($br{'alignlen'} >= $matchthreshold && $scaffold_length{$br{qid}} >= $scaffoldcutofflength) {
		push(@array , { %br } );
		# push(@${$matches{$subjectID}}, { %br } );
	    }
	    $num++;
	}

	@array = sort { $b->{alignlen} <=> $a->{alignlen} } @array;
	my $i;
	for($i = 0 ; $i < @array ; $i++) {
	    my $subjectID = $array[$i]{sid};
	    my %br = %{$array[$i]};
	    print "$subjectID -- $br{alignlen} : $array[$i]{alignlen}\n" if($debug);

	    if(!(exists $matches{$subjectID}) || @${$matches{$subjectID}} < $max_number_of_matched_alignment) {
	        push(@${$matches{$subjectID}} , {%br} );
	        $pushednum++;
	    }
	#last if($pushednum >= $max_number_of_matched_alignment);
	}
	
	print STDERR "$num hits parsed\n";
	print STDERR "$pushednum hits survived after filtering out <${matchthreshold}bp matches / <$scaffoldcutofflength scaffolds \n";
	close FH;
}
{
        while(my ($key, $value) = each %matches){
	    my $matcharray = $$value;
	    my $len = scalar (@$matcharray);
	    print STDERR "$key has $len hits\n";

	}
}

# Draw EPS
{
        print STDERR "Drawing eps ...\n";
        # for(map::iter i = matches.begin(); i != matches.end(); i++){
        while(my ($key, $value) = each %matches){
	    # print "  Output $key\n";
	    # value : *vector<BlastMatch>
	    my $matcharray = $$value;

            # map<string/*scaffold name*/, int/*x-axis offset*/> xoffset;
            my %xoffset;
	    my $maxXOffset = 0;
            {
		my %matchedlengths;
		# map<string/*scaffold name*/, int/*sum of matched lengths*/> matchedlengths
		# print "Hoge\n";
		while(my ($key,$value) = each %scaffold_length) {
		    $xoffset{$key} = -1;
		}
	        foreach( @$matcharray ){
		    $matchedlengths{$_->{qid}} += $_->{alignlen};
		}
		my @matchedlengtharray;
		while(my ($k, $v) = each %matchedlengths) {
		    push(@matchedlengtharray,  $k);
		}
		@matchedlengtharray = sort { $matchedlengths{$b} <=> $matchedlengths{$a} } @matchedlengtharray;
		my $maxPileSize = $chromosome_length{$key} * $maxpileratio;
		print STDERR "Pile scaffolds up to $maxPileSize bp at most\n";
		for(@matchedlengtharray) {
		    if($matchedlengths{$_} < $scaffoldcutoffratio * $scaffold_length{$_} || ($maxPileSize != 0 && $maxPileSize < $maxXOffset) ) {
			$xoffset{$_} = -1;
		    } else {
			$xoffset{$_} = $maxXOffset;
			$maxXOffset += $scaffold_length{$_};
		    }
		}
	    }
	    
	    # pre-calculate maximum y-size
	    my $maxy = $chromosome_length{$key} * $ratio + $ybias;
	    print "maxy = $maxy, key = '$key'\n" if($debug);

	    # pre-calculate maximum x-size
 	    # Same routine as Drawing matched lines
 	    my $maxx = $xbias + $maxXOffset * $ratio;

	    print "paper size = $maxx+$xbias x $maxy+$ybias (XxY)\n" if ($debug);

	    my $xsize = $maxx + $xbias*2;
	    my $ysize = $maxy + $ybias*2;
	    my $number = 1;
	    my $p;
	    unless($flag_nolandscape) {
	        $p = new PostScript::Simple(colour => 1,
					    landscape => 1,
					    units => "cm",
					    xsize => $xsize,
					    ysize => $ysize,
					   );
	    } else {
	        $p = new PostScript::Simple(colour => 1,
					    units => "cm",
					    xsize => $xsize,
					    ysize => $ysize,
					   );
	    }

	    print STDERR "(xsize,ysize) = ($xsize,$ysize)\n" if ($debug);
	    
	    # fill rectangle for empty area(n-sequence(s) of Chromosome)
	    {
		# read the Chromosome data from file
		print STDERR "Filling N-sequence of Chromosome\n" if($debug);
		my $IN  = Bio::SeqIO->new('-file'   => "<$genome",
					  '-format' => 'fasta'     );
		while(my $seq = $IN->next_seq) {
		    my $seqname = $seq->display_id;
		    next if($seqname ne $key);
		    my $sequence = $seq->seq;
		    my $total_length=0;
		    $p->setcolour(128,128,128);
		    while($sequence =~ m/([nN]+)/){ # for each N-sequence
			my $position = index($sequence,$1);# find the starting index of N-sequence
			my $begin = $total_length + $position;

			my $stepLength = $position + length($1);
			$total_length += $stepLength;

			my $end = $total_length;
			my $line_center = $ybias + ($begin+ $end)*$ratio/2.0;
			$p->setlinewidth(($end-$begin)*$ratio);
			$p->line($xbias+$dx,$line_center+$dy,$maxx+$dx,$line_center+$dy);
                        # cut N-sequence and its previous sequence
			$sequence = substr($sequence,$stepLength);
			print STDERR "Chromosome N-sequence : $begin : $end\n" if ($debug);
		    }
		}
	    }

	    # fill rectangle for empty area(n-sequence(s) of Scaffold
	    {
		# read the Scaffold data from file
		my $IN  = Bio::SeqIO->new('-file'   => "<$assembly",
				 '-format' => 'fasta'     );
		while(my $seq = $IN->next_seq) {
		    my $seqname = $seq->display_id;
		    next if ($xoffset{$seqname} == -1);
		    my $sequence = $seq->seq;
		    my $total_length=0;
		    $p->setcolour(128,128,128);
		    while($sequence =~ m/([nN]+)/){ # for each N-sequence
			my $position = index($sequence,$1);
			my $begin = $total_length + $position + $xoffset{$seqname};

			my $stepLength = $position + length($1);
			$total_length += $stepLength;

			my $end = $total_length + $xoffset{$seqname};
			my $line_center = $xbias + ($begin + $end)*$ratio/2.0;
			$p->setlinewidth(($end-$begin)*$ratio);
			$p->line($line_center+$dx,$ybias+$dy,$line_center+$dx,$maxy+$dy);
			$sequence = substr($sequence,$stepLength);
			print STDERR "Scaffold N-sequence($seqname,) : $begin : $end\n" if ($debug);
		    }
		}
	    }
	    

	    # vertical line showing chromosome
	    $p->setcolour(0, 0, 0);
	    $p->setlinewidth(0.1);

	    $p->setfont("$fontname", 10);

	    # draw a line to separate each scaffold and its label(name of scaffold)
	    while(my ($scaffoldname, $xoffset) = each %xoffset) {
 		next if($xoffset == -1); # this scaffold is not matched
 		$p->line($xbias + $xoffset * $ratio +$dx, $ybias +$dy, $xbias + $xoffset * $ratio +$dx, $maxy +$dy);
 		$p->text( {rotate => 270, align => 'left'},
 			  $xbias * 1.3 + $xoffset * $ratio +$dx,
 			  $ybias * 1.3 +$dy,
 			  $scaffoldname);
 	    }

	    # print Chromosome axis and it's label
 	    $p->setfont("$fontname", 16);
 	    $p->line($xbias +$dx, $ybias +$dy, $xbias +$dx, $maxy +$dy);
 	    $p->text( {rotate => 270, align => 'centre'},
 		      $xbias * 0.3 +$dx,
 		      $chromosome_length{$key} * $ratio * 0.5 + $ybias * 0.5 +$dy,
 		      $key);

	    # Vertical gauge
	    # Draw the gauge of Chromosome axis
	    $p->setcolour(0, 0, 0);
	    for(my $base = 0; ; $base += $gaugestep){
		my $y = $maxy - $base * $ratio;
		last if($ybias > $y);
		if($base % $biggaugestep == 0){
		    $p->setlinewidth(0.08);
		    $p->line($xbias +$dx, $y +$dy, $xbias * 0.8 +$dx, $y +$dy);
		} else {
		    $p->setlinewidth(0.04);
		    $p->line($xbias +$dx, $y +$dy, $xbias * 0.85 +$dx, $y +$dy);
		}
	    }

 	    # Draw matched lines
	    my %flines;
	    my $fline     = 0;
	    foreach( @$matcharray ){
		# $_ = *BlastMatch
		my $matchlen = $_->{qend} - $_->{qstart};
		$matchlen = -$matchlen if($matchlen < 0);
		my $xoff = $xoffset{$_->{qid}};
		if($xoff != -1) {
		    my $sx = ($_->{qstart} + $xoff) * $ratio + $xbias;
		    my $sy =  $_->{sstart}          * $ratio + $ybias;
		    my $ex = ($_->{qend}   + $xoff) * $ratio + $xbias;
		    my $ey =  $_->{send}            * $ratio + $ybias;
		    my $scaffoldLineWidthThick = 8000;
		    my $scaffoldLineWidthThin = 400;
		    $p->setcolour($r[5], $g[5], $b[5]);
		    $p->setlinewidth($scaffoldLineWidthThick * $ratio);
		    $p->line($sx +$dx, $sy +$dy, $ex +$dx, $ey +$dy);
		    
		    $p->setcolour(128, 255, 128);
		    $p->setlinewidth($scaffoldLineWidthThin * $ratio);
		    $p->line($sx +$dx, $sy +$dy, $ex +$dx, $ey +$dy);
		}
	    }

	    # Draw Scaffolds line and its label(Scaffolds)
	    $p->setcolour(0, 0, 0);
	    $p->setlinewidth(0.1);
	    $p->line($xbias +$dx, $maxy +$dy, $maxx +$dx, $maxy +$dy);
	    $p->setfont("$fontname", 16);
	    $p->text( {align => 'centre'},
		      $maxx * 0.5 + $xbias * 0.5 +$dx,
		      $maxy + $ybias * 0.3 +$dy,
		      'Scaffolds');
	    
	    # Draw the gauges of Scaffolds axis
	    for(my $base = 0; ; $base += $gaugestep){
		my $x = $xbias + $base * $ratio;
		last if($maxx < $x);
		if($base % $biggaugestep == 0){
		    $p->setlinewidth(0.08);
		    $p->line($x +$dx, $maxy +$dy, $x +$dx, $maxy + $ybias * 0.20 +$dy);
		} else {
		    $p->setlinewidth(0.04);
		    $p->line($x +$dx, $maxy +$dy, $x +$dx, $maxy + $ybias * 0.15 +$dy);
		}
	    }

	    #draw a scale
	    $p->setcolour(0, 0, 0);
	    $p->setlinewidth(0.08);
	    my $biggaugey = $biggaugestep * $ratio;
	    $p->line($maxx + $xbias * 0.3 +$dx, $maxy * 0.2 - $biggaugey * 0.5 +$dy,
		     $maxx + $xbias * 0.6 +$dx, $maxy * 0.2 - $biggaugey * 0.5 +$dy);
	    $p->line($maxx + $xbias * 0.3 +$dx, $maxy * 0.2 + $biggaugey * 0.5 +$dy,
		     $maxx + $xbias * 0.6 +$dx, $maxy * 0.2 + $biggaugey * 0.5 +$dy);
	    $p->setlinewidth(0.04);
	    $p->line($maxx + $xbias * 0.45 +$dx, $maxy * 0.2 - $biggaugey * 0.5 +$dy,
		     $maxx + $xbias * 0.45 +$dx, $maxy * 0.2 + $biggaugey * 0.5 +$dy);
	    $p->setfont("$fontname", 14);
	    my $kbps = int($biggaugestep / 1000);
	    $p->text( {rotate => 270},
		      $maxx + $xbias * 0.3 +$dx,
		      $maxy * 0.2 - $biggaugey * 0.5 - $ybias * 0.2 +$dy,
		      "${kbps}kbp");
	    
 	    my $outputepsfilename = "$outputepsprefix.$key.eps";
 	    $outputepsfilename =~ s/\|/-/g;
 	    $p->output("$outputepsfilename");

	    print STDERR "Drawing has fished $outputepsfilename\n" if($debug);
	    
	    &FixEPSFile($outputepsfilename) if (!$flag_nolandscape);
        }
        # }
}

sub parseblastresult {
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
    $rval{'qend'}     = $args[ 7]; # Query end position     (1-origin)v
    $rval{'sstart'}   = $args[ 8]; # Subject start position (1-origin)
    $rval{'send'}     = $args[ 9]; # Subject end position   (1-origin)
    $rval{'eval'}     = $args[10]; # e-value
    $rval{'bitscore'} = $args[11]; # bit-score
    return %rval;
}

sub FixEPSFile{
    my $filename = shift;
    my $IN;
    my $OUT;
    my $buffer;
    my $xsize;
    my $ysize;
    
    print STDERR "Fixing EPS file $filename\n" if ($debug);
    open(IN, "<$filename");
    while(<IN>){
	$buffer = $buffer . $_;
	if($_ =~ m/^%%BoundingBox:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+).*/){
	    $xsize = $3;
	    $ysize = $4;
	    print STDERR "Filesize is detected.\n" if ($debug);
	    print STDERR "X-size is $xsize, Y-size is $ysize.\n" if ($debug);
	}
	elsif($_ =~ m/^\/landscape\s+.+/){
	    my $tmp = <IN>;
	    print STDERR "$tmp" if($debug);
	    $tmp =~ s/^\s+(\d+)\s+(\d+)\s+translate(.*)/ $xsize $2 translate$3/g;
	    print STDERR "$tmp" if($debug);
	    $buffer = $buffer . $tmp;
	}
    }
    close(IN);
    open(OUT,">$filename");
    print OUT $buffer;
    close(OUT);
}

__END__

=head1 NAME

blastreport2dotplot.pl - Creating dot plot from BLAST report file

=head1 SYNOPSIS

blastreport2dotplot.pl [options] E<lt>genome FASTAE<gt> E<lt>assembly FASTAE<gt> E<lt>BLAST reportE<gt> E<lt>output EPS prefixE<gt>

Options:
   -help            show brief help message
   -man             show full documentation

=head1 OPTIONS

=over 8

=item B<-ratio>

Set a scale of base pairs. -ratio=0.01 means 1pt corresponds to 100bp.

=item B<-matchthres>

Discards BLAST matches whose match ratio is less than specified value. 
B<-matchthres=200> will discards mathes less than 200bp. The default value is set to 300bp.

=item B<-scaffoldcutoff>

Discards matches poorly aligned to finished sequence.
A scaffold that aligns to finished sequence less than this threshold in ratio is discarded.
B<-scaffoldcutoff=0.3> means if only less than 30% of a scaffold aligns to given finished sequence, it will be discarded.

=item B<-blast>

This utility usually needs BLAST output. When this option is specified, the utility runs MegaBLAST and creates BLAST output automatically. You should run BLAST by yourself because this B<-blast> option does not support flexible parameter configuration of BLAST search.

=item B<-nolandscape>

Don't rotate the image.

=item B<-fontname>

Specify the font name in the EPS file.

=item B<-debug>

For debugging purposes only.
    
=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<blastreport2dotplot.pl> will read the given BLAST report file and creates dotplot between them.
It aims for comparing WGS assembly to its finished portion of genome.
You have to give finished sequence(s) for B<genome FASTA>, and WGS assembly for B<assembly FASTA>.
Dot plot pictures are output in eps format.

=cut
