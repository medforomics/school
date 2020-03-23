#!/usr/bin/perl -w
#merge_genefusions.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'rnafile|f=s','delly|d=s','svaba|s=s','tid|t=s',
			  'pid|p=s','datadir|r=s','pindel|i=s');

open OM, "<$opt{datadir}/known_genefusions.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $known{$line} = 1;
}
close OM;

open OM, "<$opt{datadir}/panelgenes.txt" or die $!;
while (my $line = <OM>) {
    chomp($line);
    $keep{$line} = 1;
}

my %dnagf;

my @files = ('delly','svaba','pindel');

foreach my $method (@files) {
  next unless (-e $opt{$method});
  open ALG, "<$opt{$method}" or warn $!;
  my $head1 = <ALG>;
  chomp($head1);
  my ($H1,$H2,$H3,$H4,$H5,$H6,$H7,$H8,$H9,@sids) = split(/\t/,$head1);
  my %done;
  while (my $line = <ALG>) {
    chomp($line);
    my ($chr1,$p1,$chr2,$p2,$effect,$gene,$gtype,$filter,$format,@gts) = split(/\t/,$line);
    if ($filter =~ m/LOWMAPQ|LowQual/i) {
      $filter = 'FailedQC;'.$filter;
    }
    my @deschead = split(/:/,$format);
    unless ($gene =~ m/\w+/) {
	$gene = $gtype;
    }
    next unless $gene;
    next if ($gene =~ m/LINC00486/);
    $gene =~ s/&/--/g;
    my (@genes) = split(/--/,$gene);
    my $lgene = shift @genes;
    my $rgene = pop @genes;
    $chr2 =~ tr/CHR/chr/;
    unless ($chr2 =~ m/^chr/) {
	$chr2 =~ m/(chr\w+):(\d+)/;
	$chr2=$1;
	$p2 = $2;
    }
    $lbkpnt = join(":",$chr1,$p1);
    $rbkpnt = join(":",$chr2,$p2);
    unless ($rgene) {
	$rgene="intergenic";
    }
    next if ($done{$lbkpnt}{$rbkpnt});
    $done{$lbkpnt}{$rbkpnt} = 1;
    $done{$rbkpnt}{$lbkpnt} = 1;
    if ($chr1 eq $chr2) {
      $chrtype = 'INTRACHROMOSOMAL';
      $chrdistance = abs($p2 - $p1);
    }else {
      $chrtype = 'INTERCHROMOSOMAL';
      $chrdistance = join("--",$chr1,$chr2);
    }
    my %hash;
    next unless ($keep{$lgene} || $keep{$rgene} || $gtype =~ m/IG/);
    foreach my $i (0..$#sids) {
      $sids[$i] = (split(/\./,$sids[$i]))[0];
      my @gtinfo = split(/:/,$gts[$i]);
      my %gtdata;
      foreach my $i (0..$#deschead) {
	$gtdata{$deschead[$i]} = $gtinfo[$i];
      }
      if ($sids[$i] eq $opt{tid}) {
	next unless $gtdata{AO} =~ m/\d+/;
	$hash{tumor} = {FusionName=>$gene,LeftGene=>$lgene,
			LeftBreakpoint=>$lbkpnt,RightGene=>$rgene,
			RightBreakpoint=>$rbkpnt,DNAReads=>$gtdata{AO},
			Filter=>$filter,ChrType=>$chrtype,
			ChrDistance=>$chrdistance};
      }else {
	$hash{normal} = $gtdata{AO};
      }
    }
    next unless $hash{tumor};
    next if ($dnagf{$lgene}{$rgene} && $hash{tumor}{DNAReads} ne '.' && $dnagf{$lgene}{$rgene}{DNAReads} > $hash{tumor}{DNAReads});
    $dnagf{$lgene}{$rgene} = $hash{tumor};
    my %filter = map {$_=>1} split(";",$hash{tumor}{Filter});
    $ftype = 'Somatic';
    if ($hash{normal} > 10) {
      $ftype = 'Germline;NormalDNAReads:'.$hash{normal};
    }
    my @filter = sort {$a cmp $b} keys %filter;
    $dnagf{$lgene}{$rgene}{Filter} = join(";",@filter,$ftype);
  }
  close ALG;
}

open OUT, ">$opt{pid}.translocations.answer.txt" or die $!;
my @outcol= ("FusionName","LeftGene","LeftBreakpoint",
	     "LeftGeneExons","LeftStrand","RightGene",
	     "RightBreakpoint","RightGeneExons","RightStrand",
	     "RNAReads","DNAReads","FusionType","Annot",
	     'Filter','ChrType','ChrDistance');
print OUT join("\t",@outcol),"\n";

my %reported;
if ($opt{rnafile} && -e $opt{rnafile}) {
  open RNA, "<$opt{rnafile}" or die $!;
  my $header = <RNA>;
  chomp($header);
  my @colnames = split(/\t/,$header);
  while (my $line = <RNA>) {
    chomp($line);
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#colnames) {
	$hash{$colnames[$i]} = $row[$i];
    }
    if ($dnagf{$hash{LeftGene}}{$hash{RightGene}}) { 
	$hash{DNAReads} = $dnagf{$hash{LeftGene}}{$hash{RightGene}}{DNAReads};
	$reported{$hash{LeftGene}}{$hash{RightGene}} = 1;
	$reported{$hash{RightGene}}{$hash{LeftGene}} = 1;
    } elsif ($dnagf{$hash{RightGene}}{$hash{LeftGene}}) {
	$hash{DNAReads} = $dnagf{$hash{RightGene}}{$hash{LeftGene}}{DNAReads};
	$reported{$hash{LeftGene}}{$hash{RightGene}} = 1;
	$reported{$hash{RightGene}}{$hash{LeftGene}} = 1;
    }
    my @line;
    foreach (@colnames) {
	push @line, $hash{$_};
    }
    print OUT  join("\t",@line),"\n";
  }
}
foreach my $lgene (sort {$a cmp $b} keys %dnagf) {
    foreach my $rgene (sort {$a cmp $b} keys %{$dnagf{$lgene}}) {
	my $fname = $dnagf{$lgene}{$rgene}{FusionName};
	next if ($reported{$lgene}{$rgene});
	next unless ($known{$fname} ||  $dnagf{$lgene}{$rgene}{DNAReads} > 20);
	my ($lchr,$lpos) = split(/:/,$dnagf{$lgene}{$rgene}{LeftBreakpoint});
	my ($rchr,$rpos) = split(/:/,$dnagf{$lgene}{$rgene}{RightBreakpoint});
	next if ($lchr eq $rchr && abs($lpos-$rpos) < 9999);
	my @line;
	foreach (@outcol) {
	    $dnagf{$lgene}{$rgene}{$_} = '' unless $dnagf{$lgene}{$rgene}{$_};
	    push @line, $dnagf{$lgene}{$rgene}{$_};
	}
	print OUT join("\t",@line),"\n"
    }
}

close OUT;
