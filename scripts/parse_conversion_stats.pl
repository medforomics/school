#!/usr/bin/perl -w
#parse_conversion_stats.pl

use XML::LibXML;

my $dir = shift @ARGV;

my $xml = "$dir\/Stats/ConversionStats.xml";
$runname = (split(/\/+/,$xml))[-3];
my $parser = XML::LibXML->new();
my $doc = $parser->parse_file($xml);
my %info;
my %total;
foreach my $sample  ($doc->findnodes('Stats/Flowcell/Project/Sample')) {
  my $sampleName = $sample->findvalue('@name');
  next if ($sampleName eq 'all');
  $SampleMergeName = $sampleName;
  my @samplename = split(/_/,$sampleName);
  if ($samplename[-1] =~ m/^[A|B|C|D]$/ || $samplename[-1] =~ m/^[A|T|C|G]+$/) {
    pop @samplename;
    $SampleMergeName = join("_",@samplename);
  }
  my ($rawclust,$pfclust,$passfilt,$percq30,$avgqual_score,$yield,$permatch,$lane);
  foreach my $barcode ($sample->findnodes('./Barcode')) {
    my $barcodeName = $barcode->findvalue('@name');
    if ($barcodeName eq 'all') {
      foreach my $lane ($barcode->findnodes('./Lane')) {
	$lanename = $lane->findvalue('@number');
	my @raw_clct = map {$_->to_literal();} $lane->findnodes('./Tile/Raw/ClusterCount');
	$rawclust = sum(@raw_clct);
	my @pf_clct = map {$_->to_literal();} $lane->findnodes('./Tile/Pf/ClusterCount');
	$pfclust = sum(@pf_clct);
	my @raw_yield = map {$_->to_literal();} $lane->findnodes('./Tile/Raw/Read/Yield');
	$rawbp = sum(@raw_yield);
	my @pf_yield30 = map {$_->to_literal();} $lane->findnodes('./Tile/Raw/Read/YieldQ30');
	$q30 = sum(@pf_yield30);
	my @QualityScoreSum = map {$_->to_literal();} $lane->findnodes('./Tile/Raw/Read/QualityScoreSum');
	$qual = sum(@QualityScoreSum);
	$info{$SampleMergeName}{$lanename}{totclust} += $rawclust;
	$info{$SampleMergeName}{$lanename}{pfclust} += $pfclust;
	$info{$SampleMergeName}{$lanename}{yield} += $rawbp;
	$info{$SampleMergeName}{$lanename}{q30} += $q30;
	$info{$SampleMergeName}{$lanename}{totalquality} += $qual;
	$info{$SampleMergeName}{$lanename}{totalsamples} ++;
	$total{$lanename} += $rawclust;
      }
    }
  }
}
open OUT, ">$dir\/demultiplexing_stats.txt" or die $!;
#open OUT, ">$runname\.demultiplexstats.txt" or die $!;
print OUT join("\t","RunName","SampleName","Lane","Clusters",
	       "PFClusters","Percent.PFClusters",
	       "Yield","PercQ30","AvgQualScore","PercentLane"),"\n";

foreach $sname (keys %info) {
  next if ($sname =~ m/Undetermined/);
  foreach $lane (keys %{$info{$sname}}) {
    $clust = $info{$sname}{$lane}{totclust};
    $pf = $info{$sname}{$lane}{pfclust};
    $yield = $info{$sname}{$lane}{yield};
    print OUT join("\t",$runname,$sname,$lane,$clust,$pf,sprintf("%.2f",100*$pf/$clust),
		   $yield,sprintf("%.2f",100*$info{$sname}{$lane}{q30}/$yield),
		   sprintf("%.2f",$info{$sname}{$lane}{totalquality}/$yield),
		   sprintf("%.2f",100*$clust/$total{$lane})),"\n";
  }
}

sub sum{
  my @vals = @_;
  my $sum = 0;
  foreach $val (@vals) {
    $sum +=$val;
  }
  return $sum;
}
