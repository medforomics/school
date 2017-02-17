#!/usr/bin/perl -w
#parse_conversion_stats.pl

use XML::LibXML;

my $dir = shift @ARGV;

open OUT, ">$dir\/demultiplexing_stats.txt" or die $!;
print OUT join("\t","RunName","SampleName","Clusters",
	       "PFClusters","Percent.PFClusters",
	       "Yield","PercQ30","AvgQualScore"),"\n";

my $xml = "$dir\/Stats/ConversionStats.xml";
$runname = (split(/\//,$xml))[-1];
my $parser = XML::LibXML->new();
my $doc = $parser->parse_file($xml);
foreach my $sample  ($doc->findnodes('Stats/Flowcell/Project/Sample')) {
  my $sampleName = $sample->findvalue('@name');
  next if ($sampleName eq 'all' || $sampleName eq 'Undetermined');
  my ($rawclust,$pfclust,$passfilt,$percq30,$avgqual_score,$yield,$permatch);
  foreach my $barcode ($sample->findnodes('./Barcode')) {
    my $barcodeName = $barcode->findvalue('@name');
    if ($barcodeName eq 'all') {
      my @raw_clct = map {$_->to_literal();} $barcode->findnodes('./Lane/Tile/Raw/ClusterCount');
      $rawclust = sum(@raw_clct);
      my @pf_clct = map {$_->to_literal();} $barcode->findnodes('./Lane/Tile/Pf/ClusterCount');
      $pfclust = sum(@pf_clct);
      $perc_pfclust = sprintf("%.0f",100*$pfclust/$rawclust);
      
      my @raw_yield = map {$_->to_literal();} $barcode->findnodes('./Lane/Tile/Raw/Read/Yield');
      $rawbp = sum(@raw_yield);
      my @pf_yield30 = map {$_->to_literal();} $barcode->findnodes('./Lane/Tile/Raw/Read/YieldQ30');
      $percq30 = sprintf("%.2f",100*sum(@pf_yield30)/$rawbp);
      my @QualityScoreSum = map {$_->to_literal();} $barcode->findnodes('./Lane/Tile/Raw/Read/QualityScoreSum');
      $avgqual_score = sprintf("%.2f",sum(@QualityScoreSum)/$rawbp);
    }else {
      my @all_raw_clct = map {$_->to_literal();} $barcode->findnodes('./Lane/Tile/Raw/ClusterCount');
      $permatch = sum(@all_raw_clct);
    }
  }
  my $perfect_barcode = sprintf("%.0f",100*$permatch/$rawclust);
  print OUT join("\t",$runname,$sampleName,$rawclust,$pfclust,
		 $perc_pfclust,sprintf("%.2e",$rawbp),$percq30,$avgqual_score),"\n";
}

sub sum{
  my @vals = @_;
  my $sum = 0;
  foreach $val (@vals) {
    $sum +=$val;
  }
  return $sum;
}
