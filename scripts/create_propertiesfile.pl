#!/usr/bin/perl -w
#run_casava.pl

my $casedir = shift @ARGV;
my @dir = split(/\//,$casedir);
my $caseid=pop @dir;
my $seqrunid=pop @dir;


my @prop = ('run.name','dmux.conversion.stats','bait.pool',
	    'project.name','sample.alignment','sample.coverage.raw',
	    'sample.name','somatic.seq.stats','giab.snsp',
	    'somatic.translocation');


my $load_root = "/swnas/qc_nuclia";
my $demux= "$load_root/demultiplexing/$seqrunid";
my $seqqc = "$load_root/seqanalysis/$seqrunid";

my $proc = join("/",@dir,$seqrunid);
my @files = `ls $proc/*.design.txt`;
chomp(@files);
my %dtype;
my %pair;
my @samples;
foreach $f (@files) {
  open D, "<$f" or die $!;
  my $header = <D>;
  chomp($header);
  my @colnames = split(/\t/,$header);
  while (my $line = <D>) {
    chomp($line);
    my @row=split(/\t/,$line);
    my %hash;
    foreach my $j (0..$#row) {
      $hash{$colnames[$j]} = $row[$j];
    }
    if ($hash{CaseID} eq $caseid) {
      if ($hash{NormalID}) {
	$pair{TumorID} = $hash{TumorID};
	$pair{NormalID} = $hash{NormalID};
	push @samples, ($hash{TumorID}, $hash{NormalID});
      } else {
	push @samples, $hash{SampleID};
      }
    }
    foreach $id (grep/ID/, @colnames) {
      $dtype{$hash{$id}} = $hash{PanelFile};
    }
  }
}

$somatic = 0;
$somatic = 1 if ($pair{NormalID});

foreach $sampleid (@samples) {
  my %info;
  my $nasdir="$seqqc/$caseid";
  $info{'dmux.conversion.stats'} = "$demux/Stats/ConversionStats.xml";
  $info{'run.name'} = $seqrunid;
  $info{'project.name'}=$caseid;
  $info{'sample.name'}=$sampleid;
  $info{'bait.pool'} = $dtype{$sampleid};
  my $outdir;
  if ($sampleid =~ m/RNA/) {
    $outdir="rnaout";
  }else {
    $outdir="dnaout";
    if ($somatic) {
      system(qq{perl -pi -e 's/Sample_1.+/Sample_1\t$pair{TumorID}/g' $casedir/dnacallset/$caseid\.sequence.stats.txt});
      system(qq{perl -pi -e 's/Sample_2.+/Sample_2\t$pair{NormalID}/g' $casedir/dnacallset/$caseid\.sequence.stats.txt});
      $info{'somatic.seq.stats'}="$nasdir/dnacallset/$caseid\.sequence.stats.txt";
    }
  }
  $info{'sample.alignment'} = "$nasdir/$outdir/$sampleid.sequence.stats.txt";
  $info{'sample.coverage.raw'} = "$nasdir/$outdir/$sampleid\_exoncoverage.txt";
  if ($caseid =~ m/ROS1/) {
    $info{'somatic.translocation'} = "$nasdir/rnaout/$sampleid.translocations.answer.txt";
  }
  if ($caseid =~ m/GM12878/) {
    $info{'giab.snsp'} = "$nasdir/$sampleid.snsp.txt";
  }
  open OUT, ">$casedir/$sampleid\.properties" or die $!;
  foreach $prop (@prop) {
    $info{$prop} = '' unless ($info{$prop});
    print OUT join("=",$prop,$info{$prop}),"\n";
  }
}
