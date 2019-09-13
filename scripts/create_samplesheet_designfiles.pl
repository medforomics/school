#!/usr/bin/perl -w
#run_casava.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %opt = ();
my $results = GetOptions (\%opt,'help|h','input|i=s','output|o=s','umi|u=s','dout|d=s','prjid|p=s','outdir|f=s','outnf|n=s');

if ($opt{umi}) {
  open RNASS, ">$opt{umi}" or die $!;
}
open SS, "<$opt{input}" or die $!;
open SSOUT, ">$opt{output}" or die $!;
my %sampleinfo;
my %stype;
my %spairs;
my $prjid = $opt{prjid};
my $outdir = $opt{outdir};
my $outnf = $opt{outnf};
chomp $prjid;
while (my $line = <SS>){
  chomp($line);
  $line =~ s/\r//g;
  $line =~ s/ //g;
  $line =~ s/,+$//g;
  if ($line =~ m/^\[Data\]/) {
    if ($opt{umi}) {
      print SSOUT join("\n","[Settings]","ReverseComplement,0","Read2UMILength,8"),"\n";
    }
    if ($opt{umi}) {
      print RNASS $line,"\n";
    }
    print SSOUT $line,"\n";
    $header = <SS>;
    $header =~ s/\r//g;
    chomp($header);
    $header =~ s/Sample_*/Sample_/g;
    print SSOUT $header,"\n";
    if ($opt{umi}) {
      print RNASS $header,"\n";
    }    
    my @colnames = split(/,/,$header);
    while (my $line = <SS>) {
      chomp($line);
      $line =~ s/\r//g;
      $line =~ s/ //g;
      $line =~ s/,+$//g;
      my @row = split(/,/,$line);
      my %hash;
      foreach my $j (0..$#row) {
	$hash{$colnames[$j]} = $row[$j];
      }
      $hash{Sample_Project} = $hash{Project} if $hash{Project};
      $hash{Sample_Project} =~ s/\s*$//g;
      $hash{Assay} = lc($hash{Assay});
 
      $hash{Assay} = 'panel1385' if ($hash{Assay} eq 'dnaseqdevelopment');
      $hash{Assay} = 'panel1385v2' if ($hash{MergeName} =~ m/panel1385v2/);
      $hash{Assay} = 'idthemev2' if ($hash{MergeName} =~ m/IDTHemev2/);
      $hash{Assay} = 'panelrnaseq' if ($hash{MergeName} =~ m/panelrnaseq/);
      $hash{Assay} = 'wholernaseq' if ($hash{MergeName} =~ m/wholernaseq/);
      $hash{Assay} = 'solid477' if ($hash{MergeName} =~ m/solid477/i);
      $hash{Assay} = 'pancancer1505' if ($hash{MergeName} =~ m/pancancer1505/i);
      $hash{Assay} = 'panelrnaseq1527' if ($hash{MergeName} =~ m/panelrnaseq1527/i);
      $hash{Assay} = 'heme183' if ($hash{MergeName} =~ m/heme183/i);

      unless (-e "$opt{dout}/$hash{Assay}") {
	system(qq{mkdir $opt{dout}/$hash{Assay}});
      }
      my @samplename = split(/_/,$hash{Sample_Name});
      unless ($hash{Class}) {
	$hash{Class} = 'tumor';
	$hash{Class} = 'normal' if ($hash{Sample_Name} =~ m/_N_/);
      }
      $hash{SubjectID} = $hash{Sample_Project};
      unless ($hash{MergeName}) {
	$hash{MergeName} = $hash{Sample_Name};
	if ($samplename[-1] =~ m/^[A|B|C|D]$/) {
	  pop @samplename;
	  $hash{MergeName} = join("_",@samplename);
	}
      }
      if ($hash{MergeName} =~ m/T_RNA_panelrnaseq-\d+-\d+/) {
	$hash{MergeName} =~ s/T_RNA_panelrnaseq-\d+-\d+/T_RNA_panelrnaseq/;
      }
      if ($hash{Sample_Name} =~ m/T_RNA_panelrnaseq-\d+-\d+/) {
	$hash{Sample_Name} =~ s/T_RNA_panelrnaseq(-\d+-\d+)/T_RNA_panelrnaseq_Lib.$1/e;
      } 
      my $clinres = 'cases';
      $hash{VcfID} = $hash{SubjectID}."_".$prjid;
      if (($hash{Description} && $hash{Description} =~ m/research/i) ||
	  ($hash{Sample_Name} !~ m/ORD/ && $hash{SubjectID} !~ m/GM12878|ROS/)) {
	$clinres = 'researchCases';
      }
      $hash{ClinRes} = $clinres;
      unless($opt{umi}){ #unless ($umi) {
	$hash{Sample_Name} = $hash{Sample_Name}."_ClarityID-".$hash{Sample_ID};
      }
      $hash{Sample_ID} = $hash{Sample_Name};
      $stype{$hash{SubjectID}} = $hash{Case};
      $spairs{$hash{Assay}}{$hash{SubjectID}}{lc($hash{Class})}{$hash{MergeName}} = 1 unless ($hash{Assay} =~ m/rna/);
      $sampleinfo{$hash{Sample_Name}} = \%hash;
      push @{$samples{$hash{Assay}}{$hash{SubjectID}}}, $hash{Sample_Name};
      
      my @newline;
      foreach my $j (0..$#row) {
	push @newline, $hash{$colnames[$j]};
      }
      if ($opt{umi} && $hash{Index_Name} !~ m/UMI/) {
	print RNASS join(",",@newline),"\n";
      }else {
	print SSOUT join(",",@newline),"\n";
      }
    }
  } else {
    if ($opt{umi}) {
      print RNASS $line,"\n";
    }
    print SSOUT $line,"\n";
  }
}
close SSOUT;
my %inpair;
foreach $dtype (keys %spairs ){#%stype) {
  open TNPAIR, ">$opt{dout}/$dtype/design_tumor_normal.txt" or die $!;
  print TNPAIR join("\t",'PairID','VcfID','TumorID','NormalID','TumorBAM','NormalBAM',
		    'TumorGATKBAM','NormalGATKBAM'),"\n";
  foreach my $subjid (keys %{$spairs{$dtype}}) {
    my @ctypes = keys %{$spairs{$dtype}{$subjid}};
    if ($spairs{$dtype}{$subjid}{tumor} && $spairs{$dtype}{$subjid}{normal}) {
      my @tumors = keys %{$spairs{$dtype}{$subjid}{tumor}};
      my @norms = keys %{$spairs{$dtype}{$subjid}{normal}};
      my $pct = 0;
      foreach $tid (@tumors) {
	$inpair{$tid} = 1;
	foreach $nid (@norms) {
	  $inpair{$nid} = 1;
		    my $pair_id = $subjid;
	  if ($pct > 0) {
	    $pair_id .= ".$pct";
	  }
	  print TNPAIR join("\t",$pair_id,$pair_id."_".$prjid,$tid,$nid,$tid.".consensus.bam",
			    $nid.".consensus.bam",$tid.".final.bam",$nid.".final.bam"),"\n";
	  $pct ++;
	}
      }
    }
  }
  close TNPAIR;
}
foreach $dtype (keys %samples) {
  open CAS, ">$opt{dout}\/$dtype\/lnfq.sh" or die $!;
  print CAS "#!/bin/bash\n";
  open SSOUT, ">$opt{dout}\/$dtype\/design.txt" or die $!;
  open TONLY, ">$opt{dout}\/$dtype\/design_tumor_only.txt" or die $!;
  print SSOUT join("\t","SampleID",'FamilyID','FqR1','FqR2'),"\n";
  print TONLY join("\t","SampleID",'VcfID','FamilyID','BAM','GATKBAM'),"\n";
  my %thash;
  foreach $project (keys %{$samples{$dtype}}) {
    my $datadir =  "/project/PHG/PHG_Clinical/illumina/$prjid/$project/";
    foreach $samp (@{$samples{$dtype}{$project}}) {
      my %info = %{$sampleinfo{$samp}};
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $outdir\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $outdir\/$samp\.R2.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R1_*.fastq.gz $outnf\/$info{SubjectID}/fastq\/$samp\.R1.fastq.gz\n";
      print CAS "ln -s $datadir/$samp*_R2_*.fastq.gz $outnf\/$info{SubjectID}/fastq\/$samp\.R2.fastq.gz\n";
      print SSOUT join("\t",$info{MergeName},$info{SubjectID},
		       "$samp\.R1.fastq.gz","$samp\.R2.fastq.gz"),"\n"; 
      next if ($inpair{$info{MergeName}});
      if ($dtype =~ m/rna/) {
	print TONLY join("\t",$info{MergeName},$info{VcfID},$info{SubjectID},
			 $info{MergeName}.".bam",$info{MergeName}.".bam"),"\n";
      }else {
	print TONLY join("\t",$info{MergeName},$info{VcfID},$info{SubjectID},
			 $info{MergeName}.".consensus.bam",$info{MergeName}.".final.bam"),"\n";
      }
    }
  }
  close SSOUT;
  close TONLY;
  close CAS;
}
