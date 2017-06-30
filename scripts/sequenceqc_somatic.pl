#!/usr/bin/perl -w
#uploadqc.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'refdir|r=s','input|i=s','output|o=s','help|h');

open MATE, "<$opt{input}" or die $!;

### Begin File Information ###
my @stats = stat("$opt{input}");
my ($day,$month,$year) = (localtime($stats[9]))[3,4,5];
$year += 1900;
$month++;
$date = join("-",$year,sprintf("%02s",$month),sprintf("%02s",$day));
$fileowner = 's'.$stats[4];
$fileowner = $fileowner;
### End File Information ###

### Begin Version Information ###
my $source = `zgrep '#' $opt{refdir}\/cosmic.vcf.gz |grep source`;
my $cosmic_ref = (split(/=/, $source))[1];
chomp($cosmic_ref);
my $dbsnp_source = `zgrep '#' $opt{refdir}\/dbSnp.vcf.gz |grep dbSNP_BUILD_ID`;
my $dbsnp_ref = (split(/=/, $dbsnp_source))[1];
chomp($dbsnp_ref);
my $gen_ref = (split(/\//,$opt{refdir}))[-2];
my $gittag = `cd /project/PHG/PHG_Clinical/clinseq_workflows;git describe --tag`;
chomp $gittag;
### End Version Information ###

while (my $line = <MATE>) {
    chomp($line);
    my ($sam1,$pf,$sam2,$corr,$depth) = split(/\t/,$line);
    $sam1 =~ s/\.vcf//g;
    $sam2 =~ s/\.vcf//g;
    open OUT, ">$opt{output}" or die $!;
    my $status= 'PASS';	
    $status='FAIL' if($pf eq 'unmatched');
    print OUT join("\n","Sample_1\t".$sam1,"Sample_2\t".$sam2,"Correlation\t".$corr,"Depth\t".$depth,"Status\t".$status, 
		   "Somatic_Date\t".$date,"File_Owner\t".$fileowner,"Workflow_Version\t".$gittag,"Cosmic_Reference\t".$cosmic_ref,
		   "dbSnp_Reference\t".$dbsnp_ref,"Genome_Reference\t".$gen_ref),"\n";
}



