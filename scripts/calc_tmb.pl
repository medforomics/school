#!/usr/bin/perl -w
#calc_tmb.pl

my $prefix = shift @ARGV;
my $infile = shift @ARGV;

my $cdssize = 4.650646;

open IN, "<$infile" or die $!;
open TMB, ">$prefix\.tumorburden.txt" or die $!;
my %varct;
while (my $line = <IN>) {
    chomp($line);
    next unless ($line =~ m/^SN/);
    my ($sn,$id,$key,$value) = split(/\t/,$line);
    $key =~ s/number of //;
    $key =~ s/://;
    next if ($key =~ m/sites|samples/i);
    if ($key =~ m/records/) {
	$varct{'total'} = $value;
    }else {
	$varct{lc($key)} = $value;
    }
    
}
print TMB join("\t","Sample","MutationType","MutationCt","TumorMutationBurden"),"\n";
foreach $t (sort {$a cmp $b} keys %varct) {
    print TMB join("\t",$prefix,$t,$varct{$t},sprintf("%.1f",$varct{$t}/$cdssize)),"\n";
}
close TMB;
