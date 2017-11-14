#!/usr/bin/perl
#philipsCSV.pl

use strict;
use warnings;

my $pass_38 = $ARGV[0];
my $pass_37 = $ARGV[1];

chomp $pass_38;
chomp $pass_37;
open CSV38, "<$pass_38" or die $!;
open CSV37, "<$pass_37" or die $!;

my $header_38 = <CSV38>;
my $header_37 = <CSV37>;

chomp $header_38;
my ($header_chrpos, @header_info) = split(",", $header_38);

my %hash_38;
my %hash_37;
while(my $line = <CSV38>){
  chomp $line;
  my ($chr_pos,$gene,$aa,$effect,@info) = split(",", $line);
  my $info_38 = $gene.$aa.$effect;
  $hash_38{$info_38} = $chr_pos;
}
my $final_csv = $pass_38;
chomp $final_csv;
$final_csv =~ s/\.PASS/\.PASS\.38\.37/;
open OUT, ">$final_csv" or die $!;
print OUT join(",","chrpos_GRCh38","chrpos",@header_info),"\n";
while(my $line = <CSV37>){
  chomp $line;
  my ($chr_pos,$gene,$aa,$effect,@info) = split(",", $line);
  my $info_37 = $gene.$aa.$effect;
  my $final_line="";
  if($hash_38{$info_37}){
    
    print OUT join(",",$hash_38{$info_37},$chr_pos,$gene,$aa,$effect,@info),"\n";
  }
  else{
    print OUT join(",","",$chr_pos,$gene,$aa,$effect,@info),"\n";
  }
}


