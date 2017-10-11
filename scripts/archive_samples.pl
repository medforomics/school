#!/usr/bin/perl -w
#archive_samples.pl

my ($sec,$min,$hour,$mday,$mon,$year,
    $wday,$yday,$isdst) = localtime(time);

my $outdir = '/project/PHG/IR_Archival/'.join("",$year+1900,sprintf("%02s",$mon));
system("mkdir $outdir") unless (-e $outdir);

my @directories = `ls /project/PHG/PHG_Clinical/complete`;
chomp(@directories);
foreach $subjid (@directories) {
    system("tar chf $outdir\/$subjid\.tar /project/PHG/PHG_Clinical/complete/$subjid");
    system("gzip $outdir\/$subjid\.tar");
    system("rsync -L -avz /project/PHG/PHG_Clinical/complete/$subjid /project/PHG/PHG_Clinical/toarchive");
}
