#!/usr/bin/perl -w
#rerun_qc.sh

my $bam = shift @ARGV;
my $flagstat = shift @ARGV;
my $capture = shift @ARGV;

my $rnum = 4e6;
if ($capture =~ m/exome/i) {
    $rnum = 5e7;
}

$prefix = (split(/\.ontarget.bam/,$bam))[0];
open FLAG, "<$flagstat" or die $!;
my ($total, $read1ct,$read2ct,$maprate,$concorrate);
while (my $line = <FLAG>) {
    chomp($line);
    if ($line =~ m/(\d+) \+ \d+ in total/) {
	$hash{total} = $1 unless $total;
    }elsif ($line =~ m/(\d+) \+ \d+ read1/) {
	$hash{pairs} = $1;
    }elsif ($line =~ m/(\d+) \+ \d+ mapped/) {
	$hash{maprate} = 100*sprintf("%.4f",$1/$hash{total});
    }elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
	$hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
    }elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
	$hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
    }
}
$percreads = sprintf("%.4f",$rnum/$hash{total});
system(qq{sambamba view -t 30 -f bam -s $percreads -o $prefix\.subset.bam $bam});
