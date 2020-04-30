#!/usr/bin/perl -w
#update_clarity_xml.pl

my %info;
my %udf;
my @xmlfiles = @ARGV;
foreach $xmlfile (@xmlfiles) {
    chomp($xmlfile);
    open IN, "<$xmlfile" or die $!;
    my $limsid;
    my ($dirname,@ext) = split(/\./,$xmlfile);
    while (my $line = <IN>) {
	chomp($line);
	if ($line =~ m/limsid="(\w+)"/) {
	    $limsid = $1;
	    $info{$limsid}{'projectname'} = $dirname;
	}
	if ($line =~ m/<udf:field type="String" name="(.+)">(.+)<\/udf:field>/) {
	    my $name = $1;
	    my $val = $2;
	    $info{$limsid}{$name} = $val;
	    $udf{$name} = 1;
	}
    }
    close IN;
}

open OUT, ">clarity_info.txt" or die $!;
my @names = sort {$a cmp $b} keys %udf;
print OUT join("\t","LIMSID","PRJID",@names),"\n";
foreach $limsid (sort {$a cmp $b} keys %info) {
    my @line = ();
    foreach $udf (@names) {
	$info{$limsid}{$udf} = '' unless ($info{$limsid}{$udf});
	push @line, $info{$limsid}{$udf};
    }
    print OUT join("\t",$limsid,$info{$limsid}{'projectname'},@line),"\n";
}
