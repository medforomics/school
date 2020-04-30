#!/usr/bin/perl -w
#update_clarity_xml.pl

my @xmlfiles = @ARGV;
foreach $xmlfile (@xmlfiles) {
    chomp($xmlfile);
    my $dir = $xmlfile;
    $dir =~ s/\.clarity.xml//;
    my @dnadir = `ls -d $dir\/*DNA_panel1385`;
    my @rnadir = `ls -d $dir\/*panelrnaseq`;
    chomp(@dnadir);
    chomp(@rnadir);
    my @tids = grep(/T_DNA/,@dnadir);
    my @nids = grep(/N_DNA/,@dnadir);
    my $rid = shift @rnadir;
    my $tumorid = (split(/\//,$tids[0]))[-1];
    my $normalid = (split(/\//,$nids[0]))[-1];
    my $rnaid  = (split(/\//,$rid))[-1];
    my $outfile = (split(/\//,$xmlfile))[-1];
    
    open IN, "<$xmlfile" or die $!;
    open OUT, ">$outfile" or die $!;
    
    while (my $line = <IN>) {
	chomp($line);
	if ($line =~ m/<\/prj:project>/) {
	    if ($tumorid) {
		print OUT qq{    <udf:field type="String" name="Tumor DNA primary sample">$tumorid</udf:field>},"\n";
	    }if ($normalid) {
		print OUT qq{    <udf:field type="String" name="Normal DNA primary sample">$normalid</udf:field>},"\n";
	    }if ($rnaid) {
		print OUT qq{    <udf:field type="String" name="Tumor RNA primary sample">$rnaid</udf:field>},"\n";
	    }
	}
	print OUT $line,"\n";
    }
    close OUT;
    close IN;
}
