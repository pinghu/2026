#!/usr/bin/perl -w
use strict;
my $usage="$0 file long syn\n";
die $usage unless @ARGV >=1;
my ($file, $longF, $synF) = @ARGV;
if(! defined $longF){

    $longF="/media/ping/Disk2/project/GSS3049Petroleum/gene.name";
}
if(! defined $synF){

    $synF="/media/ping/Disk2/project/GSS3049Petroleum/gene.EZ";
}
my %ann;
fill_hash(\%ann, $longF);

my %acr;
fill_hash(\%acr, $synF);

#link = 'http://iip.pg.com:10080/iip-web/geneSummary.do?&annotationdata=20,56,13,74,44,11,55,50,51,49&acronym=' + gene;
open (C, "<$file")|| die "could not open $file\n";
open (B, ">$file.xls")|| die "could not open $file.xls\n";
while (my $M =<C>){
    chomp $M;
    my @tmp=split /\t/, $M;
    my $long="name";
    my $syn="EZID";
    my $probe=$tmp[0];
    if($probe =~ /(.*)_\d+$/){
	$probe=$1;
    }
    my $p=$probe;
    if (defined $ann{lc($probe)}){
	$long=$ann{lc($probe)};
    }
    if (defined $acr{lc($probe)}){
	$syn=$acr{lc($probe)};
    }
    my $x=$syn;
    
    print B "$tmp[0]\t$probe\t$x\t$long";
    for(my $i=1; $i<@tmp; $i++){
	print B "\t", $tmp[$i];
    }
    print B "\n";
}
close C;
close B;

sub fill_hash{
    my ($h, $file)=@_;
    open (X, "$file")||die "could not open $file\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (@tmp<2){next;}
	$$h{lc($tmp[0])}=lc($tmp[1]);
    }
    close X;
    return;
}
