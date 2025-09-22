#!/usr/bin/perl
use strict;
use warnings;

$|=1;

if(@ARGV!=7){
	print "perl $0 <bc_umi_read.tsv.id> <BSTViewer_project/subdata/L1_heAuto/barcodes.tsv.gz> <sample/barcode.tsv> <BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz> <in.bam> <out.bam> <L1>\n";
	exit;
}
print "$ARGV[4]\n";
print "$ARGV[5]\n";
open(IN,"$ARGV[0]");
my %hash;
while(<IN>){
	chomp;
	my ($bc,$umi,@reads)=split(/\t/,$_);
	#next unless($bc=~/^bc1_\d+-bc2_\d+-bc3_\d+/);
	foreach my $read(@reads){
		$hash{$read}{'bc'}=$bc;
		$hash{$read}{'umi'}=$umi;
	}
}
print "Loaded bc_umi_read.tsv.id\n";

open(BC,"gzip -cd $ARGV[1]|") or die $!;
my %sub_bc;
while(<BC>){
	chomp;
	$_=~s/.*_//;
	$sub_bc{$_}=1;
	#print $_,"\n";
}
print "Loaded barcodes.tsv.gz\n";

open(B,"$ARGV[2]");
my %line2bc;
while(<B>){
	chomp;
	$line2bc{$.}=$_;
}
print "Loaded barcode.tsv\n";

open(C,"gzip -cd $ARGV[3]|");
my %bc;
my %bc2clusterbc;
while(<C>){
	chomp;
	my @a=split(/\t/,$_);
	next if(not exists $sub_bc{$a[1]-1});
	$bc{$line2bc{$a[0]}}=1;
	#$bc2clusterbc{$line2bc{$a[0]}}="L7_".($a[1]-1);
	$bc2clusterbc{$line2bc{$a[0]}}=$ARGV[6]."_".($a[1]-1);
	#print $line2bc{$a[0]},"\n";
}
print "Loaded barcodes_cluster.tsv.gz\n";

open(OUT,"|samtools view -Sb > $ARGV[5]") or die "Cannot open samtools OUT pipe: $!";
print "Opened samtools OUT\n";

open(BAM,"samtools view -h $ARGV[4]|") or die "Cannot open samtools BAM pipe: $!";
print "Opened samtools BAM\n";

while(<BAM>){
	chomp;
	if(/^@/){
		print OUT "$_\n";
	}else{
		my @a=split(/\t/,$_);
		next if(not exists $hash{$a[0]});
		next if(not exists $bc{$hash{$a[0]}{'bc'}});
		#$hash{$a[0]}{'bc'}=~s/-/_/g;
		my $bc=$bc2clusterbc{$hash{$a[0]}{'bc'}};
		my $umi=$hash{$a[0]}{'bc'}."_"."$hash{$a[0]}{'umi'}";
		#print OUT "$_\tCR:Z:$hash{$a[0]}{'bc'}\tCB:Z:$hash{$a[0]}{'bc'}-1\tUR:Z:$hash{$a[0]}{'umi'}\tUB:Z:$hash{$a[0]}{'umi'}\n";
		print OUT "$_\tCR:Z:$bc\tCB:Z:$bc-1\tUR:Z:$umi\tUB:Z:$umi\n";
	}

	
}
