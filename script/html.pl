#!/usr/bin/perl

use strict;
use warnings;
use autodie;

open my $fi,'<',$ARGV[0] or die "Cannot open file: $!";

$/="|";
<$fi>;

while(<$fi>){
	chomp;
	$_=~ s/\n\,\n//g;
	$_=~ s/\,//g;
	my ($region,$type,$from,$to,$similar_type,$similar_compound,$similarity)=(split/\n/,$_)[0,1,2,3,4,5,6];
	$region=~ s/[\s]*//g;
	$type=~ s/[\s]*//g;
	$from=~ s/[\s]*//g;
	$to=~ s/[\s]*//g;
	$similar_type=~ s/[\s]*//g;
	$similar_compound=~ s/[\s]*//g;
	$similarity=~ s/[\s]*//g;
	print "$region\t$type\t$from\t$to\t$similar_type\t$similar_compound\t$similarity\n";
}

$\="|";
close $fi;