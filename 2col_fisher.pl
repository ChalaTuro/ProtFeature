#!/usr/bin/perl
#usage run_fisher_test input output
#input 2 column table of go comp... obs1, obs2 (first line equals totals for cols 1 and 2)
#output --> STDOUT: col1, col2, fisherpvalue
use Math::Round;
open(INPUT, '<', @ARGV[0]);

###################FACTORIAL CALCULATION#############################
#print "start factorial calculation...\n";
$n = 1000000;
@fact = undef;
@fact[0] = 1;
for ($i = 2; $i < ($n+1); $i++) {
	@fact[$i] = @fact[$i-1] + log($i);
}
#print "finished factorial calculation\n";
#####################################################################

$line = <INPUT>;
chomp $line;
@tabdata = split(/\t/, $line);
$tot1 = round(@tabdata[0]);
$tot2 = round(@tabdata[1]);

while (<INPUT>) {
	$line = $_;
	chomp $line;
	@tabdata = split(/\t/, $line);
	
	$a = round(@tabdata[0]);
	$b = round(@tabdata[1]);
	$c = $tot1 - $a;
	$d = $tot2 - $b;
	$n = $a + $b + $c + $d;
	unless(@tabdata[0] eq "GO_id") {
		$hyp1 = @fact[($a+$b)] + @fact[($c+$d)] + @fact[($a+$c)] + @fact[($b+$d)];
		$hyp2 = @fact[$n] + @fact[$a] + @fact[$b] + @fact[$c] + @fact[$d];
		$hyp1 = $hyp1 - $hyp2;
		$hyp1 = exp($hyp1);
	}
	#if ($hyp1 < 0.01) {
		print "$a\t$b\t$hyp1\n";
	#}
}
close FILEOUT;

