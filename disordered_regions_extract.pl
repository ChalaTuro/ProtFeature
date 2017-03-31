#! usr/bin/perl

## Script used for extracting disordered regions from disemble output, which submitted for cprofile of disordered regions separately.

##Chala Jefuka Turo

## Usage: script.pl  inputfile 

## print out disordered regions in a fastafile format (not concatenated)


open(FAS,'<', @ARGV[0]);
open(OUT1, '>', @ARGV[0] . "." . "COILS");
open(OUT2, '>', @ARGV[0] . "." . "REM465");
open(OUT3, '>', @ARGV[0] . "." . "HOTLOOPS");

$fasta = do {local $/; <FAS>};
@fasarray = split(/>/, $fasta);
%seqcount = undef;
%fashash = undef;
%occurence = undef;
$length = 0;
$threshold = 5;
$count = 0;
@alphabet = (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z,X);
#print "protein_id\tdisordered_type\taa_lengt\tdisorder_length";

foreach $aa (@alphabet) {
#	print  "\t$aa";
}

#print "\n";

foreach $fasta (@fasarray) {	
	@data = split(/\n/, $fasta, 2);
	@data[1] =~ s/\n|\t|\s+//sgi;  # seq	
	#print "@data[0]\n";  # id plus scores
	#print "@data[1]\n";   #seq
	#print "$total_length\n";
	@idarray = split(/\_/, @data[0]);
	#print "@idarray[0]\n";  ### id1     holds string "proteinid"
	#print "@idarray[1]\n"; ## id2      holds actual ids
	#print "@idarray[2]\n";  ### disorder_type_plus_scores
	$fasid1 = @idarray[0]; # for Paula data only, otherwise use next line only
	$fasid = @idarray[1];
	$seq = @data[1];
	
	#$seq =~ tr/abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/;
	$seq =~ tr/arndcqeghilkmfpstwyvbzxARNDCQEGHILKMFPSTWYVBZX/ARNDCQEGHILKMFPSTWYVBZXARNDCQEGHILKMFPSTWYVBZX/; ## aa from pfilt software
	$total_length = length($seq);
	@scores = split ( /\s+/, @idarray[2], 2);
	##print "@scores[0]\n";   # disorder classes
	##print "@scores[1]\n";  # all the range of scores				
	@val = (split /,\s+/, @scores[1]);
	$regions_disorder = undef;
	$disorder_length = undef;
	$full_ordered_length = undef;						
	foreach $range (@val) {
		unless ($range eq undef) {
			@test = split (/\-/, $range);
			$length = @test[1] - @test[0] + 1; 
			if ($length >= 5){    		#restrict the minimum disordered region to extract; 0 means no minimum
				$disordersubstr = substr ($seq, @test[0], $length -1 ); #all substers	
				$count++;			
			}
			 if (@scores[0] =~ "COILS") {
			print OUT1 ">protid_" . $fasid . "_" . @scores[0] . ".$count" . "\n" . $disordersubstr ."\n";
			}
			elsif (@scores[0] =~ "REM465") {
			print OUT2 ">protid_" . $fasid . "_" . @scores[0] . ".$count" . "\n" . $disordersubstr ."\n";		
			} 
			else { 
			print OUT3 ">protid_" . $fasid . "_" . @scores[0] . ".$count" . "\n" . $disordersubstr ."\n";
			}
																	
		}
#	print ">protid_" . $fasid . "_" . @scores[0] . ".$count" . "\n" . $disordersubstr ."\n";	
	}
}
 close FAS; 	
exit;				
			
