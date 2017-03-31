#! usr/bin/perl

## Script used for amino acid frequency analysis  within ordered and disordered regions

## Usage: disemble fileout > fileout


open(FAS,'<', @ARGV[0]);
$fasta = do {local $/; <FAS>};
@fasarray = split(/>/, $fasta);
%seqcount = undef;
%fashash = undef;
%occurence = undef;
@headers=(Tiny_total,Tiny_ordered,Tiny_disordered,Small_total,Small_ordered,Small_disordered,Aliphaic_total,
Aliphaic_ordered,Aliphaic_disordered,Aromatic_total,Aromatic_ordered,Aromatic_disordered,Nonpolar_total,
Nonpolar_ordered,Nonpolar_disordered,Polar_total,Polar_ordered,Polar_disordered,Charged_total,Charged_ordered,
Charged_disordered,positive_charge_total,positive_charge_ordered,positive_charge_disordered,negative_charge_total,
negative_charge_ordered,negative_charge_disordered,Basic_total,Basic_ordered,Basic_disordered,Acidic_total,Acidic_ordered,
Acidic_disordered,large_total,large_ordered,large_disordered,hydrophobic_Eisenberg_total,hydrophobic_Eisenberg_ordered,
hydrophobic_Eisenberg_disordered,Hydrophobic_KD_total,Hydrophobic_KD_ordered,Hydrophobic_KD_disordered,Hydrophobic_FP_total,
Hydrophobic_FP_ordered,Hydrophobic_FP_disordered, helix_total, helix_ordered, helix_disordered,turn_total,turn_ordered,
turn_disordered,sheet_total, sheet_ordered, sheet_disordered);
#$length = 0;
#$Tiny = 0;		#(A+C+G+S+T)		From EMBOSS
#$Small = 0; 		#(A+B+C+D+G+N+P+S+T+V)	From EMBOSS 
#$Aliphatic = 0;	#(A+I+L+V)		From EMBOSS
#$Aromatic = 0;		#(F+H+W+Y)		From EMBOSS
#$Nonpolar = 0;		#(A+C+F+G+I+L+M+P+V+W+Y) From EMBOSS	
#$Polar = 0;		#(D+E+H+K+N+Q+R+S+T+Z)	 From EMBOSS
#$Charged = 0;		#(B+D+E+H+K+R+Z)	From EMBOSS
#$positive_charge= 0;
#$negative_charge= 0;		
#$Basic = 0;		#(H+K+R)			
#$Acidic	 = 0;	#(B+D+E+Z)
#$helix = 0; 		# V, I, Y, F, W, L.  From http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam.ProteinAnalysis-class.html#count_amino_acids
#$turn = 0;		# N, P, G, S.  from http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam.ProteinAnalysis-class.html#count_amino_acids
#$sheet = 0;		# E, M, A, L.	from http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam.ProteinAnalysis-class.html#count_amino_acids	
#$large = 0; 		#NETVILPQHMFKWYR   From cprofile Aminoacid.rb
#$hydrophobic_Eisenberg = 0; 	#PYCGAMWLVFI From cprofile Aminoacid.rb
#$Hydrophobic_KD = 0; 	#AMCFLVI
#$Hydrophobic_FP = 0;	#HTAPYVCLFIMW 

#$Tiny_dis = 0;		
#$Small_dis = 0; 		
#$Aliphatic_dis = 0;		
#$Aromatic_dis = 0;		
#$Nonpolar_dis = 0;		
#$Polar_dis = 0;		
#$Charged_dis = 0;
#$negative_charge = 0;		
#$positive_charge_dis= 0;
#$negative_charge_dis= 0;		
#$Basic_dis = 0;					
#$Acidic_dis	 = 0;		
#$helix_dis = 0; 		
#$turn_dis = 0;		
#$sheet_dis = 0;		
#$large_dis = 0; 		
#$hydrophobic_Eisenberg_dis = 0; 	
#$Hydrophobic_KD_dis = 0; 	
#$Hydrophobic_FP_dis = 0;
#$Tiny_ord = 0;		
#$Small_ord = 0; 		
#$Aliphatic_ord = 0;		
#$Aromatic_ord = 0;		
#$Nonpolar_ord = 0;		
#$Polar_ord = 0;		
#$Charged_ord = 0;		
#$positive_charge_ord= 0;
#$negative_charge_ord= 0;		
#$Basic_ord = 0;					
#$Acidic_ord	 = 0;		
#$helix_ord = 0; 		
#$turn_ord = 0;		
#$sheet_ord = 0;		
#$large_ord = 0; 		
#$hydrophobic_Eisenberg_ord = 0; 	
#$Hydrophobic_KD_ord = 0; 	
#$Hydrophobic_FP_ord = 0;	



#@alphabet = (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z,X);

@alphabet = (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V);
#print "protein_id\tdisordered_type\tlength_total\tlength_disorderd\tlength_ordered";

print "proteinId\tdisorderType\tlength_total\tlength_disordered\tno_of_disorder_region\tlength_ordered";
foreach $aa (@alphabet) {
	print  "\t$aa" ."_total";
	print  "\t$aa" ."_disordered";
	print  "\t$aa" ."_ordered";
}

foreach $header(@headers) {
	unless ($header eq undef) {
	print "\t$header";
	}
}

print "\n";

foreach $fasta (@fasarray) {	
	@data = split(/\n/, $fasta, 2);
	@data[1] =~ s/\n|\t|\*|\s+//sgi;  # seq	
	#print "@data[0]\n";  # id plus scores
	$no_of_disorder_region = @data[0] =~ tr/-/-/;
	#print "@data[1]\n";   #seq
	#print "$total_length\n";
	@idarray = split(/\_/, @data[0]);
	#print "@idarray[0]\n";  ### id1     holds string "proteinid"
	#print "@idarray[1]\n"; ## id2      holds actual ids
	#print "@idarray[2]\n";  ### disorder_type_plus_scores
	$fasid = @idarray[1];
	$fasid =~ s/\t|\s+//sgi;
	$seq = @data[1];
	
	$seq =~ tr/arndcqeghilkmfpstwyvbzxARNDCQEGHILKMFPSTWYVBZX/ARNDCQEGHILKMFPSTWYVBZXARNDCQEGHILKMFPSTWYVBZX/; ## aa from pfilt software
	$total_length = length($seq);
	@scores = split ( /\s+/, @idarray[2], 2);
	##print "@scores[0]\n";   # disorder classes
	##print "@scores[1]\n";  # all the range of scores				
	@val = (split /,\s+/, @scores[1]);
	$regions_disorder = undef;
	$disorder_length = undef;
	$ordered_regions = undef;
	$full_ordered_length = undef;						
	foreach $range (@val) {
		unless ($range eq undef) {
			@test = split (/\-/, $range);
			$length = @test[1] - @test[0] + 1; 
			if ($length >= 0){
				$disordersubstr = substr ($seq, @test[0], $length); #all substers
				$regions_disorder .= $disordersubstr;
				$disorder_length = length($regions_disorder);
				$full_ordered_length = $total_length - $disorder_length;
				#$full_ordered_length = length($ordered_regions);
			}
																		
		}
	}
	
			
#	print "$regions_disorder\n";
	@residues = split (//,$seq );
	%occurence = undef;
		
	foreach $aa (@residues) {
		$occurence{'full'}->{$aa}++;
		#print "$aa\t$occurence{'full'}->{$aa}\n";  ## OK
		
		
	}
	
	@residues = split (//, $regions_disorder);		
	foreach $aa (@residues) {
		$occurence{'disorder'}->{$aa}++;
		
		#print "$aa   $occurence{'disorder'}->{$aa}\n";
	}
	
	
	##print  "protid_" . $fasid . "\t" . @scores[0] . "\t" .$occurence{'disorder'}->{$t};
	
	unless ($fasid eq undef) {
	print "protid_" . $fasid . "\t" . @scores[0] . "\t" . $total_length . "\t" . $disorder_length . "\t" . $no_of_disorder_region ."\t" .  $full_ordered_length; 
	} 	
	
	foreach $aa (@alphabet) {
		if ($aa eq undef || $occurence{'full'}->{$aa}==0) {
			print "\t.";
		}
		else {
			print   "\t" . $occurence{'full'}->{$aa} / $total_length ;
			
		}
		
		if ($aa eq undef || $occurence{'disorder'}->{$aa}==0) {
			print "\t.";
		}
		else {
			print   "\t" . $occurence{'disorder'}->{$aa} / $total_length;  # James second comment to divide by total aa
			
		}
		if ($full_ordered_length eq undef || $disorder_length eq undef) {
			print "\t.";
		
		}
		else {
		print "\t" . (( $occurence{'full'}->{$aa} - $occurence{'disorder'}->{$aa}) / $total_length );
		}
				
	}
		
	### COMPUTE properties in Total region
	$Tiny  =  	$seq  =~ tr/ACGSTacgst/ACGSTacgst/;
	$Small =  	$seq  =~ tr/ABCDGNPSTVabcdgnpstv/ABCDGNPSTVabcdgnpstv/;
	$Aliphaic =  	$seq  =~ tr/AILVailv/AILVailv/;
	$Aromatic =  	$seq  =~ tr/FHWYfhwy/FHWYfhwy/;
	$Nonpolar=   	$seq  =~ tr/ACFGILMPVWYacfgilmpvwy/ACFGILMPVWYacfgilmpvwy/;
	$Polar = 	$seq  =~ tr/DEHKNQRSTZdehknqrstz/DEHKNQRSTZdehknqrst/;
	$Charged = 	$seq  =~ tr/BDEHKRZbdehkrhz/BDEHKRZbdehkrhz/;
	$positive_charge = $seq  =~ tr/KRHkrh/KRHkrh/;
	$negative_charge = $seq  =~ tr/DEde/DEde/;
	$Basic =  	$seq  =~ tr/HKRhkr/HKRhkr/;
	$Acidic =  	$seq  =~ tr/BEDZbedz/BEDZbedz/;
	$large =  	$seq  =~ tr/NETVILPQHMFKWYRnetvilpqhmfkwyr/NETVILPQHMFKWYRnetvilpqhmfkwyr/;   
	$hydrophobic_Eisenberg = $seq  =~ tr/PYCGAMWLVFIpycgamwlvfi/PYCGAMWLVFIpycgamwlvfi/; 
	$Hydrophobic_KD =  $seq  =~ tr/AMCFLVIamcflvi/AMCFLVIamcflvi/; 	
	$Hydrophobic_FP = $seq  =~ tr/HTAPYVCLFIMWhtapyvclfmw/HTAPYVCLFIMWhtapyvclfmw/;	
	$helix = 	$seq =~ tr/VIYFWLviyfwl/VIYFWLviyfwl/;  
	$turn  = 	$seq =~ tr/NPGSnpgs/NPGSnpgs/; 
	$sheet = 	$seq =~ tr/EMALemal/EMALemal/; 
	
	### COMPUTE properties in disordered region
	
	$Tiny_dis  =  	$regions_disorder  =~ tr/ACGSTacgst/ACGSTacgst/;
	$Small_dis =  	$regions_disorder  =~ tr/ABCDGNPSTVabcdgnpstv/ABCDGNPSTVabcdgnpstv/;
	$Aliphaic_dis = $regions_disorder  =~ tr/AILVailv/AILVailv/;
	$Aromatic_dis = $regions_disorder  =~ tr/FHWYfhwy/FHWYfhwy/;
	$Nonpola_dis=  $regions_disorder  =~ tr/ACFGILMPVWYacfgilmpvwy/ACFGILMPVWYacfgilmpvwy/;
	$Polar_dis = 	$regions_disorder  =~ tr/DEHKNQRSTZdehknqrstz/DEHKNQRSTZdehknqrst/;
	$Charged_dis = 	$regions_disorder  =~ tr/BDEHKRZbdehkrhz/BDEHKRZbdehkrhz/;
	$positive_charge_dis = $regions_disorder  =~ tr/KRHkrh/KRHkrh/;
	$negative_charge_dis = $regions_disorder  =~ tr/DEde/DEde/;
	$Basic_dis =  	$regions_disorder  =~ tr/HKRhkr/HKRhkr/;
	$Acidic_dis =  	$regions_disorder  =~ tr/BEDZbedz/BEDZbedz/;
	$large_dis =  	$regions_disorder  =~ tr/NETVILPQHMFKWYRnetvilpqhmfkwyr/NETVILPQHMFKWYRnetvilpqhmfkwyr/;   
	$hydrophobic_Eisenberg_dis = $regions_disorder  =~ tr/PYCGAMWLVFIpycgamwlvfi/PYCGAMWLVFIpycgamwlvfi/; 
	$Hydrophobic_KD_dis =  $regions_disorder  =~ tr/AMCFLVIamcflvi/AMCFLVIamcflvi/; 	
	$Hydrophobic_FP_dis = $regions_disorder  =~ tr/HTAPYVCLFIMWhtapyvclfmw/HTAPYVCLFIMWhtapyvclfmw/;	
	$helix_dis = $regions_disorder =~ tr/VIYFWLviyfwl/VIYFWLviyfwl/;  
	$turn_dis  = $regions_disorder =~ tr/NPGSnpgs/NPGSnpgs/; 
	$sheet_dis = $regions_disorder =~ tr/EMALemal/EMALemal/; 
	
	### COMPUTE properties in Ordered region
	$Tiny_ord =	$Tiny - $Tiny_dis; 
	$Small_ord =	$Small - $Small_dis; 
	$Aliphatic_ord = $Aliphaic - $Aliphaic_dis;
	$Aromatic_ord = $Aromatic - $Aromatic_dis;
	$Nonpolar_ord =	$Nonpolar - $Nonpola_disr;
	$Polar_ord =	$Polar - $Polar_dis;
	$Charged_ord = 	$Charged  -	$Charged_dis;
	$positive_charge_ord =	$positive_charge - $positive_charge_dis;
	$negative_charge_ord =	$negative_charge - $negative_charge_dis;
	$Basic_ord =	$Basic - $Basic_dis;
	$Acidic_ord = $Acidic -	$Acidic_dis;
	$large_ord =	$large - $large_dis;
	$hydrophobic_Eisenberg_ord = $hydrophobic_Eisenberg - $hydrophobic_Eisenberg_dis;
	$Hydrophobic_KD_ord =	$Hydrophobic_KD - $Hydrophobic_KD_dis;
	$Hydrophobic_FP_ord =	$Hydrophobic_FP - $Hydrophobic_FP_dis; 
	$helix_ord = $helix - $helix_dis;  
	$turn_ord  = $turn - $turn_dis; 
	$sheet_ord = $sheet - $sheet_dis ;
	if ($total_length !=0 ) {
		print   "\t" .	$Tiny /$total_length;
		print   "\t" .  $Small / $total_length;
		print   "\t" . 	$Aliphaic / $total_length;
		print   "\t" . 	$Aromatic / $total_length ;
		print   "\t" . 	$Nonpolar/$total_length;
		print   "\t" . 	$Polar/$total_length;
		print   "\t" . 	$Charged/$total_length;
		print   "\t" . 	$positive_charge/$total_length;
		print   "\t" . 	$negative_charge/$total_length;
		print   "\t" . 	$Basic/$total_length;
		print   "\t" .	$Acidic/$total_length;
		print   "\t" . 	$large/$total_length;
		print   "\t" . 	$hydrophobic_Eisenberg/$total_length;
		print   "\t" . 	$Hydrophobic_KD/$total_length;
		print   "\t" . 	$Hydrophobic_FP/$total_length;
		print   "\t" . 	$helix/$total_length;
		print   "\t" . 	$turn/$total_length;
		print   "\t" .	$sheet/$total_length;
	}
	if ($Tiny !=0) {
		print   "\t" .	$Tiny_ord/$Tiny;
		print   "\t" .	$Tiny_dis/$Tiny;
	}	
	if ($Small !=0) {
		print   "\t" .	$Small_ord/$Small;
		print   "\t" .	$Small_dis/$Small;
	}	
	if ($Aliphaic !=0){
		print   "\t" .	$Aliphaic_ord/$Aliphaic;
		print   "\t" .	$Aliphaic_dis/$Aliphaic;
	}
	if ($Aromatic !=0 ) {
		print   "\t" .	$Aromatic_ord/$Aromatic;
		print   "\t" .	$Aromatic_dis/$Aromatic;
	}
	
		
	 if ($Nonpolar !=0) {	
		print   "\t" .	$Nonpolar_ord/$Nonpolar;
		print   "\t" .	$Nonpolar_dis/$Nonpolar;
	}
	if ($Polar!= 0) {	
		print   "\t" .	$Polar_ord/$Polar;
		print   "\t" .	$Polar_dis/$Polar;
	}
	if ($Charged !=0) {
		print   "\t" .	$Charged_ord/$Charged;
		print   "\t" .	$Charged_dis/$Charged;
	}	
	if ($positive_charge !=0) {	
		print   "\t" .	$positive_charge_ord/$positive_charge;
		print   "\t" .	$positive_charge_dis/$positive_charge;
	}
		
	if ($negative_charge !=0) {
	
		print   "\t" .	$negative_charge_ord/$negative_charge;
		print   "\t" .	$negative_charge_dis/$negative_charge;
	}	
	if($Basic !=0) {
		print   "\t" .	$Basic_ord/$Basic;
		print   "\t" .	$Basic_dis/$Basic;
	}	
	if($Acidic !=0){
		
		print   "\t" .	$Acidic_ord/$Acidic;
		print   "\t" .	$Acidic_dis/$Acidic;
	}			
	if ($large){
		
		print   "\t" .	$large_ord/$large;
		print   "\t" .	$large_dis/$large;
	}
		
	if($hydrophobic_Eisenberg !=0){
			
		print   "\t" .	$hydrophobic_Eisenberg_ord/$hydrophobic_Eisenberg;
		print   "\t" .	$hydrophobic_Eisenberg_dis/$hydrophobic_Eisenberg;
	}
	if($Hydrophobic_KD !=0) {
		
		print   "\t" .	$Hydrophobic_KD_ord/$Hydrophobic_KD;
		print   "\t" .	$Hydrophobic_KD_dis/$Hydrophobic_KD;
	}	
	if($Hydrophobic_FP !=0){
		
		print   "\t" .	$Hydrophobic_FP_ord/$Hydrophobic_FP;
		print   "\t" .	$Hydrophobic_FP_dis/$Hydrophobic_FP;
	}	
	if($helix !=0){
	
		print   "\t" .	$helix_ord/$helix;
		print   "\t" .	$helix_dis/$helix;
	}
	if ($turn !=0) {
		print   "\t" .	$turn_ord/$turn;
		print   "\t" .	$turn_dis/$turn;
	}	
	if($sheet) {
		print   "\t" .	$sheet_ord/$sheet;
		print   "\t" .	$sheet_dis/$sheet;
	}	 

print "\n";

}
	
print "\n";	
		
exit;				
			
