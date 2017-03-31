#! usr/bin/perl

##To split ProFET generated features into separate classes

##usage: script.pl filein

open(IN,'<', @ARGV[0]);
open(OUT1, '>', 'Cytoplasmic.csv' );
open(OUT2, '>', 'Omycetes .csv' );
open(OUT3, '>', 'Defensins.csv' );
open(OUT4, '>', 'Appoplastic.csv' );
open(OUT5, '>', 'All_known_effector.csv' );
open(OUT6, '>', 'MAX_lik.csv' );
open(OUT7, '>', 'Pore_forming.csv' );
open(OUT8, '>', 'All_fungal.csv' );
open(OUT9, '>', 'Hydrophobin .csv' );
open(OUT10, '>', 'AvrLM76.csv' );
open(OUT11, '>', 'AMP.csv' );

#open(OUT12, '>', 'PHIbase.csv' );     run separately with VFDB
#open(OUT13, '>', 'VFDB_fixed.csv' );

open(OUT14, '>', 'TOXA.csv' );

@line1 = undef;	

@cytoplasmic = undef;
@Omycetes = undef;
@defensins = undef;
@appoplastic =undef;
@All_known_effector= undef;
@MAX_like = undef;
@pore_forming =undef;
@All_fungal =undef;
@hydrophobin =undef;
@AvrLM76 =undef;
@AMP = undef;
@PHIbase = undef;
@VFDB = undef;
@ToxA = undef;
 
#print "IDs\tEffectorP\tProb\n";  # data header

while ($line = <IN>) {          
	#chomp $line;

	#next if $line =~ /^\s*$/;  # ignore empty lines
	
	if ($line =~ /proteinname/) {
		push(@line1, $line);
		}
		
	if ($line =~ /fungal_cytoplasmic_secreted_effector.newid/g) {
		 push(@cytoplasmic, $line);
		 }
	if ($line =~ /Omycetes_all_effector_secreted.newid/g) {
		 push(@Omycetes, $line);
		 }
	if ($line =~ /Defensin_IPR0006080_matches.newid/g) {
		 push(@defensins, $line);
		  
	  }

	if ($line =~ /fungal_appoplastic_effector_secreted.newid/g) {
		 push(@appoplastic , $line);
		 }
	if ($line =~ /All_known_secreted_effector/g) {
		 push(@All_known_effector , $line);	 
		 
	  }
	if ($line =~ /MAX_like.newid/g) {
		 push(@MAX_like , $line);	 
		 
	  }	
	  if ($line =~ /Pore_forming.newid/g) {
		 push(@pore_forming , $line);	 
		 
	  }
	  if ($line =~ /All_fungal_secreted_effectors.newid/g) {
		 push(@All_fungal , $line);	 
		 
	  }
	  if ($line =~ /hydrophobin.newid/g) {
		 push(@hydrophobin , $line);	 
		 
	  }
	  if ($line =~ /AvrLM76_like_UPPER.newid/g) {
		 push(@AvrLM76, $line);	 
		 
	  }
	  if ($line =~ /Antimicrobial_ADP3.2750.newid/g) {
		 push(@AMP, $line);	 
		 
	  }
	  
	  if ($line =~ /PHIbase2.newid/g) {
		 push(@PHIbase, $line);	 
		  }
	  
	   if ($line =~ /VFDB_fixed.newid/g) {
		 push(@VFDB, $line);	 
		  }
	    if ($line =~ /ToxA_like.fasta.newid/g) {
		 push(@ToxA, $line);	 
		  }
}

print OUT1 @line1;
print OUT1 @cytoplasmic;

print OUT2 @line1;
print OUT2 @Omycetes;
 
print OUT3 @line1;
print OUT3 @defensins; 

print OUT4 @line1;
print OUT4 @appoplastic;

print OUT5 @line1;
print OUT5 @All_known_effector;

print OUT6 @line1;
print OUT6 @MAX_like;

print OUT7 @line1;
print OUT7 @pore_forming;

print OUT8 @line1;
print OUT8 @All_fungal;

print OUT9 @line1;
print OUT9 @hydrophobin; 

print OUT10 @line1;
print OUT10 @AvrLM76; 

print OUT11 @line1;
print OUT11 @AMP; 

#print OUT12 @line1;
#print OUT12 @PHIbase;

#print OUT13 @line1;
#print OUT13 @VFDB ;

print OUT14 @line1;
print OUT14 @ToxA;
	
close IN;


