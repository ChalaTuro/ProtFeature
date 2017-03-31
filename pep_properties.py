#!/sw/bin/python3

import sys 
from Bio import SeqIO
from Bio.SeqUtils import ProtParamData  
from Bio.SeqUtils import IsoelectricPoint  
from Bio.Seq import Seq 
from Bio.Alphabet import IUPAC 
from Bio.Data import IUPACData 
from Bio.SeqUtils import molecular_weight 
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParam
from Bio.SeqRecord import SeqRecord
from __future__ import print_function 
from contextlib import suppress

file to write output

myfile = open("/home/chala/scripts/properties/303_09.txt","w") 

print ("Seq_ids", "molecular_weight","aromacity","instability","Isoelectric","helix", "turn","sheet", "gravy",file=myfile )

# import fasta file
records = list(SeqIO.parse("/media/chala/SEGEAT2/PROGRAM10/ProFET/Secretome/fasta/303.fungal_SSP.09.fa", "fasta"))

## replace unusual protein characters that the module did not accept and print analysis output.

for record in records:
    A = str(record.seq)
    B = A.replace("X", "")
    C = B.replace("B", "")
    D = C.replace("Z", "") 
    E = D.replace("J", "") 
    F = E.replace("U" , "")
    G = F.replace("O" , "")
    #Y = ProteinAnalysis(str(record.seq))
    X = ProteinAnalysis(G)
    print(record.id,X.molecular_weight(), X.aromaticity(), X.instability_index(), X.isoelectric_point(), X.secondary_structure_fraction(),X.gravy(),file=myfile) 
    
