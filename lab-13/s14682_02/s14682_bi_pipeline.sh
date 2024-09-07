

#Simple pipeline to combine the 4 scripts in the given order
python.exe cds_seq_retrieve.py
python.exe transcribe.py
python.exe translate.py
python.exe aa_seq_analyze.py


#Create two folders
mkdir intermediate_files
mkdir output

#Move the output files to the relevent folder
mv cds_seq.fasta intermediate_files/
mv mRNA_seq.fasta intermediate_files/
mv aa_seq.fasta intermediate_files/
mv aa_stats.txt output/
