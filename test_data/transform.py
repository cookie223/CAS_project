from Bio import SeqIO

sequence = SeqIO.parse('globins4.sto', 'stockholm')
SeqIO.write(sequence, 'globins4.fasta', 'fasta')