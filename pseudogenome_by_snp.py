# Generate psudo-genome based on SNP index
# Usage: pseudogenome_by_snp.py genome.fasta SNP.index out.file.basename
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

genome_filename = sys.argv[1]
snp_filename = sys.argv[2]
output_filenamebase = sys.argv[3]

# open and read snp index file
snp_file = open(snp_filename)
snps = snp_file.readlines()
snps_header = [k for k in snps if "#" in k][0].strip("#|\n").split("\t")


# open output files using a file handle
outputA_handle = open(output_filenamebase + "." + snps_header[2] + ".fasta", "w")
outputD_handle = open(output_filenamebase + "." + snps_header[3] + ".fasta", "w")

# Read and parse fasta file
for record in SeqIO.parse(genome_filename, "fasta"):
    # get corresponding snps
    snps_for_record = [k for k in snps if record.id in k]
    record_Aseq = record.seq.tomutable()
    record_Dseq = record.seq.tomutable()
    for each in snps_for_record:
        [pos,A,D] = each.strip("\n").split("\t")[1:4]
        pos = int(pos)
# assign SNP to corresponding sequence location
        record_Aseq[pos-1:pos] = A
        record_Dseq[pos-1:pos] = D
# save modified sequence with edited ID
    record_A=SeqRecord(record_Aseq, id=record.id + "_" + snps_header[2], description = '')
    record_D=SeqRecord(record_Dseq, id=record.id + "_" + snps_header[3], description = '')
# output changed seq
    SeqIO.write(record_A, outputA_handle, "fasta")
    SeqIO.write(record_D, outputD_handle, "fasta")


snp_file.close()
outputA_handle.close()
outputD_handle.close()
