import twobitreader
genome = twobitreader.TwoBitFile('hg19.2bit')
bedfile = open('hg19.bed', 'w')
o_f = open('hg_FASTA.fa', 'w')
twobitreader.twobit_reader(genome, bedfile,o_f.write)
bedfile.close()
o_f.close()


