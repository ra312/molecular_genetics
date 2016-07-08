
import sys

if len(sys.argv)!= 2:
 sys.exit( 'Enter file name:')
xml_file = sys.argv[1]


import untangle

doc = untangle.parse(xml_file)
f = open('selected.txt', 'w')
for Gene in doc.GenesList.Gene:
	print >> f, Gene['chr'],"\t", Gene['start'],"\t", Gene['end'],"\t", Gene['name'] 
f.close()
