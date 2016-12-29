import sys
import re, pprint
import os
import gzip

def get_pos (line):
	fields = line.split ('\t')
	return int (fields[0])

# read in scores, return sorted array
def get_regions (filename, outfile):
	# escape gracefully if file doesn't exist
	if not os.path.isfile (filename):
		return

	outfp = open (outfile, "w")
		
	# initialize
	out_start = 1
	out_end = -1
	in_start = -1
	in_end = -1

#	while f:	
	with open (filename) as f:
#		f.readline() # read header , only for dbNSFP
		first_line = f.readline()
		# skip comment lines
		while first_line.startswith ("#"):
			first_line = f.readline()
		pos = get_pos (first_line)
		in_start = pos
		out_end = pos -1
	
		cur_pos = pos
		prev_pos = pos

		for line in f:
			fields = line.split ('\t')
			cur_pos = get_pos (line) 
			if cur_pos == prev_pos or cur_pos == (prev_pos + 1):
				prev_pos = cur_pos
			else:
				# entered a new region
				outfp.write ( str (out_start) + "\t" + str (in_start -1) + "\tOUT\n") 
				outfp.write ( str (in_start) + "\t" + str (prev_pos)  + "\tIN\n") 
				# update
				in_start = cur_pos
				out_start = prev_pos + 1	
				prev_pos = cur_pos
	outfp.write ( str (out_start) + "\t" + str (in_start -1) + "\tOUT\n")
	outfp.write ( str (in_start) + "\t" + str (prev_pos)  + "\tIN\n")

		
	outfp.close()

##################################
# main part

if len (sys.argv) != 2: 
        print 'check_SIFTDB.py <sorted file> <outfile>'

#print sys.argv
chrom_file = sys.argv[1]
out_file = sys.argv[2]

print "looking at " + chrom_file
get_regions (chrom_file, out_file)

