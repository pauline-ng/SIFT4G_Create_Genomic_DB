import os
import sys
import re
import getopt

global REF_BASE_IDX
REF_BASE_IDX = 1 
global NEW_BASE_IDX
NEW_BASE_IDX = 2 
global REF_AA_IDX 
REF_AA_IDX = 7 
global NEW_AA_IDX 
NEW_AA_IDX =8 
global SIFT_SCORE_IDX
SIFT_SCORE_IDX = 10
global SIFT_CONF_IDX
SIFT_CONF_IDX = 11
global DBSNP_IDX 
DBSNP_IDX = 13

class SIFTPredMetric (object):
	num_damaging = 0
	num_tolerated = 0
	num_not_predicted = 0
	total = 0

	def __init__ (self):
		self.num_damaging = 0
		self.num_tolerated = 0
		self.num_not_predicted = 0
		self.total = 0

	def __repr__ (self):
		predicted_on = self.num_damaging + self.num_tolerated
		if predicted_on > 0:
#			print str (self.num_damaging) + " " + str (self.num_tolerated)
			return "%.1f" % (float(self.num_damaging) / float (predicted_on) * 100.0)	
		else:
			return "-1"

	def __str__ (self):
		predicted_on = self.num_damaging + self.num_tolerated
#		print str (self.num_damaging) + " " + str (self.num_tolerated)
                if predicted_on > 0:
#			print str (self.num_damaging) + " " + str (self.num_tolerated) 
			return "%.1f" % ( float (self.num_damaging) / float (predicted_on) * 100.0)
                else:
                        return "-1"

	def percent_predicted_on (self):
		if self.total > 0:
			return ((float (self.num_damaging) + float (self.num_tolerated))  *100)/ float (self.total) 
		else:
			return -1 
	
	def counts (self):
		predicted_on = self.num_damaging + self.num_tolerated
		return str (self.num_damaging) + "/" +  str (predicted_on) 

	def inc_dam (self):
		self.num_damaging += 1
#		print "dam value now " + str (self.num_damaging)

	def inc_tol (self):
		self.num_tolerated += 1	
#		print "tol now " + str (self.num_tolerated)

	def inc_not_pred (self):
		self.num_not_predicted += 1

	def inc_total (self):
		self.total += 1

	def inc_another_metric (self, m):
	        self.num_damaging += m.num_damaging
        	self.num_tolerated += m.num_tolerated
	        self.num_not_predicted += m.num_not_predicted
        	self.total += m.total

def output_results (outfile, DNA_base_change_syn_nonsyn_counts, DNA_base_change_aa_pred, all_nonsyn_changes, dbSNP_predictions, novel_predictions, ref_predictions, aa_change_pred):
        outfp = open (outfile, "w")
        outfp.write (outfile + "\n")
	for key, value in DNA_base_change_syn_nonsyn_counts.iteritems():
                outfp.write (key + "\t" + str (value) + "\n")
        outfp.write ("\n")
	if DNA_base_change_syn_nonsyn_counts["SYN"] == 0 and DNA_base_change_syn_nonsyn_counts["NONSYN"] == 0:
		outfp.write ("no synonymous or nonsynonymous changes, empty file?")
		outfp.close()
		return	
        ns_syn_ratio = float (DNA_base_change_syn_nonsyn_counts["NONSYN"])/float (DNA_base_change_syn_nonsyn_counts["SYN"])
        outfp.write ("nonsyn/syn ratio: %.2f\n\n" % ns_syn_ratio)

        outfp.write ("All Counts: ")
	outfp.write (str (all_nonsyn_changes) + " " + all_nonsyn_changes.counts() + "\n")
#print str (all_nonsyn_changes)

	outfp.write ("dbSNP damaging: " )
	outfp.write (str (dbSNP_predictions) + "% " +  dbSNP_predictions.counts()  +  "\n")


	outfp.write ("Not in dbSNP damaging: ")
	outfp.write (str (novel_predictions) + "% " + novel_predictions.counts() +  "\n")
	outfp.write ("\n")

	outfp.write ("Reference predictions % damaging: ")
	outfp.write (str (ref_predictions) + " " + ref_predictions.counts() +  "\n")
	outfp.write ("% predicted on: " + str (ref_predictions.percent_predicted_on()) ) 
	outfp.write ("\n")

	output_matrix_of_SIFTPredMetrics (DNA_bases, DNA_base_change_aa_pred, outfp)
	#output_matrix_of_SIFTPredMetrics (amino_acids, aa_change_pred, outfp)
	outfp.write ("#####################################\n")
	outfp.close()



def get_SIFT_prediction ( sift_score, sift_confidence):
#	print "sift score " + str (sift_score) + "confidence" + sift_confidence
	if not sift_score or not sift_confidence or sift_score == "NA" or sift_confidence == "NA":
		return ""
	if float (sift_score) < 0.05 and float (sift_confidence) < 3.5:
        	return "DAMAGING"
	elif float (sift_score) >= 0.05:
                return "TOLERATED" 
	else:
		return "NONE"

def  tabulate_dbSNP_changes (fields, dbSNP_predictions, novel_predictions, ref_predictions):
	ref_aa = fields[REF_AA_IDX]
        new_aa = fields[NEW_AA_IDX]
	ref_base = fields[REF_BASE_IDX]
	new_base = fields[NEW_BASE_IDX]
	# not scoring stops or U's
        if ref_aa == "*" or new_aa == "*" or ref_aa == "U" or new_aa == "U" or ref_aa == "" or new_aa == "" or ref_aa == "NA" or new_aa == "NA":
		return	
	# scoring synonymous
	if ref_aa == new_aa:
		if ref_base == new_base: # don't want to count codon changes to same amino acid, as then there is ovecounting and a position represented more than once
			var_to_update = ref_predictions 
		else:
			return 
	elif fields[DBSNP_IDX].startswith ("rs"): # nonsynonymous
		var_to_update = dbSNP_predictions 
		print "updating with " + fields[DBSNP_IDX]
		print ' '.join (fields)
	elif fields[DBSNP_IDX].startswith ("novel"): #nonsynynmous
		var_to_update = novel_predictions
	elif fields[DBSNP_IDX] != "ref":
		var_to_update = dbSNP_predictions
		print "updating nonstandard " + fields[DBSNP_IDX]
	else:
		print "dbsnp field " + fields[DBSNP_IDX] + "\n"
		return	# this is so ref changes aren't counted

	var_to_update.inc_total()
	pred = get_SIFT_prediction (fields[SIFT_SCORE_IDX], fields[SIFT_CONF_IDX])
	if pred == "DAMAGING": 
		var_to_update.inc_dam () 
	elif pred == "TOLERATED":
		var_to_update.inc_tol() 
	elif pred == "NONE":
		var_to_update.inc_not_pred()
	
def tabulate_aa_changes (fields, aa_change_pred):
#	print "in tabulate_aa_changes"
#	print fields
	ref_aa = fields[REF_AA_IDX]
	new_aa = fields[NEW_AA_IDX]
#	print "amino acid change " + ref_aa + "->" + new_aa
	if ref_aa == "*" or new_aa == "*" or ref_aa == "U" or new_aa == "U" or ref_aa == "" or new_aa == "" or ref_aa == "NA" or new_aa == "NA":
		return
	aa_change_pred[ref_aa][new_aa].inc_total()
	pred = get_SIFT_prediction (fields[SIFT_SCORE_IDX], fields[SIFT_CONF_IDX])
        if pred == "DAMAGING":
		aa_change_pred[ref_aa][new_aa].inc_dam ()
	elif pred == "TOLERATED":
		aa_change_pred[ref_aa][new_aa].inc_tol() 
	elif pred == "NONE":
		aa_change_pred[ref_aa][new_aa].inc_not_pred()

def output_matrix_of_SIFTPredMetrics (indices, matrix, out_fp):
	# print header
	out_fp. write ("\t") 
	for index in indices:
		out_fp.write ( index + "\t") 
		
	out_fp.write ( "\n") 

	# now print values
	for index1 in indices:
		out_fp.write (index1 + "\t") 
		for index2 in indices:
			pred = matrix[index1][index2]			 
			out_fp.write (str (pred) + "\t")
		out_fp.write ( "\n") 
	out_fp.write ("\n")

def process_chr_file (siftdb_infile, DNA_base_change_syn_nonsyn_counts, DNA_base_change_aa_pred, all_nonsyn_changes, dbSNP_predictions, novel_predictions, ref_predictions, aa_change_pred):

        siftdb_fp = open (siftdb_infile, "r")
        already_processed = set ()  # store positions that have been processed
	for line in siftdb_fp.readlines():
       		fields = line.split ('\t')
	        variant_key = ':'.join (fields[0:3])
#		print variant_key + "\n"
		if variant_key not in already_processed and len(fields) == 14 and fields[NEW_AA_IDX] != "" and fields[REF_AA_IDX] != "" and not line.startswith("#"):
#               		print "processing " + variant_key
#			print line
			tabulate_DNA_changes (fields, DNA_base_change_syn_nonsyn_counts, DNA_base_change_aa_pred, all_nonsyn_changes)
               		tabulate_dbSNP_changes (fields, dbSNP_predictions, novel_predictions, ref_predictions)
	                tabulate_aa_changes (fields, aa_change_pred)
       	        	already_processed.add (str(variant_key))
#			print already_processed
        siftdb_fp.close()

def read_metadoc (metadocfile):
	metahash = {}
	fo = open (metadocfile, "r")
	lines = fo.readlines()
	for line in lines:
		line = line.strip()
#		print "in here " + line
		if "=" in line and not line.startswith ("#"):
#			print "in here2 " + line
			var, val = line.split ("=")
			var = var.strip()
			val = val.strip()
			metahash[var] = val
	fo.close()
	return metahash
	
def tabulate_DNA_changes (fields, DNA_base_change_syn_nonsyn_counts, DNA_base_change_aa_pred, all_nonsyn_changes): 
	ref_allele = fields[REF_BASE_IDX]
	new_allele = fields[NEW_BASE_IDX]
	ref_aa = fields[REF_AA_IDX]
	new_aa = fields[NEW_AA_IDX]	
	if ref_aa != "" and new_aa != "" and ref_allele != new_allele: # don't care about reference bases
		if ref_aa == new_aa:
			DNA_base_change_syn_nonsyn_counts["SYN"] += 1
		elif ref_aa != new_aa and new_aa == "*":
			DNA_base_change_syn_nonsyn_counts["STOP-GAINED"] += 1
		elif ref_aa != new_aa: 
			DNA_base_change_syn_nonsyn_counts["NONSYN"] += 1
			all_nonsyn_changes.inc_total () 
			DNA_base_change_aa_pred[ref_allele][new_allele].inc_total() 
			pred = get_SIFT_prediction (fields[SIFT_SCORE_IDX], fields[SIFT_CONF_IDX])
		        if pred == "DAMAGING":
				all_nonsyn_changes.inc_dam () 
				DNA_base_change_aa_pred[ref_allele][new_allele].inc_dam()
			elif pred == "TOLERATED" :
				all_nonsyn_changes.inc_tol() 
				DNA_base_change_aa_pred[ref_allele][new_allele].inc_tol() 
			elif pred == "NONE":
				all_nonsyn_changes.inc_not_pred()
				DNA_base_change_aa_pred[ref_allele][new_allele].inc_not_pred()
		else:
			print "ERROR this case is not handled " + ' '.join (fields)

def initialize_metrics (DNA_base_change_syn_nonsyn_counts, DNA_base_change_aa_pred, all_nonsyn_changes, dbSNP_predictions, novel_predictions, ref_predictions, aa_change_pred):
	DNA_base_change_syn_nonsyn_counts["SYN"] = 0
	DNA_base_change_syn_nonsyn_counts["STOP-GAINED"] = 0
	DNA_base_change_syn_nonsyn_counts["NONSYN"] = 0

# main 
#initailize everything

DNA_bases = ["A", "C", "G", "T"]
amino_acids = "A C D E F G H I K L M N P Q R S T V W Y".split (" ")
#print "amino aicds " + str (len (amino_acids))
DNA_base_change_aa_pred = {} # 2x2 dictionary of A>T damaging, total, etc
for base1 in DNA_bases:
	DNA_base_change_aa_pred[base1] ={}
	for base2 in DNA_bases:
		DNA_base_change_aa_pred[base1][base2] = SIFTPredMetric ()

aa_change_pred = {} # 20x20 matrix for deleterious/tolerated
for aa1 in amino_acids:
	aa_change_pred[aa1] = {}
	for aa2 in amino_acids:
		aa_change_pred[aa1][aa2] = SIFTPredMetric()
	 
DNA_base_change_syn_nonsyn_counts = {}
DNA_base_change_syn_nonsyn_counts["SYN"] = 0 
DNA_base_change_syn_nonsyn_counts["STOP-GAINED"] = 0 
DNA_base_change_syn_nonsyn_counts["NONSYN"] = 0 

all_nonsyn_changes = SIFTPredMetric()
dbSNP_predictions = SIFTPredMetric()
novel_predictions = SIFTPredMetric()
ref_predictions = SIFTPredMetric()

#siftdb_infile = "tmp8"
#outfile = siftdb_infile + ".SIFTstats"

if len (sys.argv) != 3:
	print 'check_SIFTDB.py <infile> <outfile>'

infile = sys.argv[1]
outfile = sys.argv[2]

process_chr_file (infile, DNA_base_change_syn_nonsyn_counts, DNA_base_change_aa_pred, all_nonsyn_changes, dbSNP_predictions, novel_predictions, ref_predictions, aa_change_pred) 
output_results (outfile, DNA_base_change_syn_nonsyn_counts, DNA_base_change_aa_pred, all_nonsyn_changes, dbSNP_predictions, novel_predictions, ref_predictions, aa_change_pred)	
