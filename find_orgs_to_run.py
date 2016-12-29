import sys
import re, pprint
import os
from urllib2 import urlopen
from bs4 import BeautifulSoup
from sets import Set

def read_metadocs_files (folder, date):
	new_genomes_dict = Set()
	genomes_major_vers_dict = {} 
#	print "In metadocs"
	for metadoc_file in os.listdir (folder):
		if date in metadoc_file:
			organism, genome_build, major_genome_version = read_metadoc_file (folder, metadoc_file)
			new_genomes_dict.add (organism.lower() + "/" +  genome_build.lower() )
			genomes_major_vers_dict[organism.lower()] = major_genome_version
#			print organism + " major " +  major_genome_version
	return (new_genomes_dict, genomes_major_vers_dict)

def public_SIFT_genomes (sift_folder):

	public_db = Set()
	public_major_version =  {}

	urlpath = urlopen (sift_folder)
	string = urlpath.read().decode ('utf-8')
	#print string

	soup = BeautifulSoup (string)
	for a in soup.findAll ('a'):
		organism =  a['href']
		if organism.endswith ("/") and not "sift4g" in organism:
#			print "looking for " + organism
			subdir = urlopen (sift_folder + "/" + organism)
			string2 = subdir.read().decode ('utf-8')
			soup2 = BeautifulSoup (string2)	
			for a2 in soup2.findAll ('a'):
				build = a2['href']
		#	print organism + "\t" + a2['href']
				if build.endswith (".zip"): 
					build = build.replace (".zip", "")
#				print organism + "\t" + build	
					public_db.add (organism.lower() + build.lower())
					major_vers = get_major_genome_version (build)
#					print "public organism " + organism + " build " + build + " vers " + major_vers
					key = organism.lower().strip ("/")
					# update to the latest major version on the website
					if key in public_major_version:
						if int (major_vers) > int (public_major_version[key]):
							public_major_version[key] = major_vers
					else:
						public_major_version[key] = major_vers
	return (public_db, public_major_version)

def get_major_genome_version (genome_build ):
	major_genome_build = ""
	if genome_build != "":
		genome_build_sep = genome_build.split ("/.")[0]
#		print "genome build is " + genome_build
		major_genome_build = re.search ("\d+" , genome_build_sep).group()
#		org_name = re.search ("[a-zA-Z_]", genome_build_sep).group()
	return major_genome_build

def read_metadoc_file (folder, query):
	
	filename = folder + "/" + query
	organism = ""
	genome_build = ""
	major_genome_build = ""
	fp = open (filename, "r")
	for line in fp.readlines():
		if line.startswith ("ORG"):
			fields = line.rstrip().split("-")
			organism = fields[1]
		if line.startswith ("ORG_VERSION"):
			fields = line.rstrip().split ("-")
			genome_build = fields[-1]
			genome_build_sep = genome_build.split ("/.")[0]
			major_genome_build = get_major_genome_version (genome_build) 
			
#	print ("major " + str(major_genome_build) + " for " + str(genome_build))
	fp.close()

	return (organism, genome_build, major_genome_build)
	
##################################
# main part

folder = "metadocs/"
date = "2016-08-12"
sift_public_db = "http://sift.bii.a-star.edu.sg/sift4g/public/"

new_genomes_set, new_major_genome_dict = read_metadocs_files (folder, date)
#print "new genomes set"
#print new_genomes_set

current_genomes_set, current_major_genome_dict = public_SIFT_genomes (sift_public_db)

#print "Genomes Done"
#for genome_done in new_genomes_set.intersection (current_genomes_set):
#	print genome_done	

#print "********"
#print "Genomes To Do"
#for missing_organism in new_genomes_set.difference (current_genomes_set):
#	print missing_organism

#print "*******"
print "Metadocs To Do"
for missing_organism in new_genomes_set.difference (current_genomes_set):
        fields = missing_organism.split ("/")
	org_name = fields[-2]
	print org_name + "-" + date + ".txt"

print "***********************\n\n"

print "Major Versions to Update"
for organism, major_assembly in new_major_genome_dict.items():
	if organism in current_major_genome_dict : 
		current_major = int (current_major_genome_dict[organism])
#		print organism + " current major " + str (current_major) + " major assembly " + major_assembly
		if major_assembly != "" and current_major != "" and int (major_assembly) > current_major:
			print "Major version update for " + organism + " from " + current_major_genome_dict[organism] + " to "+ major_assembly
	else: 
		print "organism " + organism + " not in public and needs to be done "
