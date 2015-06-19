#!/usr/bin/env python

#This script downloads all NCBI Fungal genomes with protein annotation
#written by Jon Palmer palmer.jona at gmail dot com

import sys
import os
import argparse
import urllib2, shutil
import os.path
from Bio import Entrez
from xml.dom import minidom
import gzip
from Bio import SeqIO
import glob

parser=argparse.ArgumentParser(
    description='''Script searches NCBI Entrez and downloads annotated genomes ''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""")
parser.add_argument('search_term', help='Organism search term')
parser.add_argument('email', help='must tell ncbi who you are')
parser.add_argument('-n', '--num_records', default='1000', help='max number of records to download')
args=parser.parse_args()

max_num = args.num_records

#construct search term based on input
organism = args.search_term
term = '(((%s[Organism]) AND "genome sequencing"[Filter]) AND "bioproject protein"[Filter]) AND "bioproject assembly"[Filter]' % (organism)

Entrez.email = args.email   # Always tell NCBI who you are
handle = Entrez.esearch(db="bioproject", term=term, retmax=max_num)
search_result = Entrez.read(handle)
handle.close()
count = int(search_result["Count"])
print "------------------------------------------------"
print("Found %i records from BioProject." % (count))
question = raw_input('Do you want to continue?: [y/n] ')
if question == "n":
    print "Script Exited by User"
    os._exit(1)
else:
    print "------------------------------------------------"
    print "Using ELink to get assembly ID numbers."
    print "------------------------------------------------"
    gi_list = ",".join(search_result["IdList"])
    #get the link data from assembly and print list to temp file
    link_result = Entrez.read(Entrez.elink(dbfrom="bioproject", db="assembly", id=gi_list))
    fo = open('gi_list.tmp', 'a')
    for link in link_result[0]["LinkSetDb"][0]["Link"]: 
        fo.write("%s\n" % (link["Id"]))
    fo.close()

    #load in temp file to variable and get esummary for each record?
    with open('gi_list.tmp') as fi:
        link_list = fi.read().splitlines()
    fi.close()
    gi_link = ",".join(link_list)
    summary_handle = Entrez.esummary(db="assembly", id=gi_link, report="xml")
    out_handle = open('results.xml', "w")
    out_handle.write(summary_handle.read())
    summary_handle.close()
    out_handle.close()
    print "Feteched Esummary for all genomes, saved to temp file results.xml"
    print "------------------------------------------------"
    #now need to parse those records to get Assembly accession numbers and download genome
    data = minidom.parse("results.xml")
    node = data.documentElement
    ids = data.getElementsByTagName("DocumentSummary")
    for id in ids:
        accession = id.getElementsByTagName("AssemblyAccession")[0].childNodes[0].data
        name = id.getElementsByTagName("AssemblyName")[0].childNodes[0].data
        species = id.getElementsByTagName("Organism")[0].childNodes[0].data
        folder = accession + "_" + name.replace(" ", "_")
        address = folder + "_genomic.gbff.gz"
        out_name = folder + ".gbk"
        if os.path.isfile(out_name):
            print "Genome for " + species + " already downloaded and converted.  Skipping"
            break
        else:
            if os.path.isfile(address):
                print "Genome for " + species + " already downloaded, now converting."
                gbf_file = gzip.open(address, 'rb')
                gb_out = open(out_name, "w")
                gb_file = SeqIO.parse(gbf_file, "genbank")
                out_file = SeqIO.write(gb_file, gb_out, "genbank")
                gbf_file.close()
                gb_out.close()
            else:
                print("Found: " + species + ".  Downloading file:  " + address)
                try:
                    ftpfile = urllib2.urlopen("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/" + folder + "/" + address)
                    localfile = open(address, "wb")
                    shutil.copyfileobj(ftpfile, localfile)
                    localfile.close()
                    gbf_file = gzip.open(address, 'rb')
                    gb_out = open(out_name, "w")
                    gb_file = SeqIO.parse(gbf_file, "genbank")
                    out_file = SeqIO.write(gb_file, gb_out, "genbank")
                    gbf_file.close()
                    gb_out.close()
                except:
                    print "Error: file not found for " + species + " " + address
    #clean up your mess
    os.remove('results.xml')
    os.remove('gi_list.tmp')
    filelist = glob.glob("*.gbff.gz")
    for f in filelist:
        os.remove(f)
    
    total = len([name for name in os.listdir('.') if os.path.isfile(name)])
    print "------------------------------------------------"
    print ("Searched for %s:" % (term))
    print "------------------------------------------------"
    print ("Resulted in downloading and unpacking %i files" % (total))