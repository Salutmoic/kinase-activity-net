from os import listdir
from os.path import isfile, join
import sys
#parses CLFinder output directories and collects the results for further analysis
def add_commun(clusts,f,no):
	with open(f,"r") as v:
		for line in v:
			print(line)
			if "#" not in line or line != "\n":
				line = line.split(" ")
				clusts[line[0].replace(":",""),no] = line[1:]
								
def add_cliqs(clusts,f):
	with open(f,"r") as v:
		for line in v:
			if "#" not in line or line != "\n":
				
				line = line.split(" ")
				clusts[line[0].replace(":","")] = line[1:]

def write_dict(clusts,f):
	with open(f,"w") as v:
		for key in clusts:
			if type(key) is tuple:
				s = key[0] + "\t" + key[1]
			elif type(key) is str:
				s = key
			if "#" in s:
				print(s)
				continue
			for kinase in clusts[key]:
				v.write(s + "\t" + kinase.replace("\n","") + "\n")
		
if __name__ == "__main__":
	path = "CFinder-2.0.6--1448/"
	netw= sys.argv[1]
	dirs = listdir(path+netw)
	
	cliq_clusts = {}
	add_cliqs(cliq_clusts,path+netw+"/"+"cliques")
	write_dict(cliq_clusts,"out/clique-file.tsv")
	
	commun_clusts = {}

	for d in dirs:
		if 'k=' in d:
			no = d[2:]
			f = path+netw+"/" +d+"/"+ "communities" 
			 
			add_commun(commun_clusts,f,no)
	write_dict(commun_clusts,"out/commun-file.tsv")
