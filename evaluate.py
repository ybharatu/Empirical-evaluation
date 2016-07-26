import math
from dendropy import treesim
import dendropy
from dendropy import TreeList
import dendropy
import sys
import subprocess
import os

def ncr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def generateCoalescentTrees(choice, num, fout, length):
    	if choice == 1:
		sp_tree_str = """((((((((A:%f,B:%f):%f,C:%f):%f,D:%f):%f,E:%f):%f,F:%f):%f,G:%f):%f,H:%f):%f);""" % (float(length), float(length), float(length),float(length),2*float(length),float(length),3*float(length),float(length),4*float(length),float(length),5*float(length),float(length),6*float(length),float(length),7*float(length))

        	#sp_tree_str = """\
       		# [&R]  ((((((((A,B)%f,C)%f,D)%f,E)%f,F)%f,G)%f,H)%f);
        	#""" % (float(length),float(length),float(length),float(length),float(length),float(length),float(length))
    	elif choice == 2:
        	#sp_tree_str = """\
        	#[&R] (((A,B)%f,(C,D)%f)%f,((E,F)%f,(G,H)%f)%f);  
        	#""" % (float(length),float(length),float(length),float(length),float(length),float(length))
		sp_tree_str = """(((A:%f,B:%f):%f,(C:%f,D:%f):%f):%f,((E:%f,F:%f):%f,(G:%f,H:%f):%f):%f);""" % (float(length), float(length), float(length), float(length), float(length), 2*float(length),4*float(length),float(length), float(length),2*float(length),float(length), float(length), float(length),4*float(length)) 
    #print(sp_tree_str)
 	sp_tree = dendropy.Tree.get_from_string(sp_tree_str, "newick")
    	gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        	containing_taxon_namespace=sp_tree.taxon_namespace,
        	num_contained=1)
    	gene_tree_list = TreeList()

    	for i in range(num):
        	gene_tree = dendropy.simulate.treesim.contained_coalescent_tree(containing_tree=sp_tree,
        	gene_to_containing_taxon_map=gene_to_species_map)
        	dendropy.simulate.treesim.contained_coalescent_tree(containing_tree=sp_tree,
                                      gene_to_containing_taxon_map=gene_to_species_map)
        	for t in gene_tree.leaf_nodes():
            		t.taxon.label = t.taxon.label.split( )[0]
        	gene_tree_list.append(gene_tree)

   	gene_tree_list.write_to_path(fout, 'newick')

# Local Declarations
n = 8  # number of species
e = .05  # Error
f = float(sys.argv[1])  # Shortest branch length
choice = int(sys.argv[2])
estimate = 9 / (2 * (f ** 2)) * math.log(ncr(n, 4)/e*4)/ 20
m = int(estimate)

print(estimate)
print(f)
#generateCoalescentTrees(1,reps*int(m),"gene_true.tre",f)

var = 1
check = 0
sum_E = 0
reps = 50
generateCoalescentTrees(choice,reps*int(m),"gene_true.tre",f)
while(var == 1):
	print("======= m is now set to %f" %m)
	sum_E = 0
	
	str_com = "head -n " + str(int(m)*reps) + " gene_true.tre > subsampled_gene.tre"
	subprocess.call(str(str_com),shell=True)	
	str_com = "head -n " + str(int(m)*reps) + " gene_true.tre > sub_gene.tre"
	subprocess.call(str(str_com),shell=True)
	for x in range(reps):
		#generateCoalescentTrees(1, int(m),"gene.tre",f)
		str_com = "head -n " + str(int(m)) + " sub_gene.tre > gene.tre"
		subprocess.call(str_com,shell=True)
		str_com = "sed '1," + str(int(m)) + "d' subsampled_gene.tre > sub_gene.tre "
		subprocess.call(str_com,shell=True)
		subprocess.call("java -jar /home/ybharatu/ASTRAL/Astral/astral.4.10.5.jar -i gene.tre -o sp.tre",shell=True)
		if (choice == 1):
			str_com = "/home/ybharatu/global/src/shell/compare.tree.sh -s sp.tre -g sp_cat_" + str(float(f)) + ".tre > out.txt"
			subprocess.call(str_com,shell=True)
		elif (choice == 2):
			str_com = "/home/ybharatu/global/src/shell/compare.tree.sh -s sp.tre -g sp_bal_" + str(float(f)) + ".tre > out.txt"
			subprocess.call(str_com,shell=True)
		text_file = open("out.txt","r")
		lines = text_file.read().split('\n')
		nums = lines[0]
		num = nums.split(' ')
		compare = num[2]
		text_file.close()
		
		if (compare != "0"):
			sum_E = sum_E + 1
		#print("SUM of Error = %f" %int(sum_E))	
        print ("Number of error cases is %f" %sum_E)
	print(sum_E / float(reps))
	if (sum_E / float(reps)) > e:
		m = m * (1.5)
		print("3/2")
		print(m)
		#generateCoalescentTrees(1,int(m),"gene.tre",f)
		check = 1
	elif(sum_E / float(reps)) < e:
		if (check == 1):
			var = 0
		else:
			m = m / 2
			print("1/2")
			print(m)
			#generateCoalescentTrees(1,int(m),"gene.tre",f)
		
print(m)
