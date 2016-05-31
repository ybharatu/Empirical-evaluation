import math
from dendropy.simulate import treesim
import dendropy
from dendropy import TreeList
import sys
import subprocess


def ncr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


def generateCoalescentTrees(choice, num, fout, length):
    if choice == 1:
        sp_tree_str = """\
        [&R]  ((((((((A,B)%d,C)%d,D)%d,E)%d,F)%d,G)%d,H)%d);
        """ % (float(length),float(length),float(length),float(length),float(length),float(length),float(length))
    elif choice == 2:
        sp_tree_str = """\
        [&R]  ((((((((A,B)%d,C)%d,D)%d,E)%d,F)%d,G)%d,H)%d);
        """ % (float(length),float(length),float(length),float(length),float(length),float(length),float(length))

    sp_tree = dendropy.Tree.get_from_string(sp_tree_str, "newick")
    gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=sp_tree.taxon_namespace,
        num_contained=1)
    gene_tree_list = TreeList()

    for i in range(num):
        gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree,
        gene_to_containing_taxon_map=gene_to_species_map)
        treesim.contained_coalescent_tree(containing_tree=sp_tree,
                                      gene_to_containing_taxon_map=gene_to_species_map)
        for t in gene_tree.leaf_nodes():
            t.taxon.label = t.taxon.label.split( )[0]
        gene_tree_list.append(gene_tree)

    gene_tree_list.write_to_path(fout, 'newick')


def run_command(command):
    p = subprocess.Popen(command,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')

# Local Declarations
n = 8  # number of species
e = .05  # Error
type = sys.argv[1]  # Type 1 is caterpillar and type 2 is balanced
f = sys.argv[2]  # Shortest branch length
estimate = 9 / (2 * (f ** 2)) * math.log(ncr(n, 4)/e*4)

generateCoalescentTrees(1, int(estimate/2), "a", f)

print(estimate)
for output_line in run_command('java -jar astral.4.10.5.jar -i a -o out.tre'):
    print(output_line)

