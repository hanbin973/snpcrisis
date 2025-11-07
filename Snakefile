import numpy as np

# simulation settings
## effect sizes
eff_setting = "constant"
vs = 1
## basic genetics
num_individuals, num_loci = int(1e3), 1
mu = 1 / num_individuals
num_generations = num_individuals * 20

# i/o settings
prefix_trees = "trees/"
prefix_tmp = "tmp/"
src_path = "slim_src/overdominant.slim"

# Simulate under overdominant selection
_tree = prefix_trees + "a{eff_const}.{it}.trees"
_freq = prefix_tmp + "a{eff_const}.{it}.freq"
rule tree_const:
	output:
		tree=_tree,
		freq=_freq,
	shell:
		"""
		python pysrc/01-overdominant.py \
			--seed={wildcards.it} \
			--num_individuals={num_individuals} \
			--seq_length={num_loci} \
			--mu={mu} \
			--num_generations={num_generations} \
			--vs={vs} \
			--eff_setting={eff_setting} \
			--eff_var=1 \
			--eff_const={wildcards.eff_const} \
			--src_path={src_path} \
			--out_path={output.tree} \
			--frequency_path={output.freq}
		"""

its = np.arange(500) + 1
aas = np.arange(5)
rule all:
	input:
		expand(_tree, eff_const=aas, it=its)
