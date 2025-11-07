# default libs
import numpy as np
import pandas as pd
import subprocess
import argparse

# tskit libs
import tskit, msprime, pyslim


docstring = \
r"""
- Produce Burn-in population with msprime (fully neutral)
- Subsequent SLiM simulation of overdominant selection that mimics stabilizing selection
"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser(docstring)
    parser.add_argument(
            "--seed", type=int, default=973,
            help="Random seed",
            )
    parser.add_argument(
            "--num_individuals", type=int, default=1e4,
            help="Number of individuals",
            )
    parser.add_argument(
            "--seq_length", type=int, default=1,
            help="Number of loci",
            )
    parser.add_argument(
            "--mu", type=float, default=1e-2,
            help="Mutation rate",
            )
    parser.add_argument(
            "--prefix", type=str, default="",
            help="Output prefix",
            )
    parser.add_argument(
            "--out_path", type=str, default="",
            help="Output file",
            )
    parser.add_argument(
            "--id_check_path", type=str,
            help="Mutation id check script",
            )
    parser.add_argument(
            "--biallelic_path", type=str,
            help="Main simulation script",
            )
    parser.add_argument(
            "--num_generations", type=int, default=100,
            help="Number of generations in slim",
            )
    parser.add_argument(
            "--frequency_path", type=str,
            help="Path to frequency file",
            )

    args = parser.parse_args()
    
    # burn-in
    ts = msprime.sim_ancestry(
            samples = args.num_individuals,
            sequence_length = args.seq_length,
            population_size = args.num_individuals,
            random_seed = args.seed,
            )
    ts = msprime.sim_mutations(
            ts, 
            rate = args.mu,
            model = msprime.BinaryMutationModel(),
            random_seed = args.seed,
            )
    assert ts.num_mutations > 1, "No mutations in tree sequence."

    print(ts.genotype_matrix().mean(axis=1).item())

    # metadata surgery
    t = ts.dump_tables() # pull tables from ts
    pyslim.annotate_tables(t, model_type="WF", tick=1) # slim annotation
    muts = t.mutations.copy() # pull out mutations
    t.mutations.clear() # purge mutations from table

    for m, tm in zip(ts.mutations(), muts):
        ds, md = tm.derived_state, tm.metadata
        if m.derived_state == "0":
            ds, md = "", {"mutation_list": []}
        else:
            ds = "0"
            md["mutation_list"][0]["mutation_type"] = 1
        t.mutations.append(tm.replace(derived_state=ds, metadata=md))

    # write
    ts_out = t.tree_sequence()
    ts_out.dump(args.prefix + args.out_path)

    # run slim - check mutation type correct
    slim_args_0 = [
            "mamba", "run", "-n", "tskit_dev",
            "slim",
            "-d", f"sequenceLength={args.seq_length}",
            "-d", f"numIndividuals={args.num_individuals}",
            "-d", f"inTreeSequence='{args.prefix + args.out_path}'",
            f"{args.id_check_path}"
            ]
    #subprocess.run(slim_args_0)

    # run slim - biallelic simulation check
    slim_args_1 = [
            "mamba", "run", "-n", "tskit_dev",
            "slim",
            "-d", f"sequenceLength={args.seq_length}",
            "-d", f"numIndividuals={args.num_individuals}",
            "-d", f"mu={args.mu}",
            "-d", f"numTicks={args.num_generations}",
            "-d", f"inTreeSequence='{args.prefix + args.out_path}'",
            "-d", f"outTreeSequence='{args.prefix + args.out_path}'",
            "-d", f"freqFile='{args.prefix + args.frequency_path}'",
            f"{args.biallelic_path}"
            ]
    subprocess.run(slim_args_1)

    # check frequency match
    ts = tskit.load(args.prefix + args.out_path)
    freq_ts = ts.genotype_matrix().mean(axis=1)
    freq_slim = pd.read_csv(args.prefix + args.frequency_path, index_col=0)["freq_mean"]
    f1, f2 = freq_ts.item(), freq_slim.item()
    print(f1, f2, f1==f2)
