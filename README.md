This repository contains scripts for the following manuscript: Chotai, Wei and Messer (2024): https://www.biorxiv.org/content/10.1101/2024.07.26.605365v1.

Directory: `model` contains SLiM and msprime scripts for the simulation model
1. For sweeps in continuous space:
* `model/spatial_WF_grid_neutral.slim`: Neutral burn-in simulation script. Must specify a dispersal distance for parent selection.
* `model/spatial_WF_grid_neutral_panmixia.slim`: Neutral burn-in simulation script. Limited to simulating a panmictic population (uses all individuals during parent selection).
   * For additional speed-up, `reproduction()` block updated to select all parents beforehand.
* `model/spatial_WF_grid_sweep_from_neutral.slim`: Sweep simulation that runs till fixation or loss. If provided, it begins from neutral burn-in tree sequence. Must specify a dispersal distance for parent selection.
    * For better logging, updated to output sweep status directly into file.
* `model/spatial_WF_grid_sweep_from_neutral_panmixia.slim`: Sweep simulation that runs till fixation or loss. If provided, it begins from neutral burn-in tree sequence. Limited to simulating panmictic populations (uses all individuals during parent selection).
    * For additional speed-up, `reproduction()` block updated to select all parents beforehand.
    * For better logging, updated to output sweep status directly into file.
* `model/spatial_WF_grid_sweep_extend.slim`: Neutral simulation that takes in fixed sweep tree-sequence files (from `model/spatial_WF_grid_sweep_from_neutral.slim`) and runs the simulation for $10^3$ generations after fixation, periodically saving intermediate tree-sequence files. Must specify a dispersal distance for parent selection.
* `model/spatial_WF_grid_sweep_panmixia_extend.slim`: Neutral simulation that takes in fixed sweep tree-sequence files (from `model/spatial_WF_grid_sweep_from_neutral_panmixia.slim`) and runs the simulation for $10^3$ generations after fixation, periodically saving intermediate tree-sequence files. Limited to simulating panmictic populations (uses all individuals during parent selection).
    *  For additional speed-up, `reproduction()` block updated to select all parents beforehand.

2. For sweeps from standing genetic variation:
* `model/msprime_SGV_to_slim.py`: Simulates a population with mutation using msprime. Randomly picks a mutation present at the specified starting frequency and removes all other mutations. Picked mutation is assigned specified selection coefficient. Saves simulation as tree sequence.
* `model/msprime_SGV_to_slim.slim`: Loads in tree sequence from `model/msprime_SGV_to_slim.py` and runs a sweep simulation till fixation or loss.

Directory: `sampling` contains a python script for sampling and overlaying neutral mutations.
* `sampling/sampling_strategies_SLiMMut.py`: Reads in tree sequence file. Use the command-line arguments to choose from two types of sampling: `-m` is local sampling and `-r` is global sampling (as demonstrated in Figure 1B in the manuscript). Simplifies tree to include chosen samples and overlays neutral mutations. Use the command-line arguments to choose from two types of output: `--vcf` for a `.vcf` file or `--tree` for a tree sequence file.
* `sampling/sampling_strategies_SLiMMut_extend.py`: Identical to `sampling/sampling_strategies_SLiMMut.py`, but with better support to read output from `model/spatial_WF_grid_sweep_extend.slim`. 

Directory: `statistics` contains python scripts that use `tskit` and `scikit-allel` for computing statistics on tree sequence files overlaid with neutral mutations.
* `statistics/tmrca.py`: writes out TMRCA between randomly-chosen samples for each tree. [Figure 2A]
* `statistics/neutral_pi.py`: writes out pairwise $\pi$ across the whole genome. Uses [`ts.diversity`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.diversity). [Figure 2B]
*  `statistics/neutral_pi_per_genome.py`: writes out $\pi$ *within* a diploid genome within non-overlapping windows. [Figure 2B]
* `statistics/sweep_size.py`: writes out pairwise $\pi$ over non-overlapping windows within a specified region of the genome. Uses [`ts.diversity`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.diversity). [Figure 5A; Figure S7]
    * In order to decrease variance, updated to use `mode: 'branch'`
* `statistics/sfs_statistics.py`: writes out SFS and Tajima's D statistics within specified windows (whole-genome for neutral simulations).  Uses [`ts.allele_frequency_spectrum`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.allele_frequency_spectrum) and [`ts.Tajimas_D`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.Tajimas_D). [Figures 6A and 6B; Figure S5; Figure S9A]
* `statistics/haplotype_frequency_spectra.py`: writes out values for haplotype frequency spectra. Uses [`allel.HaplotypeArray.distinct_counts`](https://scikit-allel.readthedocs.io/en/stable/model/ndarray.html#allel.HaplotypeArray.distinct_counts). [Figure 6C; Figure S6]
* `statistics/haplotype_statistics.py`: writes out haplotype heterozygosity and number of haplotypes. Computation adapted from source code for [`allel.HaplotypeArray.distinct_counts`](https://scikit-allel.readthedocs.io/en/stable/model/ndarray.html#allel.HaplotypeArray.distinct_counts). [Figures 6D and 6E; Figure 7; Figure S9B]
* `statistics/mutation_age_SGV.py`: writes the age of the SGV mutation.
* Sub-directory: `statistics/supplement_statistics` contains python scripts for computing statistics to generate supplementary figures.
    * Sub-directory `statistics/supplement_statistics/euclidean_dist`: runs a SLiM model for one generation and outputs location of parental and offspring individuals, and calculates the Euclidean distance between them. [Figure S1]
    * Subdirectory `statistics/supplement_statistics/garudh_stats`: runs [`allel.HaplotypeArray.moving_garud_h`](https://scikit-allel.readthedocs.io/en/stable/stats/selection.html#allel.moving_garud_h). [Figure S8]
    * Subdirectory `statistics/supplement_statistics/neutral_stats`: runs $F_{ST}$ and PCA analyses. [Figure S3]

Directory: `plot` contains python scripts for generating figures and supplement figures in the manuscript and a zipped file `plot/data.zip` with raw-data files needed to generate the figures.
