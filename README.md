# KPOP_SELS
Kimura POpulation and Protein Structures Evolve via Linked Sites

**What is KPOP SELS?**

Kimura POpulation and Protein Structures Evolve via Linked Sites (KPOP SELS) is a python script to simulate evolution of a protein sequence based on the interaction energy of residues in a protein structure. The change in interaction energy upon residue substitution determines the fitness of the protein, which in turn determines the selection coefficient. A probability of fixation is calculated using the Kimura equation. 

**What is the motivation behind KPOP SELS?**

There are plethora of programs available to simulate evolution. However, many of these programs assume that a mutation at one site does not effect the probability of observing a mutation at another site. In our model we are calculating the change in contact energy between the wild-type residue and the mutant, which effects the protein fitness and the probability that mutation fixes. Therefore, a mutation at one site does effect the probability of observing a mutation at another site in KPOP SELS.

**How does KPOP SELS work?**

Mutations are randomly introduced along a branch in a phylogenetic tree and the corresponding contact energies for the wild type and mutant are calculated using the contact energies from Miyazawa and Jernigan 1985. The change in interaction energy upon residue mutation determines the fitness of the protein, which in turn determines the selection coefficient, and as a result the probability of fixation. If a randomly selected number between 0 and 1 is less than the probability of fixation, then the mutation becomes fixed in the population and is the new wild-type residue. 

**What can we do with KPOP SELS?**

You're only limited by your imagination. 

I mainly wrote this as a way to evaluate various evolutionary models used in likelihood, Bayesian, and machine learning approaches of ancestral sequence reconstruction. 

However, we are exploring many questions.

(1) How do different parameters like effective population size, protein structure, and tree topology affect the identity of the sequences we observe, those sequence's stability, and the types of accepted mutations over the course of evolution. All of these observations are the basis of empirical models of evolution that are standard in phylogenetics. 

(2) What assumptions limit our empirical models of evolution and how can those empirical models be improved?

(3) Are most mutations deleterious, neutral, or nearly neutral with respect to protein stability? How would this effect the observed stability of proteins we observe today?

## Implementing

**Requirements for KPOP SELS**

Dendropy: https://dendropy.org/downloading.html
```
python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git
```

NumPy: https://numpy.org/install/
```
pip install numpy
```

Biopython: https://biopython.org/wiki/Download
```
pip install biopython
```

Matplotlib: https://matplotlib.org/stable/users/installing/index.html
```
pip install matplotlib
```

PyPDB: https://pypi.org/project/pypdb/
```
pip install pypdb
```

## Running

Copy the scripts and one of the test data sets into some test folder. 

Usage
```
./KPOP_SELS_v0.1.sh -h 
```

Example
```
./KPOP_SELS_v0.1.sh -T tree_01.txt -P 4nd5
```

## Output
A histogram of all sequence dGs post selection.
4nd5_Neff_1000000_af_dGs_.png
A histogram of all sequence ddGs post selection.
4nd5_Neff_1000000_af_ddGs_.png
A histogram of all sequence ddGs before selection.
4nd5_Neff_1000000_b4_ddGs_.png
A histogram of all sequence dGs before selection.
4nd5_Neff_1000000_b4_dGs_.png
The ancestral sequences in fasta format.
4nd5_Neff_1000000_ancs_.txt
The extant sequences in fasta format.
4nd5_Neff_1000000_extants_.txt
Histograms of the selection coefficients, probabilities of accepted mutations, fixation probabilities of accepted mutations.
4nd5_Neff_1000000_slxn_probs_fix.png
