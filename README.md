# KPOP_SELS
Kimura POpulation and Protein Structures Evolve via Linked Sites

**What is KPOP SELS?**

Kimura POpulation and Protein Structures Evolve via Linked Sites (KPOP SELS) is a python script to simulate evolution of a protein sequence along a phylogenetic tree. The change in protein stability of a proposed mutation determines the change in fitness of the protein and probability of fixation, which is calculated using the Kimura equation. 

Please read the following; A thermodynamic model of protein structure evolution explains empirical amino acid substitution matrices by Norn et. al. (https://doi.org/10.1002/pro.4155). The foundation of KPOP_SELS is based upon it.

**What is the motivation behind KPOP SELS?**

There are plethora of programs available to simulate evolution. However, many of these programs assume that a mutation at one site does not affect the probability of observing a mutation at another site. In our model a mutation at one site does affect the probability of observing a mutation at another site.

**How does KPOP SELS work?**

For a given proposed mutation in a protein structure we determine the probability of fixation from the Kimura equation.

![Fig_Abstract](https://user-images.githubusercontent.com/111892527/206770574-88d90850-d4bd-4a9d-bd90-aef6264e59fb.svg)

**What can we do with KPOP SELS?**

You're only limited by your imagination. 

I mainly wrote this as a way to evaluate various evolutionary models used in likelihood, Bayesian, and machine learning approaches of ancestral sequence reconstruction. 

However, we are exploring many questions.

(1) How do different parameters like effective population size, protein structure, and tree topology affect the identity of the sequences we observe, those sequence's stability, and the types of accepted mutations over the course of evolution. All of these observations are the basis of empirical models of evolution that are standard in phylogenetics. 


![4ND5_Neff_1000000_Emut_100000_ind_dGs_miyazawa](https://user-images.githubusercontent.com/111892527/206775725-55eb2415-4d0e-4a6b-a8d0-78131786e075.png)

We can track the dG of each protein as it evolves along some branch. For the above figure, we arbitrarily track the dG for 100,000 proposed mutations and only a fraction of which actually fix into the population. We can see that sequence A gets trapped in a local minima between proposed mutations 60,000-80,000.

(2) What assumptions limit our empirical models of evolution and how can those empirical models be improved?

(3) Are most mutations deleterious, neutral, or nearly neutral with respect to protein stability? How would this affect the observed stability of proteins we observe today?

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
