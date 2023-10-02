# Benchmarking Metrics
All the code for evaluating the  statistics required to benchmark the model:

## Sampling Statistics

We would like these to be as low as possible

- **Reconstruction Accuracy** We measure how often the model can completely reconstruct a given molecule from its latent embedding
- **Validity** Percentage of chemically valid compounds.
- **Uniqueness** Percentage of unique compounds in the generated output
  
## Property Statistics

We would like these to be as low as possible

- **QED Score** stands for quantitative estimation of drug-likeness. The empirical rationale of the QED measure reflects the underlying distribution of molecular properties including molecular weight, logP, topological polar surface area, number of hydrogen bond donors and acceptors, the number of aromatic rings and rotatable bonds, and the presence of unwanted chemical functionalities.
- **LogP** The penalized logP score measures the solubility and synthetic accessibility of a compound.
- **Molecular Weight** It's the average molecular weight for all the molecules with different isotopic compositions present in a compound, each weighted for the natural abundance of the respective isotope or isotopes.
- **Synthetic Accessibility** We assign a score between 1 (easy to make) and 10 (very difficult to make) depending on the ease of development of a compound.

## Structural Statistics

We would like these to be as high as possible

- **Structural Nearest Neighbour**: Defines how similar the ECFP  fingerprints of the generated molecules are to the input molecules
- **Fragment Similarity**: How similar the fragments (defined by the BRICS algorithm) of the output molecules are to the input molecules
- **Scaffold Similarity**: How similar the scaffolds of the output molecules are to the input molecules
