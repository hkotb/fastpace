# FaSTPACE: A Fast and Scalable Tool for Peptide Alignment and Consensus Extraction

FaSTPACE (**Fa**st and **S**calable **T**ool for **P**eptide **A**lignment and **C**onsensus **E**xtraction) is a Python package implemented in C programming language.

## Table of contents
* [Description](#description)
* [Usage](#usage)
* [Development](#Development)
* [License](#license)
* [Reference](#reference)

## Description

![FaSTPACE algorithm](https://raw.githubusercontent.com/hkotb/fastpace/master/img/algorithm.png)

The core of the algorithm produces a refined global similarity matrix from which these outputs are produced. A global similarity matrix is a probabilistic representation of the similarity of the peptide to all other peptides in the dataset. The method consists of three steps: initiation, refinement and post-processing of the global similarity matrix.

## Usage

Install the pachage:
```
pip install fastpace
```

Import functions:
```
from fastpace import run_motif_discovery, rerun_motif_discovery
```

### `run_motif_discovery(peptides, weights=None, refine=1, normalization_factor=-1)`

Calculate per residue similarity scores, align peptides, and extract putative motifs.

**Parameters:**
- `peptides` (list): List of strings of standard sequence letters representing peptides in the dataset. It should contain 2 or more peptides. No sequence with less than 2 letters or longer than 32 characters is allowed.
- `weights` (list, optional, default=None): List of positive numbers representing weights. The number of weights must be equal to the number of peptides. If not provided, each peptide is assigned a weight of 1.
- `refine` (int, optional, default=1): Flag with 0 or 1. 0 means return the results after the initiation step. 1 means doing the refinement step.
- `normalization_factor` (float, optional, default=-1): The dataset size is corrected to be equal to this number in order to compare different datasets with different sizes. It must be equal to or greater than 1. By default, -1 means to use the dataset size with no normalization.

**Returns:**
- Dictionary: The dictionary object returned contains information about the discovered motif, alignment results, and similarity scores for each peptide. Below is a breakdown of the key elements in the output:
```
    .
    ├── Refinement Information:
    │   └── refinement_iterations: The number of refinement iterations performed.
    ├── Consensus Motif Information:
    │   ├── best_motif: The best consensus motif represented as a regular expression.
    │   ├── best_motif_p_val: P-value associated with the best consensus motif.
    │   ├── best_motif_significance: Significance level of the best consensus motif.
    │   ├── best_motif_num_matches: Number of matches found for the best consensus motif.
    │   └── best_motif_coverage: Coverage percentage of the best consensus motif.
    ├── Alignment Information: 
    │   ├── template: The template sequence used for alignment.
    │   └── aligned_sequences: A dictionary mapping each input sequence to its aligned counterpart.
    └── Peptide-specific Information:
        For each peptide in the dataset, the following information is provided:
        ├──  similarity_matrix: A matrix of similarity scores between the peptide and the consensus motif.
        ├──  similarity_motif: The motif represented as a regular expression for similarity.
        ├──  similarity_p_val: P-value associated with the similarity motif.
        ├──  similarity_significance: Significance level of the similarity motif.
        ├──  similarity_num_matches: Number of matches found for the similarity motif.
        ├──  similarity_coverage: Coverage percentage of the similarity motif.
        ├──  similarity_score: Cumulative similarity score for the peptide.
        ├──  matched_motif: The motif that matches the peptide best.
        ├──  matched_p_val: P-value associated with the matched motif.
        ├──  matched_significance: Significance level of the matched motif.
        ├──  matched_num_matches: Number of matches found for the matched motif.
        ├──  matched_coverage: Coverage percentage of the matched motif.
        └──  alignment_score: Cumulative alignment score for the peptide.
```
**Example:**
```python
peptides = ['ABC', 'DEF', 'GHI']
weights = [1, 2, 1]
motifs = run_motif_discovery(peptides, weights, refine=1, normalization_factor=100)
```

### `rerun_motif_discovery(original_peptides, masked_peptides, weights=None)`

Calculate per residue similarity scores and extract putative motifs after masking previously extracted motifs.

**Parameters:**
- `original_peptides` (list): List of the original peptides before masking.
- `masked_peptides` (list): List of the masked peptides. The number of masked peptides must be equal to the number of original peptides.
- `weights` (list, optional, default=None): List of positive weights.

**Returns:**
- List: List of putative motifs extracted from the aligned peptides.

**Example:**
```python
peptides = ['ABC', 'DEF', 'GHI']
weights = [1, 2, 1]
motifs = run_motif_discovery(peptides, weights, refine=1, normalization_factor=100)
```

## Development

If you clone the repository, we recommend that you install it as a development package with:
```
python setup.py develop
```

To generate distribution archive from your machine:
```
python -m build
```

To distribute binary Python extensions as wheels on Linux:
```
docker pull quay.io/pypa/manylinux2014_x86_64
docker run --rm -e PLAT=manylinux2014_x86_64 -v `pwd`:/io quay.io/pypa/manylinux2014_x86_64  /io/build-wheels.sh
```

## License

This source code is licensed under the MIT license found in the `LICENSE` file in the root directory of this source tree.

## Reference

If you find the pipeline useful in your research, we ask that you cite our paper:
```
Coming soon.
```

