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
##### examples
Coming soon

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

