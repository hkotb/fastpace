## FaSTPACE
FaSTPACE (Fast and Scalable Tool for Peptide Alignment and Consensus Extraction) is an innovative tool designed for peptide alignment and motif extraction. FaSTPACE generates a refined global similarity matrix, which serves as the basis from which these outputs are produced.

## Description

| ![FaSTPACE algorithm](https://raw.githubusercontent.com/hkotb/fastpace/master/img/algorithm.png) |
|:--:|
| **Figure 1** Schema of the steps in the FaSTPACE tool. |

The core of the algorithm produces a refined global similarity matrix from which these outputs are produced. A global similarity matrix is a probabilistic representation of the similarity of the peptide to all other peptides in the dataset. The method consists of three steps: initiation, refinement and post-processing of the global similarity matrix.

## Getting Started

### Enter peptides

| ![input](https://raw.githubusercontent.com/hkotb/fastpace/master/img/input.png) |
|:--:|
| **Figure 2** The Input Peptides Interface on the FaSTPACE Web Server. |

Load peptides or click on one of the example buttons to load a test set. Please note that input peptides must not contain non-standard amino acids.

##### Input Formats: You can enter peptides in two different formats:

- List of Peptides: Each peptide should be on a separate line. Optionally, you can include weights separated by a space.
- FASTA Format: You can also input peptides in the standard FASTA format.

### Submit
Click the "Submit" button and wait for the queue system to process your job.

### Results
After the job is completed, you will receive detailed information about the best-extracted motif and alignment, including:

- Motif Section: This section displays the following:
  
    <table align="center">
    <thead>
    <tr>
    <th align="center"><a target="_blank" rel="noopener noreferrer nofollow" href="https://raw.githubusercontent.com/hkotb/fastpace/master/img/motif.png"><img align="center" src="https://raw.githubusercontent.com/hkotb/fastpace/master/img/motif.png"></a></th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="center"><strong>Figure 3</strong> The Motif Section Interface within the FaSTPACE Web Server.</td>
    </tr>
    </tbody>
    </table>

    - Enriched Motif: Regular expression of the best-p-value motif extracted from a peptide's similarity matrix. The allowable characters for building the regular expression are listed below:
    <table align="center">
    <thead>
    <tr>
    <th align="center">Character</th><th align="center">Name</th><th align="center">Meaning</th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td>D, E, K, R, H, S, T, C, M, N, Q, F, Y, W, G, A, V, L, I, P</td><td>residue</td><td>One amino acid</td>
    </tr>
    <tr>
    <td>.</td><td>dot</td><td>Any amino acid is allowed</td>
    </tr>
    <tr>
    <td>[...]</td><td>character class</td><td>One amino acid listed is allowed</td>
    </tr>
    </tbody>
    </table>
    
    - Best Extracted Motif Scores: 
        - Sig: The significance score of the motif with the highest p-value in the dataset.
        - P-value: The highest p-value of a motif in the dataset.
        - Instances: Number of peptides matching the motif.
        - Coverage:  Number of matches Percentage in the total number of peptides in the dataset.
    - Sequence Logo: Displaying the peptide with the highest residue score in its similarity matrix. It can be downloaded in SVG and PNG formats.

- Downloads Section: This section provides the results in downloadable formats:

    <table align="center">
    <thead>
    <tr>
    <th align="center"><a target="_blank" rel="noopener noreferrer nofollow" href="https://raw.githubusercontent.com/hkotb/fastpace/master/img/downloads.png"><img align="center" src="https://raw.githubusercontent.com/hkotb/fastpace/master/img/downloads.png"></a></th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="center"><strong>Figure 4</strong> The Downloads Section Interface within the FaSTPACE Web Server.</td>
    </tr>
    </tbody>
    </table>

  - peptides JSON: Downloads a JSON file containing JSON object for each peptide with its similarity matrix and all calculated scores based on the motif extracted from its similarity matrix, matching to the best motif, and aligning to the template peptide.
  - consensus JSON: Downloads a JSON file containing a JSON object of the best motif extracted from the dataset and its calculated scores.
  - alignment JSON: Downloads a file of the alignment in JSON format.
  - alignment FASTA: Downloads a file of the alignment in FASTA format.

- Alignment Section: This section displays a table with the following columns:
  
    <table align="center">
    <thead>
    <tr>
    <th align="center"><a target="_blank" rel="noopener noreferrer nofollow" href="https://raw.githubusercontent.com/hkotb/fastpace/master/img/alignment.png"><img align="center" src="https://raw.githubusercontent.com/hkotb/fastpace/master/img/alignment.png"></a></th>
    </tr>
    </thead>
    <tbody>
    <tr>
    <td align="center"><strong>Figure 5</strong> The Alignment Section Interface within the FaSTPACE Web Server.</td>
    </tr>
    </tbody>
    </table>

    - Peptide: The alignment of the peptide to the template peptide using their global similarity matrices.
    - Alignment Score: Alignment score calculated by aligning the peptide to the template peptide using their global similarity matrices.
    - Motif Similarity Score: The p-value of the motif extracted from the peptid's global similarity matrix.
    - Motif: Regular expression of the motif extracted from the peptide's global similarity matrix.
    - Coverage: Number of peptides matching the motif extracted from this peptide (Number of matches Percentage in the total number of peptides in the dataset).

## Technology Stack
- Backend: Our backend is powered by the [Python](https://www.python.org/) programming language, with the main functionality implemented in C/C++. We utilise [SciPy](https://scipy.org/) C++ binomial distribution files as part of our open-source Python package. Additionally, we leverage [FastAPI](https://fastapi.tiangolo.com/) for our APIs, and internally, [Apache Kafka](https://kafka.apache.org/) is used for data import from internal pipelines.
- Frontend: The user interface is built using [React](https://react.dev/).

## Contact Information and Feedback
If you require technical assistance or would like to provide feedback on your experience with the website, please feel free to contact us at norman.davey@icr.ac.uk. We value your feedback and are here to assist you.
