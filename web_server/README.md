## FaSTPACE
FaSTPACE (Fast and Scalable Tool for Peptide Alignment and Consensus Extraction) is an innovative tool designed for peptide alignment and motif extraction. FaSTPACE generates a refined global similarity matrix, which serves as the basis from which these outputs are produced.

## Description

![FaSTPACE algorithm](https://raw.githubusercontent.com/hkotb/fastpace/master/img/algorithm.png)

The core of the algorithm produces a refined global similarity matrix from which these outputs are produced. A global similarity matrix is a probabilistic representation of the similarity of the peptide to all other peptides in the dataset. The method consists of three steps: initiation, refinement and post-processing of the global similarity matrix.

## Getting Started

### Enter peptides

| ![input](https://raw.githubusercontent.com/hkotb/fastpace/master/img/input.png) |
|:--:|
| **Figure 1** The Input Peptides Interface on the FaSTPACE Web Server. |

Load peptides or click on one of the example buttons to load a test set. Please note that input peptides must not contain non-standard amino acids.

### Submit
Click the "Submit" button and wait for the queue system to process your job.

### Results
After the job is completed, you will receive detailed information about the best-extracted motif and alignment, including:

- Motif Section:
  
| ![motif](https://raw.githubusercontent.com/hkotb/fastpace/master/img/motif.png) |
|:--:|
| **Figure 2** The Motif Section Interface within the FaSTPACE Web Server. |

    - Sequence Logo: Displaying the peptide with the highest residue score in its similarity matrix.
    - Regular Expression: The best-p-value motif extracted from a peptide's similarity matrix. The allowable characters for building the regular expression are listed below:
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
        - Significance.
        - p-value.
        - Number of mateched instances.
        - Coverage.
- Alignment Section: This section displays a table with the following columns:
  
| ![alignment](https://raw.githubusercontent.com/hkotb/fastpace/master/img/alignment.png) |
|:--:|
| **Figure 2** The Alignment Section Interface within the FaSTPACE Web Server. |

    - Aligned peptide.
    - Score.
    - Score.
    - Motif.
    - Score.

## Technology Stack
- Backend: Our backend is powered by the [Python](https://www.python.org/) programming language, with the main functionality implemented in C/C++. We utilise [SciPy](https://scipy.org/) C++ binomial distribution files as part of our open-source Python package. Additionally, we leverage [FastAPI](https://fastapi.tiangolo.com/) for our APIs, and internally, [Apache Kafka](https://kafka.apache.org/) is used for data import from internal pipelines.
- Frontend: The user interface is built using [React](https://react.dev/).

## Contact Information and Feedback
If you require technical assistance or would like to provide feedback on your experience with the website, please feel free to contact us at norman.davey@icr.ac.uk. We value your feedback and are here to assist you.
