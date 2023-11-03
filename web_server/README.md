# Getting Started

### Enter peptides
Load peptides or click on one of the example buttons to load a test set. Please note that input peptides must not contain non-standard amino acids.

### Submit
Click the "Submit" button and wait for the queue system to process your job.

### Results
After the job is completed, you will receive detailed information about the best-extracted motif and alignment, including:

- Motif Section:
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
