#ifndef FASTPACEUTILS_H
#define FASTPACEUTILS_H

#define PY_SSIZE_T_CLEAN
#include <Python.h>

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// GLOBAL CONSTANTS //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

extern const int BLOSUM62PSM[25][25];
extern const char LETTERS_LIST[25];

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// STRUCTURE DEFINITIONS /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// Structure representing the peptides dataset
typedef struct {
    Py_ssize_t peptides_num;
    const char** peptides_strs;
    int* peptides_lengths;
    double* peptides_weights;
    double total_weights;
    int maximum_score;
} Dataset;

// Structure representing the result of matching two sequences
typedef struct {
    int found_match;
    double match_score;
    double *sA_scores;
    double *sB_scores;
    int alignment_start;
    char *sA_matched_chars;
    char *sB_matched_chars;
} MatchResult;

// Structure representing the result of iterative similarity scores calculation
typedef struct {
    double*** peptides_scores;
    int iterations;
} IterativeSimilarityScoresResult;

// Structure representing a cell in the pssm-like scoring matrix of a peptide sequence (aminoAcidIndex-position-score)
typedef struct {
    int aa;
    int pos;
    double scr;
} AAPosScr;

// Structure representing the result of extracting motifs from a dataset of peptides based on their similarity scores
typedef struct {
    int found_motif;
    int best_motif_indx;
    double best_motif_p_val;
    double best_motif_significance;
    char** similarity_motifs;
    double* similarity_p_vals;
    double* similarity_significances;
    char** matched_motifs;
    double* matched_p_vals;
    double* matched_significances;
} MotifsResult;

// Structure representing the result of aligning a dataset of peptides to a peptide
typedef struct {
    int peptide_indx;
    int min_best_align_start;
    int max_best_align_start;
    int* best_alignment_starts;
    double* best_alignment_scores;
} AlignmentResult;

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// FUNCTION DECLARATIONS /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// Function to get the index of a letter in the letters list
int get_letter_index(char C);

// Function to parse the dataset from Python objects to C structures
Dataset parse_dataset(PyObject* peptides_list, PyObject* weights_list);

// Function pointer type for the residue weighting function which can have either of the two implementations (calculate_residue_score_using_blosum or calculate_residue_score_using_sequences_scores)
typedef double (*ResidueScoreFunc)(double**, double**, int, int, int, int);

// Function to calculate a residue score in one sequence when it matches another residue in another sequence using the BLOSUM62 PSM p-values and it is used as a residue weighting function in the start of the algorithm
double calculate_residue_score_using_blosum(double**  prev_sA_scores, double**  prev_sB_scores, int sA_letter_indx, int sB_letter_indx, int sA_sequence_indx, int sB_sequence_indx);

// Function to calculate a residue score in one sequence when it matches another residue in another sequence using the sequence's previous iteration's pssm-like similarity matrix and it is used as a residue weighting function in the subsequent iterations of the algorithm
double calculate_residue_score_using_sequences_scores(double**  prev_sA_scores, double**  prev_sB_scores, int sA_letter_indx, int sB_letter_indx, int sA_sequence_indx, int sB_sequence_indx);

// Function to match two sequences and calculate per-residue similarity scores
MatchResult match_sequences(const char* sA, const char* sB, int n, int m, double**  prev_sA_scores, double**  prev_sB_scores, ResidueScoreFunc residueScoreFunc);

// Function to free the memory allocated for the MatchResult structure
void free_match_result(MatchResult* match);

// Function to calculate similarity p-values for the dataset
double* calculate_similarity_pvals(Dataset dataset);

// Function to calculate the similarity scores for the dataset in a given iteration
double*** calculate_similarity_scores(Dataset dataset, double* pvals, double*** prev_peptides_scores, int iteration);

// Function to calculate the Frobenius norm of peptide scores
double calculate_frobenius_norm(double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num);

// Function to normalize peptide scores
void normalize_scores(double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num);

// Function to check convergence between new and old peptide scores
int check_convergence(double*** new_peptides_scores, double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num);

// Function to calculate iterative similarity scores for the dataset until convergence
IterativeSimilarityScoresResult calculate_iterative_similarity_scores(Dataset dataset);

// Compare function for sorting based on score
int compare_scores(const void* a, const void* b);

// Function to extract putative motifs from the dataset based on the similarity scores
MotifsResult extract_putative_motifs(Dataset dataset, double*** peptides_scores);

// Function to get the index of the peptide with the best residue score in the pssm-like similarity matrix
int get_best_residue_score_peptide_indx(double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num);

// Function to align dataset of peptides to a peptide by aligning the pairwise  pssm-like similarity matrices of the dataset to the pssm-like similarity matrix of the peptide
AlignmentResult align_dataset_to_peptide(Dataset dataset, double*** peptides_scores, int peptide_indx);

// Function to fill the letters objects array
void fill_letters_objects(PyObject** letters_objects);

// Function to vacate the letters objects array
void vacate_letters_objects(PyObject** letters_objects);

// Function to calculate the similarity score of a peptide given its pssm-like similarity matrix
double get_peptide_similarity_score(const char* sequence, double** peptide_scores, int peptide_length);

// Function to generate an aligned string of a sequence by adding dashes
char* generate_align_string_for_peptide(const char* sequence, int sequence_length, int best_align, int min_best_align_start, int max_best_align_start);

// Function to insert C integer item in a Python dictionary
void set_int_item_in_dict(PyObject* dict, PyObject* key, int value);

// Function to insert C double item in a Python dictionary
void set_double_item_in_dict(PyObject* dict, PyObject* key, double value);

// Function to insert C string item in a Python dictionary
void set_string_item_in_dict(PyObject* dict, PyObject* key, const char* value);

// Function to create a Python dictionary object from the result of the motif discovery
PyObject* create_result_dict(Dataset dataset, IterativeSimilarityScoresResult similarity_scores_result, MotifsResult motifs_result, AlignmentResult alignment_result, PyObject* peptides_list);

// Function to free the memory allocated for the pssm-like similarity matrices of the peptides
void free_peptides_scores(double*** peptides_scores, Py_ssize_t peptides_num);

// Function to free the memory allocated for the Dataset structure
void free_dataset(Dataset dataset);

// Function to free the memory allocated for the MotifsResult structure
void free_motifs_result(MotifsResult motifs_result);

// Function to free the memory allocated for the AlignmentResult structure
void free_alignment_result(AlignmentResult alignment_result);


#endif
