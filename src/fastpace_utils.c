#include "fastpace_utils.h"
#include <regex.h>
#include "scipy_interface.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// CONSTANTS /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

const int BLOSUM62PSM[25][25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
             F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */ 
    /*A*/  { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1,
            -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4       },
    /*R*/  {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1,
            -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4       },
    /*N*/  {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2,
            -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4       },
    /*D*/  {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3,
            -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4       },
    /*C*/  { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1,
            -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4       },
    /*Q*/  {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0,
            -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4       },
    /*E*/  {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2,
            -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4       },
    /*G*/  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3,
            -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4       },
    /*H*/  {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2,
            -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4       },
    /*I*/  {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,
             0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4       },
    /*L*/  {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,
             0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4       },
    /*K*/  {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1,
            -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4       },
    /*M*/  {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,
             0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4       },
    /*F*/  {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,
             6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4       },
    /*P*/  {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2,
            -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4       },
    /*S*/  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1,
            -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4       },
    /*T*/  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,
            -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4       },
    /*W*/  {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,
             1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4       },
    /*Y*/  {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,
             3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4       },
    /*V*/  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1,
            -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4       },
    /*B*/  {-2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3,
            -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4       },
    /*J*/  {-1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,
             0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4       },
    /*Z*/  {-1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1,
            -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4       },
    /*X*/  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4       },
    /***/  {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1       }
};

const char LETTERS_LIST[25] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*'};

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// FUNCTIONS /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int get_letter_index(char C) {
    switch(C)
    {
        case 'A':
            return 0;
        case 'R':
            return 1;
        case 'N':
            return 2;
        case 'D':
            return 3;
        case 'C':
            return 4;
        case 'Q':
            return 5;
        case 'E':
            return 6;
        case 'G':
            return 7;
        case 'H':
            return 8;
        case 'I':
            return 9;
        case 'L':
            return 10;
        case 'K':
            return 11;
        case 'M':
            return 12;
        case 'F':
            return 13;
        case 'P':
            return 14;
        case 'S':
            return 15;
        case 'T':
            return 16;
        case 'W':
            return 17;
        case 'Y':
            return 18;
        case 'V':
            return 19;
        case 'B':
            return 20;
        case 'J':
            return 21;
        case 'Z':
            return 22;
        case 'X':
            return 23;
        case '*':
            return 24;
    }
    return -1;
}

Dataset parse_dataset(PyObject* peptides_list, PyObject* weights_list) {
    Py_ssize_t peptides_num = PyList_Size(peptides_list);
    
    const char** peptides_strs = malloc(sizeof(char*)*peptides_num);
    int* peptides_lengths = malloc(sizeof(int)*peptides_num);
    double* peptides_weights = malloc(sizeof(double)*peptides_num);
    double total_weights = 0;
    int maximum_score = 0;
    int error = 0;
    for (int i = 0; i < peptides_num; i++)
    {
        PyObject* peptide = PyList_GetItem(peptides_list, i);
        peptides_strs[i] = PyUnicode_AsUTF8(peptide);
        peptides_lengths[i] = strlen(peptides_strs[i]);
        
        if (peptides_lengths[i] < 2) {
            PyErr_SetString(PyExc_Exception, "Found a sequence with less than 2 letters");
            error = 1;
            break;
        }
        if (peptides_lengths[i] > 32) {
            PyErr_SetString(PyExc_Exception, "FaSTPACE is not tested on sequences longer than 32 characters.\nIf it is required to increase this limit, the static arrays (regex_str, pos_str) in the extract_putative_motifs function should be increased accordingly or changed to dynamic arrays as their maximum size now is 256 chars.\nHowever, this will violate the SLiM definition, so it is recommended to tweak the algorithm by setting a maximum motif length to make it suitable for longer sequences.");
            error = 1;
            break;
        }
        if (weights_list != NULL && weights_list != Py_None)
        {
            PyObject* peptide_weight = PyList_GetItem(weights_list, i);
            peptides_weights[i] = PyFloat_AsDouble(peptide_weight);
            if (peptides_weights[i] <= 0)
            {
                PyErr_SetString(PyExc_Exception, "Found a non-positive weight");
                error = 1;
                break;
            }
        }
        else
        {
            peptides_weights[i] = 1;
        }
        total_weights+=peptides_weights[i];
        
        int peptide_score = 0;
        for (int j = 0; j < peptides_lengths[i]; j++)
        {
            char x = peptides_strs[i][j];
            
            int x_indx = get_letter_index(x);
            
            if (x_indx == -1) {
                // Handle the case of a non-standard sequence letter
                PyErr_SetString(PyExc_Exception, "Found a non-standard sequence letter");
                error = 1;
                break;
            }
            
            peptide_score += BLOSUM62PSM[x_indx][x_indx];
        }
        if (peptide_score > maximum_score)
        {
            maximum_score = peptide_score;
        }
        if (error)
        {
            break;
        }
    }

    if (error)
    {
        free(peptides_strs);
        free(peptides_lengths);
        free(peptides_weights);
        Dataset dataset = {
            .peptides_num = 0,
            .peptides_strs = NULL,
            .peptides_lengths = NULL,
            .peptides_weights = NULL,
            .total_weights = 0,
            .maximum_score = 0
        };
        return dataset;
    }

    Dataset dataset = {
        .peptides_num = peptides_num,
        .peptides_strs = peptides_strs,
        .peptides_lengths = peptides_lengths,
        .peptides_weights = peptides_weights,
        .total_weights = total_weights,
        .maximum_score = maximum_score
    };

    return dataset;
}

double calculate_residue_score_using_blosum(double**  prev_sA_scores, double**  prev_sB_scores, int sA_letter_indx, int sB_letter_indx, int sA_sequence_indx, int sB_sequence_indx) {
    double residue_weight = BLOSUM62PSM[sA_letter_indx][sB_letter_indx];
    return residue_weight;
}

double calculate_residue_score_using_sequences_scores(double**  prev_sA_scores, double**  prev_sB_scores, int sA_letter_indx, int sB_letter_indx, int sA_sequence_indx, int sB_sequence_indx) {
    double residue_weight = prev_sA_scores[sB_letter_indx][sA_sequence_indx]+prev_sB_scores[sA_letter_indx][sB_sequence_indx];
    return residue_weight;
}

MatchResult match_sequences(const char* sA, const char* sB, int n, int m, double**  prev_sA_scores, double**  prev_sB_scores, ResidueScoreFunc residue_score_func) {
    int start_gap = (m > 2) ? (m - 1) : 0;
    int found_match = 0;
    
    // Allocate memory for the result
    double *sA_scores = malloc(sizeof(double) * n);
    double *sB_scores = malloc(sizeof(double) * m);
    char *sA_matched_chars = malloc(sizeof(char) * n);
    char *sB_matched_chars = malloc(sizeof(char) * m);
    
    // Initialize the result arrays
    memset(sA_scores, 0, sizeof(*sA_scores) * n);
    memset(sB_scores, 0, sizeof(*sB_scores) * m);
    
    // Initialize variables for match tracking
    //for general use, it will not cause problems between different iterations because its values will be overwritten and accessed through total_match variable.
    int match_indx[n+m];
    char match_char[n+m];
    
    double best_match_score = 0;
    int best_alignment_start = 0; // Align to the start of the template sequence by default if no match is found
    
    for (int i = 0; i < n+2*start_gap-m+1; ++i) {
        int total_match = 0;
        double match_score = 0;
        for (int j = 0; j < m; ++j) {
            if (i+j < start_gap) {
                continue;
            }
            if (i+j >= n+start_gap) {
                break;
            }
            char x_i = sA[i+j-start_gap];
            char y_j = sB[j];
            
            int x_i_indx = get_letter_index(x_i);
            int y_i_indx = get_letter_index(y_j);
            
            if (x_i_indx == -1 || y_i_indx == -1) {
                // Handle the case of a non-standard sequence letter
                PyErr_SetString(PyExc_Exception, "Found a non-standard sequence letter");
                MatchResult result = { .found_match = 0, .sA_scores = NULL, .sB_scores = NULL, .sA_matched_chars = NULL, .sB_matched_chars = NULL };
                return result;
            }
            
            double residue_score = residue_score_func(prev_sA_scores, prev_sB_scores, x_i_indx, y_i_indx, i+j-start_gap, j);
            
            if (residue_score > 0) {
                match_indx[total_match] = i+j-start_gap;
                match_indx[n+total_match] = j;
                match_char[total_match] = y_j;
                match_char[n+total_match] = x_i;
                total_match++;
                match_score += residue_score;
            }
        }
        if (total_match > 1 && match_score > best_match_score) {
            found_match = 1;
            best_match_score = match_score;
            best_alignment_start = i-start_gap;
            
            memset(sA_scores, 0, sizeof(*sA_scores) * n);
            memset(sB_scores, 0, sizeof(*sB_scores) * m);
            
            for (int j = 0; j < total_match; ++j) {
                if (match_score > sA_scores[match_indx[j]]) {
                    sA_scores[match_indx[j]] = match_score;
                    sA_matched_chars[match_indx[j]] = match_char[j];
                }
                if (match_score > sB_scores[match_indx[n+j]]) {
                    sB_scores[match_indx[n+j]] = match_score;
                    sB_matched_chars[match_indx[n+j]] = match_char[n+j];
                }
            }
        }
    }

    // Create a MatchResultDouble struct and assign the result arrays
    MatchResult result = {
        .found_match = found_match,
        .match_score = best_match_score,
        .alignment_start = best_alignment_start,
        .sA_scores = sA_scores,
        .sB_scores = sB_scores,
        .sA_matched_chars = sA_matched_chars,
        .sB_matched_chars = sB_matched_chars
    };
    
    return result;
}

void free_match_result(MatchResult* match) {
    free(match->sA_scores);
    free(match->sB_scores);
    free(match->sA_matched_chars);
    free(match->sB_matched_chars);
}

double* calculate_similarity_pvals(Dataset dataset) {
    double* pvals = calloc(dataset.maximum_score + 1, sizeof(double));

    for (int i = 0; i < dataset.peptides_num; i++) {
        const char* peptide1_str = dataset.peptides_strs[i];
        int peptide1_len = dataset.peptides_lengths[i];

        for (int j = i + 1; j < dataset.peptides_num; j++) {
            const char* peptide2_str = dataset.peptides_strs[j];
            int peptide2_len = dataset.peptides_lengths[j];

            MatchResult match = match_sequences(peptide1_str, peptide2_str, peptide1_len, peptide2_len, NULL, NULL, calculate_residue_score_using_blosum);

            for (int k = 0; k < peptide1_len; ++k) {
                int res_score = match.sA_scores[k];
                pvals[res_score] += (dataset.peptides_weights[j] * dataset.total_weights) / (dataset.total_weights - dataset.peptides_weights[i]);
            }

            for (int k = 0; k < peptide2_len; ++k) {
                int res_score = match.sB_scores[k];
                pvals[res_score] += (dataset.peptides_weights[i] * dataset.total_weights) / (dataset.total_weights - dataset.peptides_weights[j]);
            }

            // Cleanup
            free_match_result(&match);
        }
    }

    double pval = 0;
    for (int s = dataset.maximum_score; s >= 0; s--) {
        pval += pvals[s];
        pvals[s] = pval;
    }

    if (pval > 0) {
        for (int s = dataset.maximum_score; s >= 0; s--) {
            pvals[s] = -1 * log(pvals[s] / pvals[0]);
        }
    }

    return pvals;
}

double*** calculate_similarity_scores(Dataset dataset, double* pvals, double*** prev_peptides_scores, int iteration) {
    double*** peptides_scores = malloc(sizeof(double**) * dataset.peptides_num);
    for (int i = 0; i < dataset.peptides_num; i++) {
        peptides_scores[i] = malloc(sizeof(double*) * 25); // 25 is the number of allowed letters
        for (int l = 0; l < 25; ++l) {
            peptides_scores[i][l] = malloc(sizeof(double) * dataset.peptides_lengths[i]);
            for (int j = 0; j < dataset.peptides_lengths[i]; j++) {
                peptides_scores[i][l][j] = 0;
            }
        }
    }

    for (int i = 0; i < dataset.peptides_num; i++) {
        const char* peptide1_str = dataset.peptides_strs[i];
        int peptide1_len = dataset.peptides_lengths[i];
        double** peptide1_scores = peptides_scores[i];

        for (int j = i + 1; j < dataset.peptides_num; j++) {
            const char* peptide2_str = dataset.peptides_strs[j];
            int peptide2_len = dataset.peptides_lengths[j];

            MatchResult match;
            if (iteration == 0) {
                match = match_sequences(peptide1_str, peptide2_str, peptide1_len, peptide2_len, NULL, NULL, calculate_residue_score_using_blosum);
            } else {
                match = match_sequences(peptide1_str, peptide2_str, peptide1_len, peptide2_len, prev_peptides_scores[i], prev_peptides_scores[j], calculate_residue_score_using_sequences_scores);
            }

            if (match.found_match > 0) {
                for (int k = 0; k < peptide1_len; ++k) {
                    void* res_score;
                    int positive_score_flag = 0;
                    if (iteration == 0) {
                        res_score = malloc(sizeof(int));
                        *(int*)res_score = match.sA_scores[k];
                        if (*(int*)res_score > 0) {
                            positive_score_flag = 1;
                        }
                    } else {
                        res_score = malloc(sizeof(double));
                        *(double*)res_score = match.sA_scores[k];
                        if (*(double*)res_score > 0) {
                            positive_score_flag = 1;
                        }
                    }
                    if (positive_score_flag) {
                        char res_match = match.sA_matched_chars[k];
                        int l = get_letter_index(res_match);

                        if (iteration == 0) {
                            peptide1_scores[l][k] += (dataset.peptides_weights[j] * dataset.total_weights * pvals[*(int*)res_score]) / (dataset.total_weights - dataset.peptides_weights[i]);
                        } else {
                            peptide1_scores[l][k] += (dataset.peptides_weights[j] * dataset.total_weights * *(double*)res_score) / (dataset.total_weights - dataset.peptides_weights[i]);
                        }
                    }
                    free(res_score);
                }

                double** peptide2_scores = peptides_scores[j];
                for (int k = 0; k < peptide2_len; ++k) {
                    void* res_score;
                    int positive_score_flag = 0;
                    if (iteration == 0) {
                        res_score = malloc(sizeof(int));
                        *(int*)res_score = match.sB_scores[k];
                        if (*(int*)res_score > 0) {
                            positive_score_flag = 1;
                        }
                    } else {
                        res_score = malloc(sizeof(double));
                        *(double*)res_score = match.sB_scores[k];
                        if (*(double*)res_score > 0) {
                            positive_score_flag = 1;
                        }
                    }
                    if (positive_score_flag) {
                        char res_match = match.sB_matched_chars[k];
                        int l = get_letter_index(res_match);

                        if (iteration == 0) {
                            peptide2_scores[l][k] += (dataset.peptides_weights[i] * dataset.total_weights * pvals[*(int*)res_score]) / (dataset.total_weights - dataset.peptides_weights[j]);
                        } else {
                            peptide2_scores[l][k] += (dataset.peptides_weights[i] * dataset.total_weights * *(double*)res_score) / (dataset.total_weights - dataset.peptides_weights[j]);
                        }
                    }
                    free(res_score);
                }
            }
            // Cleanup
            free_match_result(&match);
        }
    }

    return peptides_scores;
}

double calculate_frobenius_norm(double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num) {
    double sum_of_sq = 0;
    for (int i = 0; i < peptides_num; i++) {
        for (int l = 0; l < 25; ++l) {
            for (int j = 0; j < peptides_lengths[i]; j++) {
                sum_of_sq += (peptides_scores[i][l][j] * peptides_scores[i][l][j]);
            }
        }
    }
    return sqrt(sum_of_sq);
}

void normalize_scores(double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num) {
    double frobenius_norm = calculate_frobenius_norm(peptides_scores, peptides_lengths, peptides_num);

    if (frobenius_norm == 0) {
        return;
    }

    for (int i = 0; i < peptides_num; i++) {
        for (int l = 0; l < 25; ++l) {
            for (int j = 0; j < peptides_lengths[i]; j++) {
                peptides_scores[i][l][j] /= frobenius_norm;
            }
        }
    }
}

int check_convergence(double*** new_peptides_scores, double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num) {
    double sum_of_diff = 0;
    int sum_of_lengths = 0;
    for (int i = 0; i < peptides_num; i++) {
        int peptide_length = peptides_lengths[i];
        for (int l = 0; l < 25; ++l) {
            for (int j = 0; j < peptide_length; j++) {
                sum_of_diff += fabs(new_peptides_scores[i][l][j] - peptides_scores[i][l][j]);
            }
        }
        sum_of_lengths += peptide_length;
    }
    double convergence_threshold = 0.000001 * sum_of_lengths;
    return (sum_of_diff < convergence_threshold) ? 1 : 0;
}

IterativeSimilarityScoresResult calculate_iterative_similarity_scores(Dataset dataset) {
    double* similarity_pvals = calculate_similarity_pvals(dataset);
    double*** peptides_scores = calculate_similarity_scores(dataset, similarity_pvals, NULL, 0);
    free(similarity_pvals);
    
    normalize_scores(peptides_scores, dataset.peptides_lengths, dataset.peptides_num);
    int converged = 0;
    int itr = 0;
    do
    {
        itr++;
        double*** new_peptides_scores = calculate_similarity_scores(dataset, NULL, peptides_scores, itr);
        normalize_scores(new_peptides_scores, dataset.peptides_lengths, dataset.peptides_num);
        
        converged = check_convergence(new_peptides_scores, peptides_scores, dataset.peptides_lengths, dataset.peptides_num);
        
        for (int i = 0; i < dataset.peptides_num; i++)
        {
            for (int l = 0; l < 25; ++l)
            {
                free(peptides_scores[i][l]);
            }
            free(peptides_scores[i]);
        }
        free(peptides_scores);
        
        peptides_scores = new_peptides_scores;
    }
    while (!converged && itr < 100);

    IterativeSimilarityScoresResult result = {
        .peptides_scores = peptides_scores,
        .iterations = itr
    };

    return result;
}

int compare_scores(const void* a, const void* b) {
    const AAPosScr* obj1 = (const AAPosScr*)a;
    const AAPosScr* obj2 = (const AAPosScr*)b;
    double diff = obj2->scr - obj1->scr; // Reverse order
    if (diff < 0) {
        return -1;
    } else if (diff > 0) {
        return 1;
    } else {
        return 0;
    }
}

MotifsResult extract_putative_motifs(Dataset dataset, double*** peptides_scores) {
    // Calculate amino acid frequencies
    double aa_freq[25] = {0};
    int aa_total = 0;
    double wt_aa_freq[25] = {0};
    double wt_aa_total = 0;
    for (int i = 0; i < dataset.peptides_num; i++) {
        const char* sA = dataset.peptides_strs[i];
        for (int j = 0; j < dataset.peptides_lengths[i]; j++) {
            char x = sA[j];
            int x_indx = get_letter_index(x);
            aa_freq[x_indx] += 1;
            aa_total += 1;
            wt_aa_freq[x_indx] += dataset.peptides_weights[i];
            wt_aa_total += dataset.peptides_weights[i];
        }
    }
    // Normalize amino acid frequencies
    for (int l = 0; l < 25; ++l) {
        aa_freq[l] /= aa_total;
        wt_aa_freq[l] /= wt_aa_total;
    }
    
    // Allocate arrays to store results
    // Similarity motif for a peptide is extracted from its similarity matrix
    char** similarity_motifs = malloc(sizeof(char*)*dataset.peptides_num);
    double* similarity_p_vals = malloc(sizeof(double)*dataset.peptides_num);
    double* similarity_significances = malloc(sizeof(double)*dataset.peptides_num);
    // Matched motif for a peptide is the best similarity motif (from any peptide) that matches this peptide
    char** matched_motifs = malloc(sizeof(char*)*dataset.peptides_num);
    double* matched_p_vals = malloc(sizeof(double)*dataset.peptides_num);
    double* matched_significances = malloc(sizeof(double)*dataset.peptides_num);
    
    // Initialize matched_motifs and matched_p_vals to NULL and 1 respectively
    for (int i = 0; i < dataset.peptides_num; i++) {
        matched_motifs[i] = NULL;
        matched_p_vals[i] = 1;
        matched_significances[i] = 1;
    }

    // Variables for tracking the best motif
    int found_motif = 0;
    int best_motif_indx = -1;
    double best_motif_p_val = 1;
    double best_motif_significance = 1;
    
    // Calculate a conservative size of the motif space
    double B_L = 0;
    for (int i = 0; i < dataset.peptides_num; i++) {
        B_L += dataset.peptides_lengths[i]*20;
    }

    // Loop over all peptides
    for (int i = 0; i < dataset.peptides_num; i++) {
        // Initialize a string to store the regex motif
        char* motif = malloc(sizeof(char));
        motif[0] = '\0';
        // Start with the worst possible p-value
        double p_val = 1;
        double significance = 1;
        // Initialize an array to store the matched peptides flags for the best motif
        int* matched_flag = calloc(dataset.peptides_num, sizeof(int));
        // Get the peptide length and scores
        double** peptide_scores = peptides_scores[i];
        int peptide_length = dataset.peptides_lengths[i];
        // estimate the maximum number of entries in peptide_scores matrix to allocate memory for aa_pos_scr
        int max_entries = peptide_length*25;
        AAPosScr* aa_pos_scr = (AAPosScr*) malloc(sizeof(AAPosScr)*max_entries);
        // Loop over all positions in the peptide matrix and the amino acids at each position and store non-zero scores in aa_pos_scr while keeping track of the number of entries 
        int num_entries = 0;
        for (int j = 0; j < peptide_length; j++) {
            for (int l = 0; l < 25; ++l) {
                if (peptide_scores[l][j] > 0) {
                    aa_pos_scr[num_entries].aa = l;
                    aa_pos_scr[num_entries].pos = j;
                    aa_pos_scr[num_entries].scr = peptide_scores[l][j];
                    num_entries++;
                }
            }
        }

        // If the number of entries in aa_pos_scr is less than 2, skip the peptide
        if (num_entries < 2) {
            similarity_motifs[i] = motif;
            similarity_p_vals[i] = p_val;
            similarity_significances[i] = significance;
            continue;
        }

        // Sort aa_pos_scr by score in descending order up to num_entries to avoid uninitialized places in the array
        qsort(aa_pos_scr, num_entries, sizeof(AAPosScr), compare_scores);
        
        // Initialize a regex table matching the peptide matrix size to zeros
        int** regex_table = malloc(sizeof(int*)*25);
        for (int l = 0; l < 25; ++l) {
            regex_table[l] = calloc(peptide_length, sizeof(int));
        }
        // Update the regex table with the first entry in aa_pos_scr to start the motif
        regex_table[aa_pos_scr[0].aa][aa_pos_scr[0].pos] = 1;
        // Initialize an array to store the number of characters in each position in the regex table
        int* pos_chars_counts = calloc(peptide_length, sizeof(int));
        // Update the number of characters in the first entry's position in the regex table to 1
        pos_chars_counts[aa_pos_scr[0].pos]++;
        // Loop over all the possible number of motif components starting from 2 to the number of entries in aa_pos_scr (the maximum possible number of motif components)
        // In each iteration, add the next entry in aa_pos_scr to the regex table and update the motif and check if the p-value of the motif is improved over the previous iteration if yes, aprrove the motif, if not, stop the loop
        for (int n_residues = 2; n_residues <= num_entries; n_residues++) {
            // Update the regex table with the next entry in aa_pos_scr
            regex_table[aa_pos_scr[n_residues-1].aa][aa_pos_scr[n_residues-1].pos] = 1;
            // Update the number of characters in the next entry's position in the regex table
            pos_chars_counts[aa_pos_scr[n_residues-1].pos]++;
            // Start forming the motif
            char regex_str[256] = ".*";
            int actual_motif_length = 0;
            int motif_start = -1;
            int motif_end = -1;
            // Variable to hold the probability of the motif sequence
            double p_m = 1;
            double wt_p_m = 1;
            // Count the number of positions in the motif and update the motif start and end positions
            for (int j = 0; j < peptide_length; j++) {
                if (pos_chars_counts[j] > 0)
                {
                    if (motif_start < 0) {
                        motif_start = j;
                    }
                    motif_end = j;
                    actual_motif_length++;
                }
            }
            // If the actual motif length (number of positions) is less than 2, go to the next iteration
            if (actual_motif_length < 2) {
                continue;
            }
            // Form the regex string
            for (int j = motif_start; j < motif_end+1; j++) {
                if (pos_chars_counts[j] > 0) {
                    double p_i = 0;
                    double wt_p_i = 0;

                    char pos_str[256] = "";
                    for (int l = 0; l < 25; ++l) {
                        if (regex_table[l][j] > 0) {
                            char chStr[2];  // Create a string to hold the character
                            chStr[0] = LETTERS_LIST[l];
                            chStr[1] = '\0';  // Add the null-terminating character
                            strcat(pos_str, chStr);
                            p_i += aa_freq[l];
                            wt_p_i += wt_aa_freq[l];
                        }
                    }
                    if (pos_chars_counts[j] == 1) {
                        strcat(regex_str, pos_str);
                    }
                    else {
                        strcat(regex_str, "[");
                        strcat(regex_str, pos_str);
                        strcat(regex_str, "]");
                    }
                    p_m *= p_i;
                    wt_p_m *= wt_p_i;
                }
                else {
                    strcat(regex_str, ".");
                }
            }
            strcat(regex_str, ".*");
            
            // Calculate the number of sequences with the motif and flag them in a temporary array
            regex_t rgx;
            int ret;
            int num_sequences_with_motif = 0;
            ret = regcomp(&rgx, regex_str, REG_EXTENDED);
            if (ret) {
                printf("Failed to compile regex %s\n", regex_str);
                continue;
            }
            int* matched_flag_temp = malloc(sizeof(int)*dataset.peptides_num);
            for (int ii = 0; ii < dataset.peptides_num; ii++) {
                ret = regexec(&rgx, dataset.peptides_strs[ii], 0, NULL, 0);
                if (ret == 0) {
                    matched_flag_temp[ii] = 1;
                    num_sequences_with_motif++;
                } else {
                    matched_flag_temp[ii] = 0;
                }
            }
            regfree(&rgx);

            // Calculate the probability of the motif in the dataset
            int total_motif_length = motif_end+1-motif_start;
            double N_m = peptide_length-total_motif_length+1;
            double p_1_plus = 1-pow(1-p_m, N_m);
            double wt_p_1_plus = 1-pow(1-wt_p_m, N_m);

            double survival = 1;
            double sig = 1;
            // If the number of sequences with the motif is greater than 1 (to avoid self-matching), calculate the p-value of the motif
            if (num_sequences_with_motif > 1)
            {
                survival = calculate_binomial_pmf_plus_survival(num_sequences_with_motif, dataset.peptides_num, p_1_plus);
                double wt_survival = calculate_binomial_pmf_plus_survival(num_sequences_with_motif, dataset.peptides_num, wt_p_1_plus);
                sig = 1-pow(1-wt_survival,B_L);
            }
            // If the p-value of the motif is less than the p-value of the previous iteration, update the p-value and the motif
            // Check on p-value is less than 1 with at least 10-e6 to avoid floating point errors (legacy check, not sure if it is still needed)
            if (1-survival > 0.000001 && survival < p_val) {
                p_val = survival;
                significance = sig;
                free(motif);
                motif = malloc((strlen(regex_str) + 1) * sizeof(char));
                strcpy(motif, regex_str);
                free(matched_flag);
                matched_flag = matched_flag_temp;
            // If the p-value of the motif is not less than the p-value of the previous iteration, break the loop
            } else {
                free(matched_flag_temp);
                break;
            }
        }
        
        // Store the motif and the p-value of the motif in the results arrays
        similarity_motifs[i] = malloc((strlen(motif) + 1) * sizeof(char));
        strcpy(similarity_motifs[i], motif);
        //free(motif); // Do not free the motif here because it is used in the matched_motifs array
        similarity_p_vals[i] = p_val;
        similarity_significances[i] = significance;
        
        // Update the matched motif data for the peptides that match the motif and the motif improves their p-value of the previous matched motif
        for (int ii = 0; ii < dataset.peptides_num; ii++) {
            if (matched_flag[ii] == 1 && p_val < matched_p_vals[ii]) {
                if (matched_motifs[ii] != NULL) {
                    free(matched_motifs[ii]);
                }
                matched_motifs[ii] = malloc((strlen(motif) + 1) * sizeof(char));
                strcpy(matched_motifs[ii], motif);
                matched_p_vals[ii] = p_val;
                matched_significances[ii] = significance;
            }
        }
        free(motif);

        // If the p-value of the motif is less than the p-value of the best motif, update the best motif data
        if (p_val < best_motif_p_val) {
            found_motif = 1;
            best_motif_indx = i;
            best_motif_p_val = p_val;
            best_motif_significance = significance;
        }
        
        // Cleanup
        for (int l = 0; l < 25; ++l) {
            free(regex_table[l]);
        }
        free(matched_flag);
        free(regex_table);
        free(pos_chars_counts);
        free(aa_pos_scr);
    }

    // Create the result object
    MotifsResult motifs_result = {
        .found_motif = found_motif,
        .best_motif_indx = best_motif_indx,
        .best_motif_p_val = best_motif_p_val,
        .best_motif_significance = best_motif_significance,
        .similarity_motifs = similarity_motifs,
        .similarity_p_vals = similarity_p_vals,
        .similarity_significances = similarity_significances,
        .matched_motifs = matched_motifs,
        .matched_p_vals = matched_p_vals,
        .matched_significances = matched_significances
    };

    return motifs_result;
}

int get_best_residue_score_peptide_indx(double*** peptides_scores, int* peptides_lengths, Py_ssize_t peptides_num) {
    double best_residue_score = 0;
    int best_residue_score_peptide_indx = -1;
    for (int i = 0; i < peptides_num; i++) {
        for (int l = 0; l < 25; ++l) {
            for (int j = 0; j < peptides_lengths[i]; j++) {
                if (peptides_scores[i][l][j] > best_residue_score) {
                    best_residue_score = peptides_scores[i][l][j];
                    best_residue_score_peptide_indx = i;
                }
            }
        }
    }
    return best_residue_score_peptide_indx;
}

AlignmentResult align_dataset_to_peptide(Dataset dataset, double*** peptides_scores, int peptide_indx) {
    // Check if the peptide index is valid
    if (peptide_indx < 0 || peptide_indx >= dataset.peptides_num) {
        AlignmentResult result = { .peptide_indx = -1, .min_best_align_start = -1, .max_best_align_start = -1, .best_alignment_starts = NULL, .best_alignment_scores = NULL };
        return result;
    }

    int min_best_align_start = INT_MAX;
    int max_best_align_start = INT_MIN;
    int* best_alignment_starts = malloc(sizeof(int)*dataset.peptides_num);
    double* best_alignment_scores = malloc(sizeof(double)*dataset.peptides_num);

    const char* peptide1_str = dataset.peptides_strs[peptide_indx];
    int peptide1_len = dataset.peptides_lengths[peptide_indx];
    double** peptide1_scores = peptides_scores[peptide_indx];
    for (int i = 0; i < dataset.peptides_num; i++) {
        const char* peptide2_str = dataset.peptides_strs[i];
        int peptide2_len = dataset.peptides_lengths[i];
        double** peptide2_scores = peptides_scores[i];
        MatchResult match = match_sequences(peptide1_str, peptide2_str, peptide1_len, peptide2_len, peptide1_scores, peptide2_scores, calculate_residue_score_using_sequences_scores);
        int best_alignment_start = match.alignment_start;
        best_alignment_starts[i] = best_alignment_start;
        best_alignment_scores[i] = match.match_score;
        if (best_alignment_start < min_best_align_start) {
            min_best_align_start = best_alignment_start;
        }
        if (best_alignment_start > max_best_align_start) {
            max_best_align_start = best_alignment_start;
        }

        // Cleanup
        free_match_result(&match);
    }
    
    AlignmentResult alignment_result = {
        .peptide_indx = peptide_indx,
        .min_best_align_start = min_best_align_start,
        .max_best_align_start = max_best_align_start,
        .best_alignment_starts = best_alignment_starts,
        .best_alignment_scores = best_alignment_scores
    };

    return alignment_result;
}

void fill_letters_objects(PyObject** letters_objects) {
    for (int i = 0; i < 25; ++i) {
        letters_objects[i] = PyUnicode_FromStringAndSize(&LETTERS_LIST[i], 1);
    }
}

void vacate_letters_objects(PyObject** letters_objects) {
    for (int i = 0; i < 25; ++i) {
        Py_DECREF(letters_objects[i]);
    }
}

double get_peptide_similarity_score(const char* sequence, double** peptide_scores, int peptide_length) {
    double similarity_score = 0;
    for (int i = 0; i < peptide_length; i++) {
        char x = sequence[i];
        int x_indx = get_letter_index(x);
        similarity_score += peptide_scores[x_indx][i];
    }
    return similarity_score;
}

char* generate_align_string_for_peptide(const char* sequence, int sequence_length, int best_align_start, int min_best_align_start, int max_best_align_start) {
    // Calculate the new length after aligning with adding dashes
    int new_length = sequence_length - min_best_align_start + max_best_align_start;
    
    // Allocate memory for the new aligned sequence
    char* new_sequence = (char*)malloc((new_length + 1) * sizeof(char));
    
    if (new_sequence == NULL) {
        printf("Memory allocation failed\n");
        return NULL;
    }
    
    // Add '-' characters to the left of the string
    memset(new_sequence, '-', best_align_start - min_best_align_start);
    
    // Copy the original sequence to the new aligned sequence
    memcpy(new_sequence + best_align_start - min_best_align_start, sequence, sequence_length);
    
    // Add '-' characters to the right of the string
    memset(new_sequence + best_align_start - min_best_align_start + sequence_length, '-', new_length - (best_align_start - min_best_align_start + sequence_length));
    
    // Null-terminate the new sequence
    new_sequence[new_length] = '\0';
    
    return new_sequence;
}

void set_int_item_in_dict(PyObject* dictionary, PyObject* key, int value) {
    PyObject* long_obj = PyLong_FromLong(value);
    PyDict_SetItem(dictionary, key, long_obj);
    Py_DECREF(long_obj);
}

void set_string_item_in_dict(PyObject* dictionary, PyObject* key, const char* value) {
    PyObject* str_obj = PyUnicode_FromString(value);
    PyDict_SetItem(dictionary, key, str_obj);
    Py_DECREF(str_obj);
}

void set_float_item_in_dict(PyObject* dictionary, PyObject* key, double value) {
    PyObject* float_obj = PyFloat_FromDouble(value);
    PyDict_SetItem(dictionary, key, float_obj);
    Py_DECREF(float_obj);
}

PyObject* create_result_dict(Dataset dataset, IterativeSimilarityScoresResult similarity_scores_result, MotifsResult motifs_result, AlignmentResult alignment_result, PyObject* peptides_list) {
    PyObject* result_dict = PyDict_New();
    
    PyObject* letters_objects[25];
    fill_letters_objects(letters_objects);
    PyObject* iterations_str_obj = PyUnicode_FromString("iterations");
    PyObject* best_motif_str_obj = PyUnicode_FromString("best_motif");
    PyObject* best_motif_p_val_str_obj = PyUnicode_FromString("best_motif_p_val");
    PyObject* best_motif_significance_str_obj = PyUnicode_FromString("best_motif_significance");
    PyObject* alignment_template_str_obj = PyUnicode_FromString("alignment_template");
    PyObject* peptides_str_obj = PyUnicode_FromString("peptides");
    PyObject* scores_str_obj = PyUnicode_FromString("scores");
    PyObject* similarity_motif_str_obj = PyUnicode_FromString("similarity_motif");
    PyObject* similarity_p_val_str_obj = PyUnicode_FromString("similarity_p_val");
    PyObject* similarity_significance_str_obj = PyUnicode_FromString("similarity_significance");
    PyObject* similarity_score_str_obj = PyUnicode_FromString("similarity_score");
    PyObject* matched_motif_str_obj = PyUnicode_FromString("matched_motif");
    PyObject* matched_p_val_str_obj = PyUnicode_FromString("matched_p_val");
    PyObject* matched_significance_str_obj = PyUnicode_FromString("matched_significance");
    PyObject* aligned_sequence_str_obj = PyUnicode_FromString("aligned_sequence");
    PyObject* empty_str_obj = PyUnicode_FromString("");
    
    set_int_item_in_dict(result_dict, iterations_str_obj, similarity_scores_result.iterations);
    if (motifs_result.found_motif == 1) {
        set_string_item_in_dict(result_dict, best_motif_str_obj, motifs_result.similarity_motifs[motifs_result.best_motif_indx]);
    } else {
        PyDict_SetItem(result_dict, best_motif_str_obj, empty_str_obj);
    }
    set_float_item_in_dict(result_dict, best_motif_p_val_str_obj, motifs_result.best_motif_p_val);
    set_float_item_in_dict(result_dict, best_motif_significance_str_obj, motifs_result.best_motif_significance);

    if (alignment_result.peptide_indx == -1) {
        PyDict_SetItem(result_dict, alignment_template_str_obj, empty_str_obj);
    } else {
        PyDict_SetItem(result_dict, alignment_template_str_obj, PyList_GetItem(peptides_list, alignment_result.peptide_indx));
    }

    PyObject* peptides_dict = PyDict_New();
    PyDict_SetItem(result_dict, peptides_str_obj, peptides_dict);

    for (int i = 0; i < dataset.peptides_num; i++) {
        PyObject* peptide_i = PyList_GetItem(peptides_list, i);
        PyObject* peptide_dict = PyDict_New();
        PyDict_SetItem(peptides_dict, peptide_i, peptide_dict);
        PyObject* scores_dict = PyDict_New();
        PyDict_SetItem(peptide_dict, scores_str_obj, scores_dict);
        
        int peptide_length = dataset.peptides_lengths[i];
        double** peptide_scores = similarity_scores_result.peptides_scores[i];
        
        for (int j = 0; j < peptide_length; j++) {
            PyObject* long_obj = PyLong_FromLong(j);
            PyObject* residue_dict = PyDict_New();
            PyDict_SetItem(scores_dict, long_obj, residue_dict);
            Py_DECREF(long_obj);
            for (int l = 0; l < 25; ++l) {
                if (peptide_scores[l][j] > 0) {
                    double new_val = peptide_scores[l][j];
                    PyObject * new_val_obj = PyFloat_FromDouble(new_val);
                    PyObject* letter = letters_objects[l];
                    PyDict_SetItem(residue_dict, letter, new_val_obj);
                    Py_DECREF(new_val_obj);
                }
            }
            Py_DECREF(residue_dict);
        }
        
        set_string_item_in_dict(peptide_dict, similarity_motif_str_obj, motifs_result.similarity_motifs[i]);
        free(motifs_result.similarity_motifs[i]);
        set_float_item_in_dict(peptide_dict, similarity_p_val_str_obj, motifs_result.similarity_p_vals[i]);
        set_float_item_in_dict(peptide_dict, similarity_significance_str_obj, motifs_result.similarity_significances[i]);
        
        double peptide_similarity_score = get_peptide_similarity_score(dataset.peptides_strs[i], peptide_scores, peptide_length);
        set_float_item_in_dict(peptide_dict, similarity_score_str_obj, peptide_similarity_score);

        if (motifs_result.matched_motifs[i] != NULL) {
            set_string_item_in_dict(peptide_dict, matched_motif_str_obj, motifs_result.matched_motifs[i]);
            free(motifs_result.matched_motifs[i]);
        } else {
            PyDict_SetItem(peptide_dict, matched_motif_str_obj, empty_str_obj);
        }
        set_float_item_in_dict(peptide_dict, matched_p_val_str_obj, motifs_result.matched_p_vals[i]);
        set_float_item_in_dict(peptide_dict, matched_significance_str_obj, motifs_result.matched_significances[i]);

        if (alignment_result.peptide_indx == -1) {
            PyDict_SetItem(peptide_dict, aligned_sequence_str_obj, empty_str_obj);
        } else {
            char* aligned_sequence = generate_align_string_for_peptide(dataset.peptides_strs[i], dataset.peptides_lengths[i], alignment_result.best_alignment_starts[i], alignment_result.min_best_align_start, alignment_result.max_best_align_start);
            set_string_item_in_dict(peptide_dict, aligned_sequence_str_obj, aligned_sequence);
            free(aligned_sequence);
        }

        Py_DECREF(peptide_dict);
        Py_DECREF(scores_dict);
    }
    
    // Cleanup
    vacate_letters_objects(letters_objects);
    Py_DECREF(peptides_dict);
    Py_DECREF(peptides_str_obj);
    Py_DECREF(iterations_str_obj);
    Py_DECREF(best_motif_str_obj);
    Py_DECREF(best_motif_p_val_str_obj);
    Py_DECREF(best_motif_significance_str_obj);
    Py_DECREF(alignment_template_str_obj);
    Py_DECREF(scores_str_obj);
    Py_DECREF(similarity_motif_str_obj);
    Py_DECREF(similarity_p_val_str_obj);
    Py_DECREF(similarity_significance_str_obj);
    Py_DECREF(similarity_score_str_obj);
    Py_DECREF(matched_motif_str_obj);
    Py_DECREF(matched_p_val_str_obj);
    Py_DECREF(matched_significance_str_obj);
    Py_DECREF(aligned_sequence_str_obj);
    Py_DECREF(empty_str_obj);

    return result_dict;
}

void free_peptides_scores(double*** peptides_scores, Py_ssize_t peptides_num) {
    for (int i = 0; i < peptides_num; i++) {
        for (int l = 0; l < 25; ++l) {
            free(peptides_scores[i][l]);
        }
        free(peptides_scores[i]);
    }
    free(peptides_scores);
}

void free_dataset(Dataset dataset) {
    free(dataset.peptides_strs);
    free(dataset.peptides_lengths);
    free(dataset.peptides_weights);
}

void free_motifs_result(MotifsResult motifs_result) {
    free(motifs_result.similarity_motifs);
    free(motifs_result.similarity_p_vals);
    free(motifs_result.similarity_significances);
    free(motifs_result.matched_motifs);
    free(motifs_result.matched_p_vals);
    free(motifs_result.matched_significances);
}

void free_alignment_result(AlignmentResult alignment_result) {
    if (alignment_result.peptide_indx == -1) {
        return;
    }
    free(alignment_result.best_alignment_starts);
    free(alignment_result.best_alignment_scores);
}
