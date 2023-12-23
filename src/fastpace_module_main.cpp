#include "fastpace_utils.hpp"

static int validate_peptides_list(PyObject* peptides_list) {
    // Check if peptides_list is empty or less than 2 peptides
    Py_ssize_t peptides_num = PyList_Size(peptides_list);
    if (peptides_num < 2) {
        PyErr_SetString(PyExc_Exception, "The list of peptides has less than 2 peptides");
        return 1;
    }

    // Check if peptides_list is a list of strings
    for (Py_ssize_t i = 0; i < peptides_num; i++) {
        PyObject *item = PyList_GetItem(peptides_list, i);
        if (!PyUnicode_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "First argument must be a list of strings");
            return 1;
        }
    }

    return 0;
}

static int validate_weights_list(PyObject* weights_list, Py_ssize_t peptides_num) {
    // Check if weights_list is a list
    if (!PyList_Check(weights_list)) {
        PyErr_SetString(PyExc_TypeError, "Second argument must be a list of numbers");
        return 1;
    }

    // Check if weights_list has the same length as peptides_list
    Py_ssize_t num_weights = PyList_Size(weights_list);
    if (num_weights != peptides_num) {
        PyErr_SetString(PyExc_Exception, "The number of weights must be equal to the number of peptides");
        return 1;
    }

    // Check if weights_list is a list of numbers (doubles or integers)
    for (Py_ssize_t i = 0; i < num_weights; i++) {
        PyObject *item = PyList_GetItem(weights_list, i);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "Second argument must be a list of numbers");
            return 1;
        }
    }

    return 0;
}

static PyObject* get_result_dict(PyObject* peptides_list, PyObject* weights_list, double* original_similarity_pvals, AAFreq original_aa_freq, int refine, int normalization_factor){
    Dataset dataset = parse_dataset(peptides_list, weights_list);

    // Check if error occured while parsing the dataset
    if (dataset.peptides_num == 0) {
        return NULL;
    }

    double* similarity_pvals;
    if (original_similarity_pvals == NULL) {
        similarity_pvals = calculate_similarity_pvals(dataset);
    } else {
        similarity_pvals = original_similarity_pvals;
    }
    IterativeSimilarityScoresResult similarity_scores_result = calculate_iterative_similarity_scores(dataset, similarity_pvals, refine);
    AAFreq aa_freq;
    if (original_aa_freq.aa_freq == NULL) {
        aa_freq = calculate_aa_freq(dataset);
    } else {
        aa_freq = original_aa_freq;
    }
    MotifsResult motifs_result = extract_putative_motifs(dataset, similarity_scores_result.peptides_scores, aa_freq.aa_freq, aa_freq.wt_aa_freq, normalization_factor);
    int best_residue_score_peptide_indx = get_best_residue_score_peptide_indx(similarity_scores_result.peptides_scores, dataset.peptides_lengths, dataset.peptides_num);
    AlignmentResult alignment_result = align_dataset_to_peptide(dataset, similarity_scores_result.peptides_scores, best_residue_score_peptide_indx);
    PyObject* result_dict = create_result_dict(dataset, similarity_scores_result, motifs_result, alignment_result, peptides_list);

    // Cleanup
    free(similarity_pvals);
    free_peptides_scores(similarity_scores_result.peptides_scores, dataset.peptides_num);
    free_dataset(dataset);
    free_motifs_result(motifs_result);
    free_alignment_result(alignment_result);
    free_aa_freq(aa_freq);
    
    return result_dict;
}

static PyObject* run_motif_discovery(PyObject* self, PyObject* args, PyObject* keywords) {
    PyObject* peptides_list;
    PyObject* weights_list = NULL;
    int refine = 1;
    int normalization_factor = -1;
    
    // Check if the function is called with the correct arguments
    static char *kwlist[] = {(char*) "peptides", (char*) "weights", (char*) "refine", (char*) "normalization_factor", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O!|Opi", kwlist, &PyList_Type, &peptides_list, &weights_list, &refine, &normalization_factor)) {
        PyErr_SetString(PyExc_Exception, "This function takes a list of peptides, and optionally a list of weights, refine flag, and a normalization factor");
        return NULL;
    }

    // Validate peptides_list
    if (validate_peptides_list(peptides_list)) {
        return NULL;
    }

    // If weights_list is provided
    if (weights_list != NULL && weights_list != Py_None) {
        Py_ssize_t peptides_num = PyList_Size(peptides_list);
        // Validate weights_list
        if (validate_weights_list(weights_list, peptides_num)) {
            return NULL;
        }
    }

    PyObject* result_dict = get_result_dict(peptides_list, weights_list, NULL, (AAFreq){NULL, NULL}, refine, normalization_factor);

    return result_dict;
}

static PyObject* rerun_motif_discovery(PyObject* self, PyObject* args, PyObject* keywords) {
    PyObject* original_peptides_list;
    PyObject* masked_peptides_list;
    PyObject* weights_list = NULL;
    
    // Check if the function is called with the correct arguments
    static char *kwlist[] = {(char*) "original_peptides", (char*) "masked_peptides", (char*) "weights", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O!O!|O", kwlist, &PyList_Type, &original_peptides_list, &PyList_Type, &masked_peptides_list, &weights_list)) {
        PyErr_SetString(PyExc_Exception, "This function takes a list of the original peptides, a list of the masked peptides, and optionally a list of weights");
        return NULL;
    }

    // Validate original_peptides_list
    if (validate_peptides_list(original_peptides_list)) {
        return NULL;
    }

    // Validate masked_peptides_list
    if (validate_peptides_list(masked_peptides_list)) {
        return NULL;
    }

    // Check if original_peptides_list and masked_peptides_list have the same length
    Py_ssize_t original_peptides_num = PyList_Size(original_peptides_list);
    Py_ssize_t masked_peptides_num = PyList_Size(masked_peptides_list);
    if (original_peptides_num != masked_peptides_num) {
        PyErr_SetString(PyExc_Exception, "The number of original peptides must be equal to the number of masked peptides");
        return NULL;
    }

    // If weights_list is provided
    if (weights_list != NULL && weights_list != Py_None) {
        Py_ssize_t peptides_num = PyList_Size(original_peptides_list);
        // Validate weights_list
        if (validate_weights_list(weights_list, peptides_num)) {
            return NULL;
        }
    }

    Dataset dataset = parse_dataset(original_peptides_list, weights_list);

    // Check if error occured while parsing the dataset
    if (dataset.peptides_num == 0) {
        return NULL;
    }

    double* similarity_pvals = calculate_similarity_pvals(dataset);
    AAFreq aa_freq = calculate_aa_freq(dataset);
    PyObject* result_dict = get_result_dict(masked_peptides_list, weights_list, similarity_pvals, aa_freq, 1, -1);
    
    return result_dict;
}

static char run_motif_discovery_docs[] = "run_motif_discovery(peptides, weights, refine, normalization_factor): Calculate per residue similarity scores, align peptides and extract putative motifs.\n";
static char rerun_motif_discovery_docs[] = "rerun_motif_discovery(original_peptides, masked_peptides, weights): Calculate per residue similarity scores, extract putative motifs after masking previously extracted motifs.\n";

static PyMethodDef module_funcs[] = {
    {
        "run_motif_discovery",
        (PyCFunction)run_motif_discovery,
        METH_VARARGS | METH_KEYWORDS,
        run_motif_discovery_docs
    },
    {
        "rerun_motif_discovery",
        (PyCFunction)rerun_motif_discovery,
        METH_VARARGS | METH_KEYWORDS,
        rerun_motif_discovery_docs
    },
    {NULL}
};

static char module_doc[] = "C extension module implementing Brute Force Approach to do Pattern Matching between peptides, and generate global similarity matrices for each peptide. These matrices are used in producing alignments and extracting motifs from the peptides. \n";

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "fastpace",
    module_doc,
    -1,
    module_funcs
};

PyMODINIT_FUNC PyInit_fastpace(void) {
    return PyModule_Create(&moduledef);
}
