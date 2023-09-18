#include "fastpace_utils.h"

static PyObject* run_motif_discovery(PyObject* self, PyObject* args, PyObject* keywords) {
    PyObject* peptides_list;
    PyObject* weights_list = NULL;
    
    // Check if the function is called with the correct arguments
    if (!PyArg_ParseTuple(args, "O!|O", &PyList_Type, &peptides_list, &weights_list)) {
        PyErr_SetString(PyExc_Exception, "This function takes a list of sequences in the first argument, and optionally a list of weights in the second argument");
        return NULL;
    }

    // Check if peptides_list is empty or less than 2 sequences
    Py_ssize_t peptides_num = PyList_Size(peptides_list);
    if (peptides_num < 2) {
        PyErr_SetString(PyExc_Exception, "The list of sequences has less than 2 sequences");
        return NULL;
    }

    // Check if peptides_list is a list of strings
    for (Py_ssize_t i = 0; i < peptides_num; i++) {
        PyObject *item = PyList_GetItem(peptides_list, i);
        if (!PyUnicode_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "First argument must be a list of strings");
            return NULL;
        }
    }

    // If weights_list is provided
    if (weights_list != NULL && weights_list != Py_None) {
        // Check if weights_list is a list
        if (!PyList_Check(weights_list)) {
            PyErr_SetString(PyExc_TypeError, "Second argument must be a list of numbers");
            return NULL;
        }

        // Check if weights_list has the same length as peptides_list
        Py_ssize_t num_weights = PyList_Size(weights_list);
        if (num_weights != peptides_num) {
            PyErr_SetString(PyExc_Exception, "The number of weights must be equal to the number of sequences");
            return NULL;
        }

        // Check if weights_list is a list of numbers (doubles or integers)
        for (Py_ssize_t i = 0; i < num_weights; i++) {
            PyObject *item = PyList_GetItem(weights_list, i);
            if (!PyNumber_Check(item)) {
                PyErr_SetString(PyExc_TypeError, "Second argument must be a list of numbers");
                return NULL;
            }
        }
    }

    Dataset dataset = parse_dataset(peptides_list, weights_list);

    // Check if error occured while parsing the dataset
    if (dataset.peptides_num == 0) {
        return NULL;
    }

    IterativeSimilarityScoresResult similarity_scores_result = calculate_iterative_similarity_scores(dataset);
    MotifsResult motifs_result = extract_putative_motifs(dataset, similarity_scores_result.peptides_scores);
    int best_residue_score_peptide_indx = get_best_residue_score_peptide_indx(similarity_scores_result.peptides_scores, dataset.peptides_lengths, dataset.peptides_num);
    AlignmentResult alignment_result = align_dataset_to_peptide(dataset, similarity_scores_result.peptides_scores, best_residue_score_peptide_indx);
    PyObject* result_dict = create_result_dict(dataset, similarity_scores_result, motifs_result, alignment_result, peptides_list);

    // Cleanup
    free_peptides_scores(similarity_scores_result.peptides_scores, dataset.peptides_num);
    free_dataset(dataset);
    free_motifs_result(motifs_result);
    free_alignment_result(alignment_result);
    
    return result_dict;
}

static char run_motif_discovery_docs[] = "run_motif_discovery(sequences, weights): Calculate per residue similarity scores, extract putative motifs.\n";

static PyMethodDef module_funcs[] = {
    {
        "run_motif_discovery",
        (PyCFunction)run_motif_discovery,
        METH_VARARGS | METH_KEYWORDS,
        run_motif_discovery_docs
    },
    {NULL}
};

static char module_doc[] = "C extension module implementing Brute Force Approach to do Pattern Matching between sequences, and generate pssm-like similarity matrices for each peptide. These matrices are used in extracting motifs from the sequences. \n";

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
