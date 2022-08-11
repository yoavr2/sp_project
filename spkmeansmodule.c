#include "spkmeans.c"
#define PY_SSIZE_T_CLEAN  /* For all # variants of unit formats (s#, y#, etc.) use Py_ssize_t rather than int. */
#include <Python.h>       /* MUST include <Python.h>, this implies inclusion of the following standard headers:
                             <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */


/*============================================================================
==============================================================================
== THIS WAS MY KMEANS IMPLEMENTATION, THE FUNCTIONS THEMSELVES NOT THE API! ==
============================ - END UP HERE - =================================
==============================================================================*/

static PyObject* GetList(double **mat, int *mat_size)
{
    int m,n, i, j;
    m = mat_size[0];
    n = mat_size[1];
    PyObject* python_mat = PyList_New(m);
    for (i = 0; i < m; ++i)
    {
        PyObject* python_row = PyList_New(n);
        for (j = 0; j < n; j++){
            double r = mat[i][j];
            PyObject* python_float = Py_BuildValue("d", r);
            PyList_SetItem(python_row, j, python_float);
        }
        PyList_SetItem(python_mat, i, python_row);
    }
    return python_mat;
}

/*
 * This actually defines the kmeans function using a wrapper C API function
 * The wrapping function needs a PyObject* self argument.
 * This is a requirement for all functions and methods in the C API.
 * It has input PyObject *args from Python.
 */
static PyObject* wam_capi(PyObject *self, PyObject *args)
{//should add the epsilon
    //at the python code we should send the max_iter or the default value if we don't get one from the user
    char *file_name;



    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "s", &file_name)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

/* This builds the answer ("d" = Convert a C double to a Python floating point number) back into a python object */
    /*return Py_BuildValue("d", test_ctopy(d));*/ /*  Py_BuildValue(...) returns a PyObject*  */

    return GetList(wam(file_name), mat_size(file_name));
}





/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef capiMethods[] = {
    {"test",                   /* the Python method name that will be used */
      (PyCFunction) test_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function *///NEED TO CHANGE DESCRIPTION
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};


/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}