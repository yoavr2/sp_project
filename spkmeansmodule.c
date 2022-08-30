#include "spkmeans.c"
#define PY_SSIZE_T_CLEAN  /* For all # variants of unit formats (s#, y#, etc.) use Py_ssize_t rather than int. */
#include <Python.h>       /* MUST include <Python.h>, this implies inclusion of the following standard headers:
                             <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */


/*============================================================================
==============================================================================
== THIS WAS MY KMEANS IMPLEMENTATION, THE FUNCTIONS THEMSELVES NOT THE API! ==
============================ - END UP HERE - =================================
==============================================================================*/

static PyObject* GetMat(double **mat, int *size){
    int m,n, i, j;
    PyObject* python_mat;
    m = size[0];
    n = size[1];
    python_mat = PyList_New(m);
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

double **c_mat(PyObject *float_mat){
    PyObject *row, *item;
    double **mat;
    int rows, columns,i ,j;


    rows = (int)PyObject_Length(float_mat);
    if (rows < 0){
        return NULL;
    }
    row = PyList_GetItem(float_mat, 0);
    columns = (int)PyObject_Length(row);
    mat = (double **) malloc(sizeof(double *) * rows);
    for (i = 0; i < rows; i++){
        mat[i] = (double *) malloc(sizeof(double) * columns);
    }
    for (i = 0; i < rows; i++){
        row = PyList_GetItem(float_mat, i);
        for (j = 0; j< columns; j++){
            item = PyList_GetItem(row, j);
            if (!PyFloat_Check(item))
                mat[i][j] = 0.0;
            mat[i][j] = PyFloat_AsDouble(item);
        }
    }

    return mat;

}

int *c_mat_size(PyObject *float_mat){
    int rows, columns, *size;
    PyObject *row;

    rows = (int)PyObject_Length(float_mat);
    row = PyList_GetItem(float_mat, 0);
    columns = (int)PyObject_Length(row);

    size = (int *) malloc(sizeof(int) * 2);
    size[0] = rows;
    size[1] = columns;

    return size;
}

/*
 * This actually defines the kmeans function using a wrapper C API function
 * The wrapping function needs a PyObject* self argument.
 * This is a requirement for all functions and methods in the C API.
 * It has input PyObject *args from Python.
 */
static PyObject* wam_capi(PyObject *self, PyObject *args)
{
    double **mat;
    int *size, *res_size;
    PyObject *float_mat;

    if (!PyArg_ParseTuple(args, "O", &float_mat)){
        return NULL;
    }

    mat = c_mat(float_mat);
    size = c_mat_size(float_mat);
    res_size = (int *) malloc(sizeof(int) * 2);
    res_size[0] = size[0];
    res_size[1] = size[0];

    return GetMat(wam(mat, size[0], size[1]), res_size);
}

static PyObject* ddg_capi(PyObject *self, PyObject *args)
{
    double **mat;
    int *size, *res_size;
    PyObject *float_mat;

    if (!PyArg_ParseTuple(args, "O", &float_mat)){
        return NULL;
    }

    mat = c_mat(float_mat);
    size = c_mat_size(float_mat);
    res_size = (int *) malloc(sizeof(int) * 2);
    res_size[0] = size[0];
    res_size[1] = size[0];

    return GetMat(ddg(mat, size[0], size[1]), res_size);
}

static PyObject* lnorm_capi(PyObject *self, PyObject *args)
{
    double **mat;
    int *size, *res_size;
    PyObject *float_mat;

    if (!PyArg_ParseTuple(args, "O", &float_mat)){
        return NULL;
    }

    mat = c_mat(float_mat);
    size = c_mat_size(float_mat);
    res_size = (int *) malloc(sizeof(int) * 2);
    res_size[0] = size[0];
    res_size[1] = size[0];

    return GetMat(lnorm(mat, size[0], size[1]), res_size);
}

static PyObject* jacobi_capi(PyObject *self, PyObject *args)
{
    double **mat;
    int *size;
    PyObject *float_mat;

    if (!PyArg_ParseTuple(args, "O", &float_mat)){
        return NULL;
    }

    mat = c_mat(float_mat);
    size = c_mat_size(float_mat);
    size[0] = size[0] + 1;

    return GetMat(jacobi(mat, size[1]), size);
}

static PyObject* heuristic_capi(PyObject *self, PyObject *args){
    int k;
    double **mat;
    int *size;
    PyObject *float_mat;

    if(!PyArg_ParseTuple(args, "O", &float_mat)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    mat = c_mat(float_mat);
    size = c_mat_size(float_mat);
    k = heuristic(mat, size[0], size[1]);

    return Py_BuildValue("i", k);

}

static PyObject* spk_capi(PyObject *self, PyObject *args)
{
    double **mat, **res_mat;
    int *size, *res_size;
    int k, max_iter;
    double epsilon;
    PyObject *float_mat;

    if(!PyArg_ParseTuple(args, "iidO", &k, &max_iter, &epsilon, &float_mat)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    mat = c_mat(float_mat);
    size = c_mat_size(float_mat);

    if (k == 0){
        k = heuristic(mat, size[0], size[1]);
    }
    res_mat = spk(mat, size[0], size[1], k, max_iter, epsilon);

    res_size = (int *) malloc(sizeof(int) * 2);
    res_size[0] = k;
    res_size[1] = size[1];

    return GetMat(res_mat, res_size);

    //return Py_BuildValue("i", 0);

}


/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef capiMethods[] = {
    {"wam",                   /* the Python method name that will be used */
      (PyCFunction) wam_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function *///NEED TO CHANGE DESCRIPTION
    {"ddg",                   /* the Python method name that will be used */
      (PyCFunction) ddg_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function *///NEED TO CHANGE DESCRIPTION
    {"lnorm",                   /* the Python method name that will be used */
      (PyCFunction) lnorm_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function *///NEED TO CHANGE DESCRIPTION
    {"jacobi",                   /* the Python method name that will be used */
      (PyCFunction) jacobi_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function *///NEED TO CHANGE DESCRIPTION
    {"heuristic",                   /* the Python method name that will be used */
      (PyCFunction) heuristic_capi, /* the C-function that implements the Python function and returns static PyObject*  */
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