#include "cPythonUtility.h"

/// Write a string to Python console
void print(char* str)
{
    PyObject* sysmod = PyImport_ImportModule("sys");
    PyObject* pyout = PyObject_GetAttrString(sysmod, "stdout");
    PyObject* result = PyObject_CallMethod(pyout, "write", "s", str);
    Py_XDECREF(result);
    Py_XDECREF(pyout);
    Py_XDECREF(sysmod);
}

/// Format a string and write it to the Python console
void print(char *format, float v1)
{
	char buffer[100];
	sprintf(buffer, format, v1);
	print(buffer);
}

void print(char *format, float v1, float v2)
{
	char buffer[100];
	sprintf(buffer, format, v1, v2);
	print(buffer);
}

void print(char *format, float v1, float v2, float v3)
{
	char buffer[100];
	sprintf(buffer, format, v1, v2, v3);
	print(buffer);
}

void print(char *format, float v1, float v2, float v3, float v4)
{
	char buffer[100];
	sprintf(buffer, format, v1, v2, v3, v4);
	print(buffer);
}

// Print an array to the Python console
void print(int *a, int n)
{
	print("(");
	print("%d", a[0]);
	for (int i = 1; i < n; i++)
		print(", %d", a[i]);
	print(")\n");
}
/*
template <class T>
T min(T a, T b) { return a < b ? a : b; }

template <class T>
T max(T a, T b) { return a < b ? b : a; }

template <class T>
T max(T *a, T n)
{
	T current_max = a[0];
	for (T i = 1; i < n; i++)
		current_max = max(a[i], current_max);

	return current_max;
}*/

double abs(double x) { return x < 0 ? -x : x; }

int min(int a, int b) { return a < b ? a : b; }

int max(int a, int b) { return a < b ? b : a; }

int max(int *a, int n)
{
	int current_max = a[0];
	for (int i = 1; i < n; i++)
		current_max = max(a[i], current_max);

	return current_max;
}




/// Take a numpy array and convert it to a c array of integers.
/// Currently only converts 1D arrays.
int* CopyPyArrayObjectToIntArray(PyArrayObject *a)
{
	int n = a->dimensions[0];
	int *copy = new int[n];

	for (int i = 0; i < n; i++)
		copy[i] = GetInt(a, i);

	return copy;
}

/// Take a numpy array and convert it to a c array of doubles.
/// Currently only converts 1D arrays.
int boob = 0;
double* CopyPyArrayObjectToDoubleArray(PyArrayObject *a)
{
	int n = a->dimensions[0];
	double *copy = new double[n];

	for (int i = 0; i < n; i++)
		copy[i] = GetDouble(a, i);

	//delete[] copy;return copy;
	return copy;
}

/// Check that a PyArrayObject has a certain type and dimensions.
bool CheckArray (PyArrayObject* a, int type, int num_dims, int* dims, string name)
{
	// Check number of dimensions
	if (a->nd != num_dims) { PyErr_SetString (PyExc_ValueError, (name + " has wrong number of dimensions").c_str()); return false; }
	
	// Check type
	if (a->descr->type_num != type) { PyErr_SetString (PyExc_ValueError, (name + " has wrong type").c_str()); return false; }

	// Check dimensions
	for (int i = 0; i < num_dims; i++)
		if (a->dimensions[i] != dims[i] && dims[i] != 0) {
			PyErr_SetString (PyExc_ValueError, (name + " has wrong size").c_str());
			return false;
		}

	return true;
}

/// Check that a PyArrayObject has a certain type and numer of dimensions, and check the size of the first dimension.
bool CheckArray (PyArrayObject* a, int type, int num_dims, int len, string name)
{
	// Check number of dimensions
	if (a->nd != num_dims) { PyErr_SetString (PyExc_ValueError, (name + " has wrong number of dimensions").c_str()); return false; }
	
	// Check type
	if (a->descr->type_num != type) { PyErr_SetString (PyExc_ValueError, (name + " has wrong type").c_str()); return false; }
	
	// Check the size of the first dimension.
	if (a->dimensions[0] != len) {PyErr_SetString (PyExc_ValueError, (name + " has wrong size").c_str()); return false; }

	return true;
}
