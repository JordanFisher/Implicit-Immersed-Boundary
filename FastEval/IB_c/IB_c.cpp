#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include<Python.h>
#include"arrayobject.h"

#include"cPythonUtility.h"
#include"cDerive.h"

const double pi = 3.14159265358979323;
const int DeltaHalfWidth[] = {1, 1, 0};
const int DeltaWidth[] = {3, 3, 1};
double invh;

/// Calculate Delta(r) given h and the type of Delta function to use. invh must equal 1 / h.
double Delta (double h, double invh, double r, int DeltaType)
{
	static double x, absx;

	switch (DeltaType)
	{
		case 0:
			if (abs(r) < 2 * h)
				return (1. + cos(pi * r / (2*h))) / (4. * h);
			else
				return 0;

		case 1:
			x = r / h;
			absx = abs(x);
			if (absx <= 2)
				if (absx <= 1)
					return .125 * (3. - 2 * absx + sqrt(1. + 4. * absx - 4. * x * x)) / h;
				else
					return .125 * (5. - 2 * absx - sqrt(-7. + 12. * absx - 4 * x * x)) / h;
			else
				return 0;

		case 2:
			absx = abs(r * invh);
			if (absx < 1)
				return (1. - absx) * invh;
			else
				return 0;
	}

	PyErr_SetString (PyExc_ValueError, "Unknown delta type");

	return NULL;
}

/// Perform a linear interpolation between v0 and v1. The variable a should be a value between 0 and 1.
double LinearInterpolation (double a, double v0, double v1)
{
	return (1 - a) * v0 + a * v1;
}

/// Perform a bilinear interpolation between v00, v01, v10, v11. The variables a and b should be between 0 and 1.
double BiLinearInterpolation (double a, double b, double v00, double v01, double v10, double v11)
{
	return (1 - a) * LinearInterpolation (b, v00, v01) + a * LinearInterpolation (b, v10, v11);
}

/// Perform a trilinear interpolation between v000, v001, v010, v011, v100, v101, v110, v111. The variables a, b, and c should be between 0 and 1.
double TriLinearInterpolation (double a, double b, double c, double v000, double v001, double v010, double v011, double v100, double v101, double v110, double v111)
{
	return (1 - a) * BiLinearInterpolation (b, c, v000, v001, v010, v011) + a * BiLinearInterpolation (b, c, v100, v101, v110, v111);
}

int mod (int a, int m) {
	return (a % m + m) % m;
}

/// Construct the fluid matrix that encodes the interaction between all fiber points
/// X is an Nb x 2 or Nb x 3 array storing the fiber point positions
/// G is the discrete Greens function of the Spread/Interpolated Stokeslet
/// A is 2Nb x 2Nb or 3Nb x 3Nb array that will store the fluid matrix
/// band, when smaller than Nb, forces A to be zero for entries more than band away from the diagonal
/// Dim is the number of spatial dimensions of the underlying fluid domain
void ComputeMatrix(int Dim, PyArrayObject *X, int* N, int Nb, double* h, double hb, PyArrayObject *G, PyArrayObject *A, int band, bool MatrixOfTensors)
{
	// Note: The comments for this function assume Dim == 3. The code works when Dim == 2 as well.

	int j, k, _j, _k, _l;
	int i1, i2;
	double a, b, c, a1b1c1, a1b1c0, a1b0c1, a1b0c0, a0b1c1, a0b1c0, a0b0c1, a0b0c0;

	// If band width is negative set to a very wide width
	if (band < 0) band = max(N, Dim);

	// Loop through all the fiber points, X_j
	for (j = 0; j < Nb; j++) {
		// Loop through all the fiber points, X_k
		for (k = j; k < Nb; k++) {
			// Calculate the distance in grid cells (a,b,c) between X_j and X_k
			// Continue to next iteration if the distance is more than band
			a = (GetDouble(X,j,0) - GetDouble(X,k,0)) / h[0];
			if (abs(a) > band) continue;
			b = (GetDouble(X,j,1) - GetDouble(X,k,1)) / h[1];
			if (abs(b) > band) continue;
			if (Dim == 3) {
				c = (GetDouble(X,j,2) - GetDouble(X,k,2)) / h[2];
				if (abs(c) > band) continue; }

			// Find the Eulerian intersection corresponding to (a,b,c), rounding down where needed
			_j = (int) floor(a);
			_k = (int) floor(b);
			if (Dim == 3)
				_l = (int) floor(c);

			// Store in (a,b,c) the fractional position within the Eulerian cell (_j,_k,_l)
			a -= (double) _j;
			b -= (double) _k;
			if (Dim == 3)
				c -= (double) _l;

			// Add N to (_j,_k,_l) so that modulo arithmetic works correctly
			_j = _j + N[0];
			_k = _k + N[1];
			if (Dim == 3)
				_l = _l + N[2];

			// We will look up the value G(_j + a, _k + b, _l + c), using trilinear interpolation between grid cells when (a,b,c) != 0
			// First calculate the scalars needed for trilinear interpolation. These can be reused for the 9 components of the tensor relating X_j and X_k
			if (Dim == 2)
			{
				a1b1c0 = (1. - a) * (1. - b) * hb;
				a1b0c0 = (1. - a) * b * hb;
				a0b1c0 = a * (1. - b) * hb;
				a0b0c0 = a * b * hb;
			}
			else
			{
				a1b1c1 = (1. - a) * (1. - b) * (1. - c) * hb;
				a1b1c0 = (1. - a) * (1. - b) * c * hb;
				a1b0c1 = (1. - a) * b * (1. - c) * hb;
				a1b0c0 = (1. - a) * b * c * hb;
				a0b1c1 = a * (1. - b) * (1. - c) * hb;
				a0b1c0 = a * (1. - b) * c * hb;
				a0b0c1 = a * b * (1. - c) * hb;
				a0b0c0 = a * b * c * hb;
			}

			// Lookup the 9 components of the tensor relating X_j and X_k
			for (i1 = 0; i1 < Dim; i1++)
			{
				for (i2 = 0; i2 < Dim; i2++)
				{
					// Lookup the (i1,i2)-th component of the tensor
					double val;
					if (Dim == 2) val =
						// Bilinear interpolation
						a1b1c0 * GetDouble(G, _j     % N[0],  _k    % N[1], i2, i1) + 
						a1b0c0 * GetDouble(G, _j     % N[0], (_k+1) % N[1], i2, i1) + 
						a0b1c0 * GetDouble(G, (_j+1) % N[0],  _k    % N[1], i2, i1) + 
						a0b0c0 * GetDouble(G, (_j+1) % N[0], (_k+1) % N[1], i2, i1);
					else val =
						// Trilinear interpolation
						a1b1c1 * GetDouble(G,  _j    % N[0],  _k    % N[1],  _l    % N[2], i2, i1) + 
						a1b1c0 * GetDouble(G,  _j    % N[0],  _k    % N[1], (_l+1) % N[2], i2, i1) + 
						a1b0c1 * GetDouble(G,  _j    % N[0], (_k+1) % N[1],  _l    % N[2], i2, i1) + 
						a1b0c0 * GetDouble(G,  _j    % N[0], (_k+1) % N[1], (_l+1) % N[2], i2, i1) + 
						a0b1c1 * GetDouble(G, (_j+1) % N[0],  _k    % N[1],  _l    % N[2], i2, i1) + 
						a0b1c0 * GetDouble(G, (_j+1) % N[0],  _k    % N[1], (_l+1) % N[2], i2, i1) + 
						a0b0c1 * GetDouble(G, (_j+1) % N[0], (_k+1) % N[1],  _l    % N[2], i2, i1) + 
						a0b0c0 * GetDouble(G, (_j+1) % N[0], (_k+1) % N[1], (_l+1) % N[2], i2, i1);

					// Store the tensor in the matrix A
					if (MatrixOfTensors)
						GetDouble(A, j, k, i1, i2) = val;
					else
						GetDouble(A, j + i1 * Nb, k + i2 * Nb) = val;

					// The matrix A is block symmetric, copy values as needed
					if (j != k)
					{
						if (MatrixOfTensors)
							GetDouble(A, k, j, i1, i2) = val;
						else
							GetDouble(A, k + i1 * Nb, j + i2 * Nb) = val;
					}
				}
			}
		}
	}
}

/// Parses and error checks arguments before sending on to ComputeMatrix
static PyObject* cComputeMatrix (PyObject *self, PyObject *args)
{
	PyArrayObject *X, *G, *A, *N, *h;
	int Nb;
	double hb;
	int band;
	
	bool MatrixOfTensors = false; // When true the matrix A is stored as an Nb x Nb array of 2 x 2 or 3 x 3 tensors

	if (!PyArg_ParseTuple(args, "OOiOdOOi", &X, &N, &Nb, &h, &hb, &G, &A, &band)) {
		PyErr_SetString(PyExc_ValueError, "C++ prototype doesn't match Python input");

		return NULL;
	}

	// Get the number of spatial dimensions of the underlying fluid domain
	int Dim = N->dimensions[0];

	// Check the arrays N, h
	if (!CheckArray(N, PyArray_LONG, 1, Dim, "N")) return NULL;
	if (!CheckArray(h, PyArray_DOUBLE, 1, Dim, "h")) return NULL;
	
	int *_N = CopyPyArrayObjectToIntArray(N);
	double *_h = CopyPyArrayObjectToDoubleArray(h);
	
	// Calculate the shape of G
	int *Gdims = new int[Dim + 2];
	for (int i = 0; i < Dim; i++)
		Gdims[i] = _N[i];
	Gdims[Dim] = Gdims[Dim+1] = Dim;

	// Check the arrays X, U, A
	int Xdims[] = {Nb, Dim};
	if (!CheckArray(X, PyArray_DOUBLE, 2, Xdims, "X")) return NULL;
	if (!CheckArray(G, PyArray_DOUBLE, Dim + 2, Gdims, "G")) return NULL;
	
	if (MatrixOfTensors)
	{
		int Adims[] = {Nb, Nb, Dim, Dim};
		if (!CheckArray(A, PyArray_DOUBLE, 4, Adims, "A")) return NULL;
	}
	else
	{
		int Adims[] = {Dim * Nb, Dim * Nb};
		if (!CheckArray(A, PyArray_DOUBLE, 2, Adims, "A")) return NULL;
	}

	ComputeMatrix(Dim, X, _N, Nb, _h, hb, G, A, band, MatrixOfTensors);

	delete[] _N;
  delete[] _h;

	Py_INCREF(Py_None);
	return Py_None;
}


void WholeGridSpread_2D(PyArrayObject *u, int N[2], double h[2], PyArrayObject *U, int range, int DeltaType)
{
	int j, k, _j, _k, r;
	double _delt, delt, x, y;

	double invh[] = { double(N[0]), double(N[1]) };
	double h2 = h[0] * h[1];

	// Initialize the interpolated velocity field U to zero.
	for (j = 0; j < N[0]; j++) {
		for (k = 0; k < N[1]; k++) {
			GetDouble(U,j,k,0) = GetDouble(U,j,k,1) = 0;
		}
	}

	// Get the Delta width and half width
	int hw = DeltaHalfWidth[DeltaType];
	int w = DeltaWidth[DeltaType];

	if (range > N[0]) range = N[0];

	// Loop through all the Eulerian points within a distance of range from the origin
	for (_j = 0; _j < range; _j++) {
		for (_k = 0; _k < range; _k++) {
			x = _j * h[0];
			y = _k * h[1];
			
			for (j = _j - hw; j <= _j + hw; j++) {
				_delt = Delta(h[0], invh[0], j * h[0] - x, DeltaType) * h2;
				for (k = _k - hw; k <= _k + hw; k++) {
					delt = _delt * Delta(h[1], invh[0], k * h[1] - y, DeltaType);
					GetDouble(U,_j,_k,0) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), 0) * delt;
					GetDouble(U,_j,_k,1) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), 1) * delt;
				}
			}
		}
	}
}

void WholeGridSpread_3D(PyArrayObject *u, int N[3], double h[3], PyArrayObject *U, int range, int DeltaType)
{
	int j, k, l, _j, _k, _l, r;
	double __delt, _delt, delt, x, y, z;

	double invh[] = { double(N[0]), double(N[1]) };
	double h3 = h[0] * h[1] * h[2];

	// Initialize the interpolated velocity field U to zero.
	for (j = 0; j < N[0]; j++) {
		for (k = 0; k < N[1]; k++) {
			for (l = 0; l < N[2]; l++) {
				GetDouble(U,j,k,l,0) = GetDouble(U,j,k,l,1) = GetDouble(U,j,k,l,2) = 0;
			}
		}
	}

	// Get the Delta width and half width
	int hw = DeltaHalfWidth[DeltaType];
	int w = DeltaWidth[DeltaType];

	if (range > N[0]) range = N[0];

	// Loop through all the Eulerian points within a distance of range from the origin
	for (_j = 0; _j < range; _j++) {
		for (_k = 0; _k < range; _k++) {
			for (_l = 0; _l < range; _l++) {
				x = _j * h[0];
				y = _k * h[1];
				z = _l * h[2];

				for (j = _j - hw; j <= _j + hw; j++) {
					__delt = Delta(h[0], invh[0], j * h[0] - x, DeltaType) * h3;
					for (k = _k - hw; k <= _k + hw; k++) {
						_delt = __delt * Delta(h[1], invh[0], k * h[1] - y, DeltaType);
						for (l = _l - hw; l <= _l + hw; l++) {
							delt = _delt * Delta(h[2], invh[0], l * h[2] - z, DeltaType);
							GetDouble(U,_j,_k,_l,0) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 0) * delt;
							GetDouble(U,_j,_k,_l,1) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 1) * delt;
							GetDouble(U,_j,_k,_l,2) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 2) * delt;
						}
					}
				}
			}
		}
	}
}

static PyObject* cWholeGridSpread(PyObject *self, PyObject *args)
{
	PyArrayObject *U, *u, *N, *h;
	int DeltaType;
	int range;

	if (!PyArg_ParseTuple(args, "OOOOii", &u, &N, &h, &U, &range, &DeltaType)) {
		PyErr_SetString(PyExc_ValueError, "C++ prototype doesn't match Python input");

		return NULL;
	}

	// Get the number of spatial dimensions of the underlying fluid domain
	int Dim = N->dimensions[0];

	// Check the arrays N, h
	if (!CheckArray(N, PyArray_LONG, 1, Dim, "N")) return NULL;
	if (!CheckArray(h, PyArray_DOUBLE, 1, Dim, "h")) return NULL;
	
	int *_N = CopyPyArrayObjectToIntArray(N);
	double *_h = CopyPyArrayObjectToDoubleArray(h);
	
	// Calculate the shape of u and U
	int *ushape = new int[Dim + 1];
	for (int i = 0; i < Dim; i++)
		ushape[i] = _N[i];
	ushape[Dim] = Dim;

	// Check the arrays u, U
	if (!CheckArray(u, PyArray_DOUBLE, Dim + 1, ushape, "u")) return NULL;
	if (!CheckArray(U, PyArray_DOUBLE, Dim + 1, ushape, "U")) return NULL;

	// Perform the spreading
	if (Dim == 2)
		WholeGridSpread_2D(u, _N, _h, U, range, DeltaType);
	else if (Dim == 3)
		WholeGridSpread_3D(u, _N, _h, U, range, DeltaType);
	else
	{
		PyErr_SetString(PyExc_ValueError, "N must be a numpy array of length 2 or 3.");

		return NULL;
	}

	delete[] _N;
  delete[] _h;

	Py_INCREF(Py_None);
	return Py_None;
}


void FiberToGrid_2D(int N[2], double h[2], int Nb, double hb, PyArrayObject *X, PyArrayObject *F, PyArrayObject *f, int DeltaType)
{
	double invh[] = { double(N[0]), double(N[1]) };

	int i, j, k, jMin, jMax, kMin, kMax;
	double _delt, delt;

	// Initialize the force field f to zero.
	for (j = 0; j < N[0]; j++) {
		for (k = 0; k < N[1]; k++) {
				GetDouble(f,j,k,0) = 0.;
		}
	}

	// Get the Delta width and half width
	int hw = DeltaHalfWidth[DeltaType];
	int w = DeltaWidth[DeltaType];

	// Loop through each fiber point
	for (i = 0; i < Nb; i++)
	{
		// Form a neighborhood around a given fiber point
		jMin = floor(GetDouble(X,i,0) / h[0] - hw);
		jMax = jMin + w;
		kMin = floor(GetDouble(X,i,1) / h[1] - hw);
		kMax = kMin + w;	

		// Loop through all Eulerian points in the neighborhood
		for (j = jMin; j <= jMax; j++) {
			_delt = Delta(h[0], invh[0], j * h[0] - GetDouble(X,i,0), DeltaType) * hb;
			for (k = kMin; k <= kMax; k++) {
				// Calculate the 2D delta function of the difference between the Eulerian point and the fiber point
				delt = _delt * Delta(h[1], invh[1], k * h[1] - GetDouble(X,i,1), DeltaType);

				// Add the delta function value times the fiber force to the Eulerian fiber force field
				GetDouble(f, mod(j,N[0]), mod(k,N[1]), 0) += GetDouble(F,i,0) * delt;
				GetDouble(f, mod(j,N[0]), mod(k,N[1]), 1) += GetDouble(F,i,1) * delt;
			}
		}
	}
}

void FiberToGrid_3D(int N[3], double h[3], int Nb, double hb, PyArrayObject *X, PyArrayObject *F, PyArrayObject *f, int DeltaType)
{
	double invh[] = { double(N[0]), double(N[1]), double(N[2]) };

	int i, j, k, l, jMin, jMax, kMin, kMax, lMin, lMax;
	double __delt, _delt, delt;

	// Initialize the force field f to zero.
	for (j = 0; j < N[0]; j++) {
		for (k = 0; k < N[1]; k++) {
			for (l = 0; l < N[2]; l++) {
				GetDouble(f,j,k,l,0) = 0.;
			}
		}
	}

	// Get the Delta width and half width
	int hw = DeltaHalfWidth[DeltaType];
	int w = DeltaWidth[DeltaType];

	// Loop through each fiber point
	for (i = 0; i < Nb; i++)
	{
		// Form a neighborhood around a given fiber point
		jMin = floor(GetDouble(X,i,0) / h[0] - hw);
		jMax = jMin + w;
		kMin = floor(GetDouble(X,i,1) / h[1] - hw);
		kMax = kMin + w;
		lMin = floor(GetDouble(X,i,2) / h[2] - hw);
		lMax = lMin + w;	

		// Loop through all Eulerian points in the neighborhood
		for (j = jMin; j <= jMax; j++) {
			__delt = Delta(h[0], invh[0], j * h[0] - GetDouble(X,i,0), DeltaType) * hb;
			for (k = kMin; k <= kMax; k++) {
				_delt = __delt * Delta(h[1], invh[1], k * h[1] - GetDouble(X,i,1), DeltaType);
				for (l = lMin; l <= lMax; l++) {
					// Calculate the 3D delta function of the difference between the Eulerian point and the fiber point
					delt = _delt * Delta(h[2], invh[2], l * h[2] - GetDouble(X,i,2), DeltaType);

					// Add the delta function value times the fiber force to the Eulerian fiber force field
					GetDouble(f, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 0) += GetDouble(F,i,0) * delt;
					GetDouble(f, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 1) += GetDouble(F,i,1) * delt;
					GetDouble(f, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 2) += GetDouble(F,i,2) * delt;
				}
			}
		}
	}
}

static PyObject* cFiberToGrid (PyObject *self, PyObject *args)
{
	PyArrayObject *X, *F, *f, *N, *h;
	int Nb;
	double hb;
	int DeltaType;

	if (!PyArg_ParseTuple(args, "OOidOOOi", &N, &h, &Nb, &hb, &X, &F, &f, &DeltaType)) {
		PyErr_SetString(PyExc_ValueError, "C++ prototype doesn't match Python input");

		return NULL;
	}

	// Get the number of spatial dimensions of the underlying fluid domain
	int Dim = N->dimensions[0];

	// Check the arrays N, h
	if (!CheckArray(N, PyArray_LONG, 1, Dim, "N")) return NULL;
	if (!CheckArray(h, PyArray_DOUBLE, 1, Dim, "h")) return NULL;
	
	int *_N = CopyPyArrayObjectToIntArray(N);
	double *_h = CopyPyArrayObjectToDoubleArray(h);
	
	// Calculate the shape of f
	int *fshape = new int[Dim + 1];
	for (int i = 0; i < Dim; i++)
		fshape[i] = _N[i];
	fshape[Dim] = Dim;

	// Check the arrays X, F, f
	int dims[] = {Nb, Dim};
	if (!CheckArray(X, PyArray_DOUBLE, 2, dims, "X")) return NULL;
	if (!CheckArray(F, PyArray_DOUBLE, 2, dims, "F")) return NULL;
	if (!CheckArray(f, PyArray_DOUBLE, Dim + 1, fshape, "f")) return NULL;

	// Perform the spreading
	if (Dim == 2)
		FiberToGrid_2D(_N, _h, Nb, hb, X, F, f, DeltaType);
	else if (Dim == 3)
		FiberToGrid_3D(_N, _h, Nb, hb, X, F, f, DeltaType);
	else
	{
		PyErr_SetString(PyExc_ValueError, "N must be a numpy array of length 2 or 3.");

		return NULL;
	}

	delete[] _N;
  delete[] _h;

	Py_INCREF(Py_None);
	return Py_None;
}

void GridToFiber_2D(int N[2], double h[2], int Nb, double hb, PyArrayObject *X, PyArrayObject *u, PyArrayObject *U, int DeltaType)
{
	int i, j, k, jMin, jMax, kMin, kMax;
	double _delt, delt;

	double h2 = h[0] * h[1];
	double invh[] = { double(N[0]), double(N[1]) };

	// Initialize the fiber velocity U to zero.
	for (i = 0; i < Nb; i++)
		GetDouble(U,i,0) = GetDouble(U,i,1) = 0;

	// Get the Delta width and half width
	int hw = DeltaHalfWidth[DeltaType];
	int w = DeltaWidth[DeltaType];

	// Loop through each fiber point
	for (i = 0; i < Nb; i++) {
		jMin = floor(GetDouble(X,i,0) / h[0] - hw);
		jMax = jMin + w;
		kMin = floor(GetDouble(X,i,1) / h[1] - hw);
		kMax = kMin + w;
	
		// Loop through all Eulerian points in the neighborhood
		for (j = jMin; j <= jMax; j++) {
			_delt = Delta(h[0], invh[0], j * h[0] - GetDouble(X,i,0), DeltaType) * h2;
			for (k = kMin; k <= kMax; k++) {
				// Calculate the 2D delta function of the difference between the Eulerian point and the fiber point
				delt = _delt * Delta(h[1], invh[1], k * h[1] - GetDouble(X,i,1), DeltaType);

				// Add the delta function value times the Eulerian velocity to the fiber velocity
				GetDouble(U,i,0) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), 0) * delt;
				GetDouble(U,i,1) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), 1) * delt;
			}
		}
	}
}

void GridToFiber_3D(int N[3], double h[3], int Nb, double hb, PyArrayObject *X, PyArrayObject *u, PyArrayObject *U, int DeltaType)
{
	int i, j, k, l, jMin, jMax, kMin, kMax, lMin, lMax;
	double __delt, _delt, delt;

	double h3 = h[0] * h[1] * h[2];
	double invh[] = { double(N[0]), double(N[1]), double(N[2]) };

	// Initialize the fiber velocity U to zero.
	for (i = 0; i < Nb; i++)
		GetDouble(U,i,0) = GetDouble(U,i,1) = GetDouble(U,i,2) = 0;

	// Get the Delta width and half width
	int hw = DeltaHalfWidth[DeltaType];
	int w = DeltaWidth[DeltaType];

	// Loop through each fiber point
	for (i = 0; i < Nb; i++) {
		jMin = floor(GetDouble(X,i,0) / h[0] - hw);
		jMax = jMin + w;
		kMin = floor(GetDouble(X,i,1) / h[1] - hw);
		kMax = kMin + w;
		lMin = floor(GetDouble(X,i,2) / h[2] - hw);
		lMax = lMin + w;
	
		// Loop through all Eulerian points in the neighborhood
		for (j = jMin; j <= jMax; j++) {
			__delt = Delta(h[0], invh[0], j * h[0] - GetDouble(X,i,0), DeltaType) * h3;
			for (k = kMin; k <= kMax; k++) {
				_delt = __delt * Delta(h[1], invh[1], k * h[1] - GetDouble(X,i,1), DeltaType);
				for (l = lMin; l <= lMax; l++) {
					// Calculate the 3D delta function of the difference between the Eulerian point and the fiber point
					delt = _delt * Delta(h[2], invh[2], l * h[2] - GetDouble(X,i,2), DeltaType);

					// Add the delta function value times the Eulerian velocity to the fiber velocity
					GetDouble(U,i,0) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 0) * delt;
					GetDouble(U,i,1) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 1) * delt;
					GetDouble(U,i,2) += GetDouble(u, mod(j,N[0]), mod(k,N[1]), mod(l,N[2]), 2) * delt;
				}
			}
		}
	}
}

static PyObject* cGridToFiber(PyObject *self, PyObject *args)
{
	PyArrayObject *X, *u, *U, *N, *h;
	int Nb;
	double hb;
	int DeltaType;

	if (!PyArg_ParseTuple(args, "OOidOOOi", &N, &h, &Nb, &hb, &X, &u, &U, &DeltaType)) {
		PyErr_SetString(PyExc_ValueError, "C++ prototype doesn't match Python input");

		return NULL;
	}

	// Get the number of spatial dimensions of the underlying fluid domain
	int Dim = N->dimensions[0];

	// Check the arrays N, h
	if (!CheckArray(N, PyArray_LONG, 1, Dim, "N")) return NULL;
	if (!CheckArray(h, PyArray_DOUBLE, 1, Dim, "h")) return NULL;
	
	int *_N = CopyPyArrayObjectToIntArray(N);
	double *_h = CopyPyArrayObjectToDoubleArray(h);

	// Calculate the shape of u
	int *ushape = new int[Dim + 1];
	for (int i = 0; i < Dim; i++)
		ushape[i] = _N[i];
	ushape[Dim] = Dim;

	// Check the arrays X, U, u
	int dims[] = {Nb, Dim};
	if (!CheckArray(X, PyArray_DOUBLE, 2, dims, "X")) return NULL;
	if (!CheckArray(U, PyArray_DOUBLE, 2, dims, "U")) return NULL;
	if (!CheckArray(u, PyArray_DOUBLE, Dim + 1, ushape, "u")) return NULL;

	// Perform the interpolation
	if (Dim == 2)
		GridToFiber_2D(_N, _h, Nb, hb, X, u, U, DeltaType);
	else if (Dim == 3)
		GridToFiber_3D(_N, _h, Nb, hb, X, u, U, DeltaType);
	else
	{
		PyErr_SetString(PyExc_ValueError, "N must be a numpy array of length 2 or 3.");

		return NULL;
	}

	delete[] _N;
  delete[] _h;

	Py_INCREF(Py_None);
	return Py_None;
}

static PyMethodDef IB_c_methods[] = {
	{"ComputeMatrix", cComputeMatrix, METH_VARARGS, "ComputeMatrix (X, N, Nb, h, hb, G, A, band)"},
	{"Eno2", cEno2, METH_VARARGS, "Eno2 (N, h, u, phi, f)"},
	{"Eno3", cEno3, METH_VARARGS, "Eno3 (N, h, u, phi, f)"},
	{"Convect_Centered", cConvect_Centered, METH_VARARGS, "Eno2 (N, h, u, phi, f)"},
	{"FiberToGrid", cFiberToGrid, METH_VARARGS, "FiberToGrid (N, h, Nb, hb, X, F, f, DeltaType)"},
	{"GridToFiber", cGridToFiber, METH_VARARGS, "GridToFiber (N, h, Nb, hb, X, u, U, DeltaType)"},
	{"WholeGridSpread", cWholeGridSpread, METH_VARARGS, "WholeGridSpread (u, N, h, U, range, DeltaType)"},
	{NULL, NULL}
};


PyMODINIT_FUNC
initIB_c(void)
{
	Py_InitModule("IB_c", IB_c_methods);
}
