#ifndef CDERIVE_H
#define CDERIVE_H

#include<Python.h>
#include"arrayobject.h"

#include"cPythonUtility.h"

void eno2(int n, double h, PyArrayObject* u, PyArrayObject* phi, PyArrayObject* f);
void eno3(int n, double h, PyArrayObject* u, PyArrayObject* phi, PyArrayObject* f);
void _eno2(int N, double h, double* u, double* phi, PyArrayObject* f);
void _eno3(int N, double h, double* u, double* phi, PyArrayObject* f);

/// Parses and error checks arguments before sending on to eno2
static PyObject* cEno2 (PyObject *self, PyObject *args)
{
	PyArrayObject *u, *phi, *f;
	double h;
	int n;

	if (!PyArg_ParseTuple(args, "idOOO", &n, &h, &u, &phi, &f)) {
		PyErr_SetString(PyExc_ValueError, "C++ prototype doesn't match Python input");

		return NULL;
	}

	if (!CheckArray(u, PyArray_DOUBLE, 1, n, "X")) return NULL;
	if (!CheckArray(phi, PyArray_DOUBLE, 1, n, "phi")) return NULL;
	if (!CheckArray(f, PyArray_DOUBLE, 1, n, "f")) return NULL;
	
	eno2(n, h, u, phi, f);

	Py_INCREF(Py_None);
	return Py_None;
}

/// Parses and error checks arguments before sending on to eno3
static PyObject* cEno3 (PyObject *self, PyObject *args)
{
	PyArrayObject *u, *phi, *f;
	double h;
	int n;

	if (!PyArg_ParseTuple(args, "idOOO", &n, &h, &u, &phi, &f)) {
		PyErr_SetString(PyExc_ValueError, "C++ prototype doesn't match Python input");

		return NULL;
	}

	if (!CheckArray(u, PyArray_DOUBLE, 1, n, "X")) return NULL;
	if (!CheckArray(phi, PyArray_DOUBLE, 1, n, "phi")) return NULL;
	if (!CheckArray(f, PyArray_DOUBLE, 1, n, "f")) return NULL;
	
	eno3(n, h, u, phi, f);

	Py_INCREF(Py_None);
	return Py_None;
}

void Convect_Centered(int n, double h, PyArrayObject* u, PyArrayObject* phi, PyArrayObject* f);
void _Convect_Centered(int N, double h, double* u, double* phi, PyArrayObject* f);

/// Parses and error checks arguments before sending on to eno2
static PyObject* cConvect_Centered (PyObject *self, PyObject *args)
{
	PyArrayObject *u, *phi, *f;
	double h;
	int n;

	if (!PyArg_ParseTuple(args, "idOOO", &n, &h, &u, &phi, &f)) {
		PyErr_SetString(PyExc_ValueError, "C++ prototype doesn't match Python input");

		return NULL;
	}

	if (!CheckArray(u, PyArray_DOUBLE, 1, n, "X")) return NULL;
	if (!CheckArray(phi, PyArray_DOUBLE, 1, n, "phi")) return NULL;
	if (!CheckArray(f, PyArray_DOUBLE, 1, n, "f")) return NULL;
	
	Convect_Centered(n, h, u, phi, f);

	Py_INCREF(Py_None);
	return Py_None;
}

#endif