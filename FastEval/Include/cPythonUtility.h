#ifndef CPYTHONUTILITY_H
#define CPYTHONUTILITY_H

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include<Python.h>
#include"arrayobject.h"

/// Write a string to Python console
void print(char* str);

/// Format a string and write it to the Python console
void print(char *format, float v1);
void print(char *format, float v1, float v2);
void print(char *format, float v1, float v2, float v3);
void print(char *format, float v1, float v2, float v3, float v4);

// Print an array to the Python console
void print(int *a, int n);

double abs(double x);
int min(int a, int b);
int max(int a, int b);
int max(int *a, int n);

// Array accessors.
// Why are these needed?
//		Python arrays may have non-uniform strides, so we can't use normal C array syntax.
// Usage:
//		If u is a PyArrayObject (a Python array) of type double:
//			GetDouble(u,j,k) = 1,		not		u[j,k] = 1
//		If u is a PyArrayObject (a Python array) of type int:
//			GetInt(u,j,k) = 1,			not		u[j,k] = 1

inline int& GetInt(PyArrayObject* list, int j) {
	return *(int *)(list->data + j * list->strides[0]); };
inline int& GetInt(PyArrayObject* list, int j, int k) {
	return *(int *)(list->data + j * list->strides[0] + k * list->strides[1]); };
inline int& GetInt(PyArrayObject* list, int j, int k, int l) {
	return *(int *)(list->data + j * list->strides[0] + k * list->strides[1] + l * list->strides[2]); };
inline int& GetInt(PyArrayObject* list, int j, int k, int l, int m) {
	return *(int *)(list->data + j * list->strides[0] + k * list->strides[1] + l * list->strides[2] + m * list->strides[3]); };
inline double& GetDouble(PyArrayObject* list, int j) {
	return *(double *)(list->data + j * list->strides[0]); };
inline double& GetDouble(PyArrayObject* list, int j, int k) {
	return *(double *)(list->data + j * list->strides[0] + k * list->strides[1]); };
inline double& GetDouble(PyArrayObject* list, int j, int k, int l) {
	return *(double *)(list->data + j * list->strides[0] + k * list->strides[1] + l * list->strides[2]); };
inline double& GetDouble(PyArrayObject* list, int j, int k, int l, int m) {
	return *(double *)(list->data + j * list->strides[0] + k * list->strides[1] + l * list->strides[2] + m * list->strides[3]); };
inline double& GetDouble(PyArrayObject* list, int j, int k, int l, int m, int n) {
	return *(double *)(list->data + j * list->strides[0] + k * list->strides[1] + l * list->strides[2] + m * list->strides[3] + n * list->strides[4]); };
inline double& GetDouble(PyArrayObject* list, int j, int k, int l, int m, int n, int o, int p, int q) {
	return *(double *)(list->data + j * list->strides[0] + k * list->strides[1] + l * list->strides[2] + m * list->strides[3] + n * list->strides[4] + o * list->strides[5] + p * list->strides[6] + q * list->strides[7]); };

inline double& SafeGetDouble(PyArrayObject* list, int i, int j, int k, int N) {
	return GetDouble(list, i >= 0 ? i : N + i, j >= 0 ? j : N + j, k >= 0 ? k : N + k); }


/// Take a numpy array and convert it to a c array of integers.
/// Currently only converts 1D arrays.
int* CopyPyArrayObjectToIntArray(PyArrayObject *a);

/// Take a numpy array and convert it to a c array of doubles.
/// Currently only converts 1D arrays.
double* CopyPyArrayObjectToDoubleArray(PyArrayObject *a);

/// Check that a PyArrayObject has a certain type and dimensions.
bool CheckArray (PyArrayObject* a, int type, int num_dims, int* dims, string name);

/// Check that a PyArrayObject has a certain type and numer of dimensions, and check the size of the first dimension.
bool CheckArray (PyArrayObject* a, int type, int num_dims, int len, string name);

#endif