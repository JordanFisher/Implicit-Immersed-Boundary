#ifndef VECTOR3_H
#define VECTOR3_H

#include<Python.h>
#include "arrayobject.h"

#include "cPythonUtility.h"

struct Vector3 {
	double x, y, z;
	
	Vector3();
	Vector3(double _x, double _y, double _z);
	void Abs();
};

struct Tensor33 {
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;

	Tensor33();
	void operator=(double b);
	void MultAdd(Tensor33 &A, Tensor33 &B);
	void MultAdd(Tensor33 &A, Vector3 v);
	double* get(int i);
};

void operator+=(Tensor33 &a, Tensor33 &b);
void operator*=(Tensor33 &a, Tensor33 &b);
bool operator==(Vector3 a, Vector3 b);
Vector3 operator+(Vector3 a, Vector3 b);
Vector3 operator-(Vector3 a, Vector3 b);
Vector3 operator*(double a, Vector3 b);
Vector3 operator*(Vector3 a, double b);
Vector3 operator*(Vector3 a, int b);
Vector3 operator/(Vector3 a, double b);
void operator+=(Vector3 &a, Vector3 b);
void operator-=(Vector3 &a, Vector3 b);
void operator*=(Tensor33 &a, Vector3 &b);

Vector3 GetVector3(PyArrayObject* list, int i);

void PrintVector3 (Vector3 &vec);

Vector3 PyArrayToVector3 (PyArrayObject *a);

PyObject* Vector3ToTuple (Vector3 &vec);

#endif