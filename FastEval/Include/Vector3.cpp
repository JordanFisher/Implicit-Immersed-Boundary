#include "Vector3.h"

Vector3::Vector3() { x = y = z = 0; }
Vector3::Vector3(double _x, double _y, double _z) { x = _x; y = _y; z = _z; }
void Vector3::Abs() { x = abs(x); y = abs(y); z = abs(z); }

Tensor33::Tensor33() { xx = xy = xz = yx = yy = yz = zx = zy = zz = 0; }
void Tensor33::operator=(double b)
{
	xx = b;
	xy = b;
	xz = b;
	yx = b;
	yy = b;
	yz = b;
	zx = b;
	zy = b;
	zz = b;
}

void Tensor33::MultAdd(Tensor33 &A, Tensor33 &B)
{
	xx += A.xx * B.xx;
	xy += A.xy * B.xy;
	xz += A.xz * B.xz;
	yx += A.yx * B.yx;
	yy += A.yy * B.yy;
	yz += A.yz * B.yz;
	zx += A.zx * B.zx;
	zy += A.zy * B.zy;
	zz += A.zz * B.zz;
}

void Tensor33::MultAdd(Tensor33 &A, Vector3 v)
{
	xx += A.xx * v.x;
	xy += A.xy * v.x;
	xz += A.xz * v.x;
	yx += A.yx * v.y;
	yy += A.yy * v.y;
	yz += A.yz * v.y;
	zx += A.zx * v.z;
	zy += A.zy * v.z;
	zz += A.zz * v.z;
}

double* Tensor33::get(int i)
{
	switch(i) {
		case 0: return &xx;
		case 1: return &xy;
		case 2: return &xz;
		case 3: return &yx;
		case 4: return &yy;
		case 5: return &yz;
		case 6: return &zx;
		case 7: return &zy;
		case 8: return &zz;
	}

	return 0;
}


void operator+=(Tensor33 &a, Tensor33 &b)
{
	a.xx += b.xx;
	a.xy += b.xy;
	a.xz += b.xz;
	a.yx += b.yx;
	a.yy += b.yy;
	a.yz += b.yz;
	a.zx += b.zx;
	a.zy += b.zy;
	a.zz += b.zz;
}

void operator*=(Tensor33 &a, Tensor33 &b)
{
	a.xx *= b.xx;
	a.xy *= b.xy;
	a.xz *= b.xz;
	a.yx *= b.yx;
	a.yy *= b.yy;
	a.yz *= b.yz;
	a.zx *= b.zx;
	a.zy *= b.zy;
	a.zz *= b.zz;
}



bool operator==(Vector3 a, Vector3 b)
{
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

Vector3 operator+(Vector3 a, Vector3 b)
{
	return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector3 operator-(Vector3 a, Vector3 b)
{
	return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector3 operator*(double a, Vector3 b)
{
	return Vector3(a * b.x, a * b.y, a * b.z);
}

Vector3 operator*(Vector3 a, double b)
{
	return Vector3(a.x * b, a.y * b, a.z * b);
}

Vector3 operator*(Vector3 a, int b)
{
	return Vector3(a.x * b, a.y * b, a.z * b);
}

Vector3 operator/(Vector3 a, double b)
{
	return Vector3(a.x / b, a.y / b, a.z / b);
}

 void operator+=(Vector3 &a, Vector3 b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

void operator-=(Vector3 &a, Vector3 b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

void operator*=(Tensor33 &a, Vector3 &b)
{
	a.xx *= b.x;
	a.xy *= b.x;
	a.xz *= b.x;
	a.yx *= b.y;
	a.yy *= b.y;
	a.yz *= b.y;
	a.zx *= b.z;
	a.zy *= b.z;
	a.zz *= b.z;
}

Vector3 GetVector3(PyArrayObject* list, int i) {
	return Vector3(GetDouble(list,i,0), GetDouble(list,i,1), GetDouble(list,i,2)); };

void PrintVector3 (Vector3 &vec)
{
	fprintf(stderr, "(%f %f %f)\n", vec.x, vec.y, vec.z);
}

Vector3 PyArrayToVector3 (PyArrayObject *a)
{
	Vector3 vec;
	vec.x = GetDouble(a, 0);
	vec.y = GetDouble(a, 1);
	vec.z = GetDouble(a, 2);

	return vec;
}

PyObject* Vector3ToTuple (Vector3 &vec)
{
	return Py_BuildValue("(ddd)", vec.x, vec.y, vec.z);
}