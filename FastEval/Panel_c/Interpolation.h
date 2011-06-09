#include<Python.h>
#include "arrayobject.h"

#include "cPythonUtility.h"
#include "Vector3.h"


/// Calculate the value of Field at x, snapping x to the nearest Eulerian grid point
double DirectLookup(PyArrayObject *Field, Vector3 &x)
{
	int i, j, k;
	i = x.x > 0 ? (int)(x.x + .5) : (int)(x.x - .5);
	j = x.y > 0 ? (int)(x.y + .5) : (int)(x.y - .5);
	k = x.z > 0 ? (int)(x.z + .5) : (int)(x.z - .5);
	int N = Field->dimensions[0];

	return SafeGetDouble(Field, i, j, k, N);
}

/// Calculate a trilinear interpolation of Field at position x
/// The Field is assumed to be periodic, and the tensor value to lookup is assumed to be symmetric
void Periodic_Lerp(Tensor33 &output, PyArrayObject *Field, Vector3 &x)
{
	int i, j, k;
	i = x.x >= 0 ? (int)(x.x) : (int)(x.x)-1;
	j = x.y >= 0 ? (int)(x.y) : (int)(x.y)-1;
	k = x.z >= 0 ? (int)(x.z) : (int)(x.z)-1;
	int N = Field->dimensions[2];

	double v1 = x.x - i, v2 = x.y - j, v3 = x.z - k;
	double w1 = 1 - v1, w2 = 1 - v2, w3 = 1 - v3;

	// Use periodicity of Field for negative indices
	if (i < 0) i += N;
	if (j < 0) j += N;
	if (k < 0) k += N;

	// Get the strides for Field
	int stride_a = Field->strides[3], stride_b = Field->strides[4], stride_0 = Field->strides[0], stride_1 = Field->strides[1], stride_2 = Field->strides[2];
	
	// Get the address for Field[i,j,k]
	char *index_000, *index_100, *index_010, *index_001, *index_110, *index_011, *index_101, *index_111;
	index_000 = (Field)->data + i * stride_0 + j * stride_1 + k * stride_2;

	// Get the addresses for Field[i+1,j,k], Field[i,j+1,k], Field[i,j,k+1]
	index_100 = index_000 + stride_0;
	index_010 = index_000 + stride_1;
	index_001 = index_000 + stride_2;
	
	// Get the addresses for Field[i,j+1,k+1], Field[i+1,j+1,k], Field[i+1,j,k+1]
	index_011 = index_001 + stride_1;
	index_110 = index_010 + stride_0;
	index_101 = index_001 + stride_0;

	// Get the address for Field[i+1,j+1,k+1]
	index_111 = index_101 + stride_1;

	// Use periodicity of Field for large indices
	if (i == N-1)
	{
		index_100 -= N * stride_0;
		index_110 -= N * stride_0;
		index_101 -= N * stride_0;
		index_111 -= N * stride_0;
	}

	if (j == N-1)
	{
		index_010 -= N * stride_1;
		index_110 -= N * stride_1;
		index_011 -= N * stride_1;
		index_111 -= N * stride_1;
	}

	if (k == N-1)
	{
		index_001 -= N * stride_2;
		index_011 -= N * stride_2;
		index_101 -= N * stride_2;
		index_111 -= N * stride_2;
	}

	// Calculate the scalar values needed for trilinear interpolation
	// These values can be reused for the 9 separate components of the tensor we are looking up
	double w3w2w1, w3w2v1, w3v2w1, w3v2v1, v3w2w1, v3w2v1, v3v2w1, v3v2v1;
	
	w3w2w1 = w3w2v1 = w3 * w2;
	w3w2w1 *= w1;
	w3w2v1 *= v1;
	w3v2w1 = w3v2v1 = w3 * v2;
	w3v2w1 *= w1;
	w3v2v1 *= v1;
	v3w2w1 = v3w2v1 = v3 * w2;
	v3w2w1 *= w1;
	v3w2v1 *= v1;
	v3v2w1 = v3v2v1 = v3 * v2;
	v3v2w1 *= w1;
	v3v2v1 *= v1;


	// We need an additional address offset to go from the address at Field[j,k,l] to the ab component at Field[j,k,l,a,b]
	int offset;

	// Calculate the 9 components, assuming the tensor is symmetric
	offset = 0 * stride_a + 0 * stride_b;
	output.xx = 
					w3w2w1 * (*(double*)(index_000 + offset)) + w3w2v1 * (*(double*)(index_100 + offset)) + 
					w3v2w1 * (*(double*)(index_010 + offset)) + w3v2v1 * (*(double*)(index_110 + offset)) + 
					v3w2w1 * (*(double*)(index_001 + offset)) + v3w2v1 * (*(double*)(index_101 + offset)) +
					v3v2w1 * (*(double*)(index_011 + offset)) + v3v2v1 * (*(double*)(index_111 + offset));

	offset = 0 * stride_a + 1 * stride_b;
	output.xy = output.yx = 
					w3w2w1 * (*(double*)(index_000 + offset)) + w3w2v1 * (*(double*)(index_100 + offset)) + 
					w3v2w1 * (*(double*)(index_010 + offset)) + w3v2v1 * (*(double*)(index_110 + offset)) + 
					v3w2w1 * (*(double*)(index_001 + offset)) + v3w2v1 * (*(double*)(index_101 + offset)) +
					v3v2w1 * (*(double*)(index_011 + offset)) + v3v2v1 * (*(double*)(index_111 + offset));

	offset = 0 * stride_a + 2 * stride_b;
	output.xz = output.zx = 
					w3w2w1 * (*(double*)(index_000 + offset)) + w3w2v1 * (*(double*)(index_100 + offset)) + 
					w3v2w1 * (*(double*)(index_010 + offset)) + w3v2v1 * (*(double*)(index_110 + offset)) + 
					v3w2w1 * (*(double*)(index_001 + offset)) + v3w2v1 * (*(double*)(index_101 + offset)) +
					v3v2w1 * (*(double*)(index_011 + offset)) + v3v2v1 * (*(double*)(index_111 + offset));

	offset = 1 * stride_a + 1 * stride_b;
	output.yy = 
					w3w2w1 * (*(double*)(index_000 + offset)) + w3w2v1 * (*(double*)(index_100 + offset)) + 
					w3v2w1 * (*(double*)(index_010 + offset)) + w3v2v1 * (*(double*)(index_110 + offset)) + 
					v3w2w1 * (*(double*)(index_001 + offset)) + v3w2v1 * (*(double*)(index_101 + offset)) +
					v3v2w1 * (*(double*)(index_011 + offset)) + v3v2v1 * (*(double*)(index_111 + offset));

	offset = 1 * stride_a + 2 * stride_b;
	output.yz = output.zy = 
					w3w2w1 * (*(double*)(index_000 + offset)) + w3w2v1 * (*(double*)(index_100 + offset)) + 
					w3v2w1 * (*(double*)(index_010 + offset)) + w3v2v1 * (*(double*)(index_110 + offset)) + 
					v3w2w1 * (*(double*)(index_001 + offset)) + v3w2v1 * (*(double*)(index_101 + offset)) +
					v3v2w1 * (*(double*)(index_011 + offset)) + v3v2v1 * (*(double*)(index_111 + offset));

	offset = 2 * stride_a + 2 * stride_b;
	output.zz = output.zz = 
					w3w2w1 * (*(double*)(index_000 + offset)) + w3w2v1 * (*(double*)(index_100 + offset)) + 
					w3v2w1 * (*(double*)(index_010 + offset)) + w3v2v1 * (*(double*)(index_110 + offset)) + 
					v3w2w1 * (*(double*)(index_001 + offset)) + v3w2v1 * (*(double*)(index_101 + offset)) +
					v3v2w1 * (*(double*)(index_011 + offset)) + v3v2v1 * (*(double*)(index_111 + offset));

	/*
	for (int a = 0; a < 3; a++)
	{
		for (int b = 0; b < 3; b++)
		{
			offset = a * stride_a + b * stride_b;

			*output.get(3 * a + b) = 
						  w3w2w1 * (*(double*)(index_000 + offset)) + w3w2v1 * (*(double*)(index_100 + offset)) + 
						  w3v2w1 * (*(double*)(index_010 + offset)) + w3v2v1 * (*(double*)(index_110 + offset)) + 
						  v3w2w1 * (*(double*)(index_001 + offset)) + v3w2v1 * (*(double*)(index_101 + offset)) +
						  v3v2w1 * (*(double*)(index_011 + offset)) + v3v2v1 * (*(double*)(index_111 + offset));
		}
	}*/
}

void Lerp(Tensor33 &output, PyArrayObject *Field_xx, PyArrayObject *Field_xy, Vector3 &x)
{
	int i, j, k;
	i = (int)(x.x);
	j = (int)(x.y);
	k = (int)(x.z);
	int N = Field_xx->dimensions[2];

	// Mirror negative indices
	i = abs(i);
	j = abs(j);
	k = abs(k);

	// Calculate interpolation coefficients
	double v1 = x.x - i, v2 = x.y - j, v3 = x.z - k;
	double w1 = 1 - v1, w2 = 1 - v2, w3 = 1 - v3;
	double w3w2w1, w3w2v1, w3v2w1, w3v2v1, v3w2w1, v3w2v1, v3v2w1, v3v2v1;
	
	w3w2w1 = w3w2v1 = w3 * w2;
	w3w2w1 *= w1;
	w3w2v1 *= v1;
	w3v2w1 = w3v2v1 = w3 * v2;
	w3v2w1 *= w1;
	w3v2v1 *= v1;
	v3w2w1 = v3w2v1 = v3 * w2;
	v3w2w1 *= w1;
	v3w2v1 *= v1;
	v3v2w1 = v3v2v1 = v3 * v2;
	v3v2w1 *= w1;
	v3v2v1 *= v1;

	// Calculate array indices
	int stride_xx_0 = Field_xx->strides[0], stride_xx_1 = Field_xx->strides[1], stride_xx_2 = Field_xx->strides[2];
	int stride_xy_0 = Field_xy->strides[0], stride_xy_1 = Field_xy->strides[1], stride_xy_2 = Field_xy->strides[2];
	char *index_000, *index_100, *index_010, *index_001, *index_110, *index_011, *index_101, *index_111;

	// Lookup Field_xx(x)
	index_000 = (Field_xx)->data + i * stride_xx_0 + j * stride_xx_1 + k * stride_xx_2;

	index_100 = index_000 + stride_xx_0;
	index_010 = index_000 + stride_xx_1;
	index_001 = index_000 + stride_xx_2;
	
	index_011 = index_001 + stride_xx_1;
	index_110 = index_010 + stride_xx_0;
	index_101 = index_001 + stride_xx_0;

	index_111 = index_101 + stride_xx_1;

	output.xx = w3w2w1 * (*(double*)index_000) + w3w2v1 * (*(double*)index_100) + 
				w3v2w1 * (*(double*)index_010) + w3v2v1 * (*(double*)index_110) + 
				v3w2w1 * (*(double*)index_001) + v3w2v1 * (*(double*)index_101) +
				v3v2w1 * (*(double*)index_011) + v3v2v1 * (*(double*)index_111);
	

	// Lookup Field_xy(x)
	index_000 = (Field_xy)->data + i * stride_xy_0 + j * stride_xy_1 + k * stride_xy_2;

	index_100 = index_000 + stride_xy_0;
	index_010 = index_000 + stride_xy_1;
	index_001 = index_000 + stride_xy_2;
	
	index_011 = index_001 + stride_xy_1;
	index_110 = index_010 + stride_xy_0;
	index_101 = index_001 + stride_xy_0;

	index_111 = index_101 + stride_xy_1;

	/*
	if (debug)
	{
		fprintf(stderr, "35 lerp = %f\n", *(double*)index_000);
		fprintf(stderr, "35 lerp = %f\n", *(double*)index_001);
		fprintf(stderr, "35 lerp = %f\n", *(double*)index_010);
		fprintf(stderr, "35 lerp = %f\n", *(double*)index_011);
		fprintf(stderr, "35 lerp = %f\n", *(double*)index_100);
		fprintf(stderr, "35 lerp = %f\n", *(double*)index_101);
		fprintf(stderr, "35 lerp = %f\n", *(double*)index_110);
		fprintf(stderr, "35 lerp = %f %f\n", *(double*)index_111, v3v2v1);
	}*/

	output.yx =
	output.xy = w3w2w1 * (*(double*)index_000) + w3w2v1 * (*(double*)index_100) + 
				w3v2w1 * (*(double*)index_010) + w3v2v1 * (*(double*)index_110) + 
				v3w2w1 * (*(double*)index_001) + v3w2v1 * (*(double*)index_101) +
				v3v2w1 * (*(double*)index_011) + v3v2v1 * (*(double*)index_111);
	


	// Lookup Field_xz(x), using Field_xy with y and z coordinates swapped
	index_000 = (Field_xy)->data + i * stride_xy_0 + k * stride_xy_1 + j * stride_xy_2;

	index_100 = index_000 + stride_xy_0;
	index_010 = index_000 + stride_xy_1;
	index_001 = index_000 + stride_xy_2;
	
	index_011 = index_001 + stride_xy_1;
	index_110 = index_010 + stride_xy_0;
	index_101 = index_001 + stride_xy_0;

	index_111 = index_101 + stride_xy_1;

	output.zx =
	output.xz = w3w2w1 * (*(double*)index_000) + w3w2v1 * (*(double*)index_100) + 
				v3w2w1 * (*(double*)index_010) + v3w2v1 * (*(double*)index_110) + 
				w3v2w1 * (*(double*)index_001) + w3v2v1 * (*(double*)index_101) +
				v3v2w1 * (*(double*)index_011) + v3v2v1 * (*(double*)index_111);


				
	// Lookup Field_yy(x), using Field_xx with x and y coordinates swapped
	index_000 = (Field_xx)->data + j * stride_xx_0 + i * stride_xx_1 + k * stride_xx_2;

	index_100 = index_000 + stride_xx_0;
	index_010 = index_000 + stride_xx_1;
	index_001 = index_000 + stride_xx_2;
	
	index_011 = index_001 + stride_xx_1;
	index_110 = index_010 + stride_xx_0;
	index_101 = index_001 + stride_xx_0;

	index_111 = index_101 + stride_xx_1;

	output.yy = w3w2w1 * (*(double*)index_000) + w3v2w1 * (*(double*)index_100) + 
				w3w2v1 * (*(double*)index_010) + w3v2v1 * (*(double*)index_110) + 
				v3w2w1 * (*(double*)index_001) + v3v2w1 * (*(double*)index_101) +
				v3w2v1 * (*(double*)index_011) + v3v2v1 * (*(double*)index_111);



	// Lookup Field_yz(x), using Field_xy with (x,y,z) permuted to (y,z,x)
	index_000 = (Field_xy)->data + j * stride_xy_0 + k * stride_xy_1 + i * stride_xy_2;

	index_100 = index_000 + stride_xy_0;
	index_010 = index_000 + stride_xy_1;
	index_001 = index_000 + stride_xy_2;
	
	index_011 = index_001 + stride_xy_1;
	index_110 = index_010 + stride_xy_0;
	index_101 = index_001 + stride_xy_0;

	index_111 = index_101 + stride_xy_1;

	output.zy = 
	output.yz = w3w2w1 * (*(double*)index_000) + w3v2w1 * (*(double*)index_100) + 
				v3w2w1 * (*(double*)index_010) + v3v2w1 * (*(double*)index_110) + 
				w3w2v1 * (*(double*)index_001) + w3v2v1 * (*(double*)index_101) +
				v3w2v1 * (*(double*)index_011) + v3v2v1 * (*(double*)index_111);


	// Lookup Field_zz(x), using Field_xx with x and z coordinates swapped
	index_000 = (Field_xx)->data + k * stride_xx_0 + j * stride_xx_1 + i * stride_xx_2;

	index_100 = index_000 + stride_xx_0;
	index_010 = index_000 + stride_xx_1;
	index_001 = index_000 + stride_xx_2;
	
	index_011 = index_001 + stride_xx_1;
	index_110 = index_010 + stride_xx_0;
	index_101 = index_001 + stride_xx_0;

	index_111 = index_101 + stride_xx_1;

	output.zz = w3w2w1 * (*(double*)index_000) + v3w2w1 * (*(double*)index_100) + 
				w3v2w1 * (*(double*)index_010) + v3v2w1 * (*(double*)index_110) + 
				w3w2v1 * (*(double*)index_001) + v3w2v1 * (*(double*)index_101) +
				w3v2v1 * (*(double*)index_011) + v3v2v1 * (*(double*)index_111);
}