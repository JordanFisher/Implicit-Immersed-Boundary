#include "cDerive.h"
#include"cPythonUtility.h"

void eno2(int n, double h, PyArrayObject* u, PyArrayObject* phi, PyArrayObject* f)
{
	double* _u = CopyPyArrayObjectToDoubleArray(u);
	double* _phi = CopyPyArrayObjectToDoubleArray(phi);

	_eno2(n, h, _u, _phi, f);

	delete[] _u; delete[] _phi;
}

void eno3(int n, double h, PyArrayObject* u, PyArrayObject* phi, PyArrayObject* f)
{
	double* _u = CopyPyArrayObjectToDoubleArray(u);
	double* _phi = CopyPyArrayObjectToDoubleArray(phi);

	_eno3(n, h, _u, _phi, f);

	delete[] _u; delete[] _phi;
}

/// Second order Essentially Non-Oscillatory advection
/// Output: f[j] = u(j) * d/dx phi(j)
void _eno2(int N, double h, double* u, double* phi, PyArrayObject* f)
{
	// Chosen indices, depending on the values of u and the slope of phi
	int k1, k2;
	double hb;

	double s1, s2;
	double* flux = new double[N];

	for (int j = 0; j < N; j++)
	{
		// Decide for upwinding
        if (u[j] >= 0)
		{
			k1 = j;
			hb = h;
		}
        else
		{
			k1 = j+1;
			hb = -h;
		}
        
		//k1 = j+1;
		//flux[j] = phi[(k1+N)%N];
		
		// Find smallest slope
        s1 = (phi[(k1+1+N)%N] - phi[(k1+N)%N]) / h;
        s2 = (phi[(k1+N)%N] - phi[(k1-1+N)%N]) / h;
        if (abs(s1) >= abs(s2))
           k2 = k1 - 1;
        else
           k2 = k1;
        
        flux[(j+N)%N] = phi[(k1+N)%N] + 0.5 * hb * (phi[(k2+1+N)%N] - phi[(k2+N)%N]) / h;
	}

	for (int j = 0; j < N; j++)
        GetDouble(f, j) = u[(j+N)%N] * (flux[(j+N)%N] - flux[(j-1+N)%N]) / h;

	// Cleanup
	delete[] flux;
}


void _eno3(int N, double h, double* u, double* phi, PyArrayObject* f)
{
	// Chosen indices, depending on the values of u and the slope of phi
	int k, ks;

	double* D1 = new double[N];
	double* D2 = new double[N];
	double* D3 = new double[N];

	double h2 = h*h;
	double h3 = h2*h;

	double c, cs;

	for (int j = 0; j < N; j++)
		D1[j] = (phi[(j+1+N)%N] - phi[j]) / h;
	for (int j = 0; j < N; j++)
		D2[j] = (D1[j] - D1[(j-1+N)%N]) / (2*h);
	for (int j = 0; j < N; j++)
		D3[j] = (D2[(j+1+N)%N] - D2[j]) / (3*h);

	for (int j = 0; j < N; j++)
	{
		// Decide for upwinding
        if (u[j] >= 0)
			k = j-1;
        else
			k = j;
        

		if (abs(D2[(k+N)%N]) <= abs(D2[(k+1+N)%N]))
		{
			c=D2[(k+N)%N];
			ks=k-1;
		}
		else
		{
			c=D2[(k+1+N)%N];
			ks=k;
		}


        if(abs(D3[(ks+N)%N]) <= abs(D3[(ks+1+N)%N]))
			cs=D3[(ks+N)%N];
        else
			cs=D3[(ks+1+N)%N];

        float der = D1[(k+N)%N] + c * (2 * (j-k) - 1) * h
                          + cs * (3 * pow((j-ks),2.0) - 6 * (j-ks) + 2) * h2;
        
		GetDouble(f, j) = u[j] * der;
	}

	// Cleanup
	delete[] D1; delete[] D2; delete[] D3;
}



void Convect_Centered(int n, double h, PyArrayObject* u, PyArrayObject* phi, PyArrayObject* f)
{
	double* _u = CopyPyArrayObjectToDoubleArray(u);
	double* _phi = CopyPyArrayObjectToDoubleArray(phi);

	_Convect_Centered(n, h, _u, _phi, f);

	delete _u, _phi;
}

/// Second order Essentially Non-Oscillatory advection
/// Output: f[j] = u(j) * d/dx phi(j)
void _Convect_Centered(int N, double h, double* u, double* phi, PyArrayObject* f)
{
	for (int j = 0; j < N; j++)
        GetDouble(f, j) = u[(j+N)%N] * (phi[(j+1+N)%N] - phi[(j-1+N)%N]) / (2*h);
}