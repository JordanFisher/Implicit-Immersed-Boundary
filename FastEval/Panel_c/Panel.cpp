#include<Python.h>
#include "structmember.h"
#include "arrayobject.h"

#include "cPythonUtility.h"
#include "Vector3.h"

#include "Interpolation.h"
#include "Decomposition.h"
#include "DecompGroup.h"

// Returns true if Point is inside the box with bottom left coordinate BL and top right coordinate TR
// NOTE: set up for the periodic domain [0,1]^3
bool PointIsInBox(Vector3 Point, Vector3 BL, Vector3 TR)
{
	Point -= BL;
	TR -= BL;
	BL -= BL;

    Point.x -= floor(Point.x);
	if (Point.x < BL.x || Point.x >= TR.x) return false;

    Point.y -= floor(Point.y);    
	if (Point.y < BL.y || Point.y >= TR.y) return false;    
    
    Point.z -= floor(Point.z);	
    if (Point.z < BL.z || Point.z >= TR.z) return false;
    
	return true;
}

struct Panel;
struct PointData;

// Stores the data necessary to run a FastEval on a single point
struct PointData {
	int NumFarPanels, NumNearPanels;
	Panel **Near;	// The list of panels that must considered when evaluating the force at the point
	Tensor33 **Far;
	Tensor33 **V, **A;

	int Index;

	/// No memory is destroyed or allocated, but indices are reset.
	void Clean()
	{
		NumFarPanels = NumNearPanels = 0;
	}

	void Initialize(int Depth)
	{
		Clean();

		Far = new Tensor33*[(Depth - 2) * 27 * 8];
		V = new Tensor33*[(Depth - 2) * 27 * 8];
		Near = new Panel*[40*3];
		A = new Tensor33*[40*3];
	}

	void InitializeFrom(PointData data)
	{
		NumFarPanels = data.NumFarPanels;
		NumNearPanels = data.NumNearPanels;

		Far = new Tensor33*[NumFarPanels];
		V = new Tensor33*[NumFarPanels];
		Near = new Panel*[NumNearPanels];
		A = new Tensor33*[NumNearPanels];

		for (int i = 0; i < NumFarPanels; i++)
		{
			Far[i] = data.Far[i];
			V[i] = data.V[i];
		}

		for (int i = 0; i < NumNearPanels; i++)
		{
			Near[i] = data.Near[i];
			A[i] = data.A[i];
		}
	}

	/// Shallow deallocation.
	/// Arrays are deleted, but objects comprising arrays are not deallocated.
	void shallow_dealloc()
	{
		delete[] A;
		delete[] V;

		delete[] Far;
		delete[] Near;
	}

	/// Full deallocation.
	void dealloc()
	{
		for (int i = 0; i < NumNearPanels; i++)
			delete[] A[i];

		for (int i = 0; i < NumFarPanels; i++)
			delete[] V[i];

		shallow_dealloc();
	}
};

struct Panel {
    PyObject_HEAD
	Panel **Child;
	bool Childless;
    int Level;

	int MinPointsToHaveChildren;

	PyArrayObject *X;
	int *Index, N;

	Vector3 Center, Size, BL, TR, WellSeparated_BL, WellSeparated_TR;
	Vector3 Corners[8];

	int Terms;
	Tensor33 **UField_FarField, **UField_U;

	// For parent panel only
	PointData *PointDataList;
};

static PyObject *Panel_BL(Panel* self) { return Vector3ToTuple(self->BL); }
static PyObject *Panel_TR(Panel* self) { return Vector3ToTuple(self->TR); }
static PyObject *Panel_WellSeparated_TR(Panel* self) { return Vector3ToTuple(self->WellSeparated_TR); }
static PyObject *Panel_WellSeparated_BL(Panel* self) { return Vector3ToTuple(self->WellSeparated_BL); }
static PyObject *Panel_Center(Panel* self) { return Vector3ToTuple(self->Center); }
static PyObject *Panel_Size(Panel* self) { return Vector3ToTuple(self->Size); }

static void
Panel_dealloc(Panel* self)
{
	Py_DECREF(self->X);

	if (!self->Childless)
	{
		for (int i = 0; i < 8; i++)
			Py_XDECREF(self->Child[i]);
		delete[] self->Child;
	}

	if (self->Index != NULL)
		delete[] self->Index;

	if (self->PointDataList != NULL)
	{
		for (int i = 0; i < self->N; i++)
			self->PointDataList[i].dealloc();
		delete[] self->PointDataList;
	}
	
	if (self->UField_U != NULL)
	{
		for (int i = 0; i < 8; i++)
			delete[] self->UField_U[i];
		delete[] self->UField_U;
	}
	if (self->UField_FarField != NULL)
	{
		for (int i = 0; i < 8; i++)
			delete[] self->UField_FarField[i];
		delete[] self->UField_FarField;
	}

    self->ob_type->tp_free((PyObject*)self);
	//fprintf(stderr, "Panel Deallocated\n");
}

static PyObject *
Panel_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Panel *self;

    self = (Panel *)type->tp_alloc(type, 0);
    if (self != NULL)
	{    
		self->Childless = true;
        self->Level = 0;
    }

    return (PyObject *)self;
}

void
Panel_init_body(Panel *self, int MinPointsToHaveChildren, int MaxDepth, Panel *Parent, PyArrayObject* X, Vector3 BL, Vector3 TR, int Level)
{
	self->X = X;
	Py_INCREF(self->X);

	self->MinPointsToHaveChildren = MinPointsToHaveChildren;

	self->BL = BL;
	self->TR = TR;

	self->Center = (BL + TR) / 2.;
	self->Size = TR - BL;

    self->WellSeparated_BL = BL - self->Size;
	self->WellSeparated_TR = TR + self->Size;

	self->Corners[0] = self->BL;
	self->Corners[1] = self->BL; self->Corners[1].x = self->TR.x;
	self->Corners[2] = self->BL; self->Corners[2].y = self->TR.y;
	self->Corners[3] = self->Corners[1]; self->Corners[3].y = self->TR.y;
	self->Corners[4] = self->Corners[0]; self->Corners[4].z = self->TR.z;
	self->Corners[5] = self->Corners[1]; self->Corners[5].z = self->TR.z;
	self->Corners[6] = self->Corners[2]; self->Corners[6].z = self->TR.z;
	self->Corners[7] = self->Corners[3]; self->Corners[7].z = self->TR.z;

	self->Level = Level;

    if (Parent == NULL)
	{
		self->N = X->dimensions[0];

		// Let Index be an array of integers from 0 to N-1
        self->Index = new int[self->N];
		for (int i = 0; i < self->N; i++)
			self->Index[i] = i;
	}
    else
	{
		// Let Index be an array large enough to potentially store the indices of all the points in the parent panel
        self->N = 0;
        self->Index = new int[Parent->N];

        for (int i = 0; i < Parent->N; i++)
		{
            // Check to see if i-th point in the parent panel is in this panel
            Vector3 x = GetVector3(X, Parent->Index[i]);
            if (PointIsInBox(x, self->BL, self->TR))
			{
                self->Index[self->N] = Parent->Index[i];
                self->N++;
			}
		}
	}

	// If we have enough points create children panels
	if (self->N > MinPointsToHaveChildren && self->Level + 1 < MaxDepth)
	{
		self->Childless = false;
		self->Child = new Panel*[8];
		
		Vector3 Center = self->Center;

		for (int i = 0; i < 8; i++)
			self->Child[i] = (Panel*)Panel_new(self->ob_type, NULL, NULL);

        Panel_init_body(self->Child[0], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(BL.x,     BL.y,     BL.z), Vector3(Center.x, Center.y, Center.z), self->Level + 1);
        Panel_init_body(self->Child[1], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(Center.x, BL.y,     BL.z), Vector3(TR.x,     Center.y, Center.z), self->Level + 1);
        Panel_init_body(self->Child[2], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(BL.x,     Center.y, BL.z), Vector3(Center.x, TR.y,     Center.z), self->Level + 1);
        Panel_init_body(self->Child[3], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(Center.x, Center.y, BL.z), Vector3(TR.x,     TR.y,     Center.z), self->Level + 1);

        Panel_init_body(self->Child[4], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(BL.x,     BL.y,     Center.z), Vector3(Center.x, Center.y, TR.z), self->Level + 1);
        Panel_init_body(self->Child[5], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(Center.x, BL.y,     Center.z), Vector3(TR.x,     Center.y, TR.z), self->Level + 1);
        Panel_init_body(self->Child[6], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(BL.x,     Center.y, Center.z), Vector3(Center.x, TR.y,     TR.z), self->Level + 1);
        Panel_init_body(self->Child[7], MinPointsToHaveChildren, MaxDepth, self, X, Vector3(Center.x, Center.y, Center.z), Vector3(TR.x,     TR.y,     TR.z), self->Level + 1);
	}
	else
		self->Childless = true;
}

static int
Panel_init(Panel *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"MinPointsToHaveChildren", "MaxDepth", "X", "BL", "TR", "Level", NULL};

	PyArrayObject *X, *BL, *TR;
	int Level, MinPointsToHaveChildren, MaxDepth;

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|iiOOOi", kwlist, 
                                      &MinPointsToHaveChildren, &MaxDepth, &X, &BL, &TR, &Level))
        return -1; 

	Panel_init_body(self, MinPointsToHaveChildren, MaxDepth, NULL, X, PyArrayToVector3(BL), PyArrayToVector3(TR), Level);

    return 0;
}

static PyMemberDef Panel_members[] = {
	{"Childless", T_BOOL, offsetof(Panel, Childless), 0,
     "True if Panel has no children"},
	{"Level", T_INT, offsetof(Panel, Level), 0,
     "Panel Level"},
	{"N", T_INT, offsetof(Panel, N), 0,
     "The number of points in the panel"},
	{"X", T_OBJECT_EX, offsetof(Panel, X), 0,
     "The main Nb x 3 array, storing the vector locations of all the immersed structures"},

	{NULL}  /* Sentinel */
};

static PyObject *
Panel_Child(Panel* self, PyObject *args)
{
	int Index;
	if (!PyArg_ParseTuple(args, "i", &Index)) {
		PyErr_SetString(PyExc_ValueError, "Function Child expects an integer argument");

		return NULL;
	}

	if (self->Child == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Panel is childless");

		return NULL;
	}

	if (Index < 0 || Index >= 8)
	{
		PyErr_SetString(PyExc_ValueError, "Index must satisfy 0 <= Index < 8");

		return NULL;
	}

	PyObject *child = (PyObject*)self->Child[Index];

	if (child == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Child Panel does not exist");

		return NULL;
	}
	
	Py_INCREF(child);
	
	return child;
}


void
Panel_CalcUV_body(Panel* self, int Terms, DecompGroup *Lookup)
{
	self->Terms = Terms;

	PyObject *_Decomp = Lookup->UField[self->Level];
    if (_Decomp == Py_None)
	{
        self->UField_FarField = NULL;
        self->UField_U = NULL;
	}
    else
	{
		Decomposition *Decomp = (Decomposition*)_Decomp;

		int N = self->N;
		int Terms = self->Terms;
        self->UField_FarField = new Tensor33*[8];
        self->UField_U = new Tensor33*[8];
		for (int quadrant = 0; quadrant < 8; quadrant++)
		{
	        self->UField_FarField[quadrant] = new Tensor33[Terms];
	        self->UField_U[quadrant] = new Tensor33[N * Terms];
		}
    

		Vector3 dif, vec;
        for (int i = 0; i < N; i++)
		{
			vec = GetVector3(self->X, self->Index[i]);

			for (int quadrant = 0; quadrant < 8; quadrant++)
			{
				// Calculate the difference vector between the source and the panel's corner,
				// then convert that vector to the lookup coordinates in [0,N]^3, assuming the computational domain is [0,1]^3				
				dif = (vec - self->Corners[quadrant]) * Lookup->N;
				dif.Abs();

	            for (int l = 0; l < self->Terms; l++)
				{
					Lerp(self->UField_U[quadrant][i * Terms + l], Decomp->AB_xx[l], Decomp->AB_xy[l], dif);

					// Modidfy the sign of the expansion coefficients according to the anti-symmetries of the Greens function G
					Tensor33 *A = &self->UField_U[quadrant][i * Terms + l];
					if (quadrant == 1 || quadrant == 2 || quadrant == 5 || quadrant == 6)
					{
						A->xy *= -1;
						A->yx = A->xy;
					}
					if (quadrant == 1 || quadrant == 4 || quadrant == 3 || quadrant == 6)
					{
						A->xz *= -1;
						A->zx = A->xz;
					}
					if (quadrant == 2 || quadrant == 4 || quadrant == 3 || quadrant == 5)
					{
						A->yz *= -1;
						A->zy = A->yz;
					}
				}
			}
		}
	}

	if (!self->Childless)
		for (int i = 0; i < 8; i++)
				Panel_CalcUV_body(self->Child[i], self->Terms, Lookup);
}

static PyObject *
Panel_CalcUV(Panel* self, PyObject *args)
{
	DecompGroup *Lookup; // The DecompGroup storing all the rank 1 decompositions
	if (!PyArg_ParseTuple(args, "iO", &self->Terms, &Lookup)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function CalcUV");

		return NULL;
	}

	Panel_CalcUV_body(self, self->Terms, Lookup);

	Py_RETURN_NONE;
}

void
Panel_CalcSeriesTerms_body(Panel* self, PyArrayObject* F)
{
    if (self->UField_U != NULL)
	{
		int Terms = self->Terms;

		for (int i = 0; i < Terms; i++)
			for (int quadrant = 0; quadrant < 8; quadrant++)
				self->UField_FarField[quadrant][i] = 0;
        
		//Tensor33 *UField_U_Address = self->UField_U;
		Vector3 f;
		for (int i = 0; i < self->N; i++)
		{
            f = GetVector3(F, self->Index[i]);
			for (int quadrant = 0; quadrant < 8; quadrant++)
            for (int l = 0; l < Terms; l++)
			{
				/*
				self->UField_FarField[quadrant][l].xx += f.x * (*UField_U_Address).xx;
				self->UField_FarField[quadrant][l].xy += f.x * (*UField_U_Address).xy;
				self->UField_FarField[quadrant][l].xz += f.x * (*UField_U_Address).xz;
				self->UField_FarField[quadrant][l].yx += f.y * (*UField_U_Address).yx;
				self->UField_FarField[quadrant][l].yy += f.y * (*UField_U_Address).yy;
				self->UField_FarField[quadrant][l].yz += f.y * (*UField_U_Address).yz;
				self->UField_FarField[quadrant][l].zx += f.z * (*UField_U_Address).zx;
				self->UField_FarField[quadrant][l].zy += f.z * (*UField_U_Address).zy;
				self->UField_FarField[quadrant][l].zz += f.z * (*UField_U_Address).zz;
				
				UField_U_Address++;
				*/

				Tensor33 U = self->UField_U[quadrant][i * Terms + l];
				Tensor33* H = &self->UField_FarField[quadrant][l];

				H->xx += f.x * U.xx;
				H->xy += f.x * U.xy;
				H->xz += f.x * U.xz;
				H->yx += f.y * U.yx;
				H->yy += f.y * U.yy;
				H->yz += f.y * U.yz;
				H->zx += f.z * U.zx;
				H->zy += f.z * U.zy;
				H->zz += f.z * U.zz;
				/*
				self->UField_FarField[quadrant][l].xx += f.x * U.xx;
				self->UField_FarField[quadrant][l].xy += f.x * U.xy;
				self->UField_FarField[quadrant][l].xz += f.x * U.xz;
				self->UField_FarField[quadrant][l].yx += f.y * U.yx;
				self->UField_FarField[quadrant][l].yy += f.y * U.yy;
				self->UField_FarField[quadrant][l].yz += f.y * U.yz;
				self->UField_FarField[quadrant][l].zx += f.z * U.zx;
				self->UField_FarField[quadrant][l].zy += f.z * U.zy;
				self->UField_FarField[quadrant][l].zz += f.z * U.zz;*/
			}

			//if (f > 0)
			//	fprintf(stderr, "   ! %d %d %d %f\n", self->Index[i], i, self->N, self->UField_U[i * Terms + 0]);
		}

		//for (int l = 0; l < Terms; l++)
		//	fprintf(stderr, "%f\n", self->UField_FarField[l]);
	}

	if (!self->Childless)
		for (int i = 0; i < 8; i++)
				Panel_CalcSeriesTerms_body(self->Child[i], F);
}

static PyObject *
Panel_CalcSeriesTerms(Panel* self, PyObject *args)
{
	PyArrayObject *F; // The vector force distribution on the immersed structure X
	if (!PyArg_ParseTuple(args, "O", &F)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function CalcSeriesTerms");

		return NULL;
	}

	Panel_CalcSeriesTerms_body(self, F);

	Py_RETURN_NONE;
}

Tensor33
Panel_EvalPoint_body(Panel* P, Vector3 &x, PyArrayObject *UField, DecompGroup *Lookup, PyArrayObject *X, PyArrayObject *F)
{
	Tensor33 Sum, Hold;
	Sum = 0;

	if (P->N == 0)
        return Sum;

    // Check to see if we are well Separated from the panel P
    if (!PointIsInBox(x, P->WellSeparated_BL, P->WellSeparated_TR))
	{
		//fprintf(stderr, "Far sum\n");
		Vector3 dif = x - P->Center;
		
		//fprintf(stderr, "x = ");
		//PrintVector3(x);

		//fprintf(stderr, "dif = ");
		//PrintVector3(dif);

		Vector3 x_periodic = x;
		if (dif.x > .5) x_periodic.x -= 1;
		else if (dif.x < -.5) x_periodic.x += 1;
		if (dif.y > .5) x_periodic.y -= 1;
		else if (dif.y < -.5) x_periodic.y += 1;
		if (dif.z > .5) x_periodic.z -= 1;
		else if (dif.z < -.5) x_periodic.z += 1;
		dif = x_periodic - P->Center;

		//fprintf(stderr, "dif = ");
		//PrintVector3(dif);

		// Determine which quadrant x lies in
		int quadrant = 0;
		if (dif.x < 0) quadrant += 1;
		if (dif.y < 0) quadrant += 2;
		if (dif.z < 0) quadrant += 4;

		//fprintf(stderr, "quadrant = %d\n", quadrant);

		dif = (x_periodic - P->Corners[quadrant]) * Lookup->N;
		dif.Abs();

		//fprintf(stderr, "dif = ");
		//PrintVector3(dif);

		//fprintf(stderr, "### dif vector = ");	
		//PrintVector3(dif);

		Decomposition *Decomp = (Decomposition*)Lookup->UField[P->Level];
		//Tensor33 *Address = P->UField_FarField[0];
        for (int l = 0; l < P->Terms; l++)
		{
			// Want to calculate: Sum += P->UField_FarField[l] * Lerp(Decomp->V[l], dif);
			
			Lerp(Hold, Decomp->AB_xx[l], Decomp->AB_xy[l], dif);
			// Optomized version of: Hold *= P->UField_FarField[0][l];
			//Hold *= (*(Address++));
			Hold *= P->UField_FarField[quadrant][l];
            Sum += Hold;
		}
	}
    else
	{		
        // Check to see if P has children
        if (!P->Childless)
		{
			//fprintf(stderr, "To children\n");
			for (int i = 0; i < 8; i++)
			{
				if (P->Child[i]->N > 0)
					Sum += Panel_EvalPoint_body(P->Child[i], x, UField, Lookup, X, F);
			}
		}
        else
		{	
			//fprintf(stderr, "Direct sum\n");
            // The panel is too close, do a direct summation
			int N = P->N;
			int Index;
			for (int i = 0; i < N; i++)
			{
				Index = P->Index[i];
				Vector3 y = (x - GetVector3(X, Index)) * Lookup->N;
				//PrintVector3(y);
				Periodic_Lerp(Hold, UField, y);
				//fprintf(stderr, "Lookup = %f\n", Hold.xy);
				Hold *= GetVector3(F, Index);
                Sum += Hold;
				//fprintf(stderr, "Direct Sum = %f\n", Sum.xy);
			}
		}
	}

    return Sum;
}

// Evaluate the force at x generated by a force distribution F on the immersed structure X
static PyObject *
Panel_EvalPoint(Panel* self, PyObject *args)
{
	PyArrayObject *x;      // The 3-vector where we want to evaluate the force field
	PyArrayObject *UField; // The Green's function
	DecompGroup *Lookup;   // The DecompGroup storing all the lookup tables for the rank 1 decompositions
	PyArrayObject *X;      // The Nb x 3 array storing the immersed structure X
	PyArrayObject *F;      // The vector force distribution on the immersed structure X
	if (!PyArg_ParseTuple(args, "OOOOO", &x, &UField, &Lookup, &X, &F)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function EvalPoint");

		return NULL;
	}

	Tensor33 Sum = Panel_EvalPoint_body(self, PyArrayToVector3(x), UField, Lookup, X, F);

	return Py_BuildValue("d", Sum.xx);
}



void Panel_PreEval_Point(Panel* P, Vector3 &x, PyArrayObject *UField, DecompGroup *Lookup, PyArrayObject *X, PointData &Data)
{
	if (P->N == 0)
        return;

    // Check to see if we are well Separated from the panel P
	if (P->N > P->MinPointsToHaveChildren && !PointIsInBox(x, P->WellSeparated_BL, P->WellSeparated_TR))
	{
		//fprintf(stderr, "Well seperated N = %d\n", P->N);

		//Data.Far[Data.NumFarPanels] = P;		
		Data.V[Data.NumFarPanels] = new Tensor33[P->Terms];


		Vector3 dif = x - P->Center;

		/*
		fprintf(stderr, "x = ");
		PrintVector3(x);
		fprintf(stderr, "dif = ");
		PrintVector3(dif);
		*/
		Vector3 x_periodic = x;
		if (dif.x > .5) x_periodic.x -= 1;
		else if (dif.x < -.5) x_periodic.x += 1;
		if (dif.y > .5) x_periodic.y -= 1;
		else if (dif.y < -.5) x_periodic.y += 1;
		if (dif.z > .5) x_periodic.z -= 1;
		else if (dif.z < -.5) x_periodic.z += 1;
		dif = x_periodic - P->Center;

		//fprintf(stderr, "dif = ");
		//PrintVector3(dif);

		// Determine which quadrant x lies in
		int quadrant = 0;
		if (dif.x < 0) quadrant += 1;
		if (dif.y < 0) quadrant += 2;
		if (dif.z < 0) quadrant += 4;

		// Associate this point with the correct far field expansion
		Data.Far[Data.NumFarPanels] = P->UField_FarField[quadrant];

		dif = (x_periodic - P->Corners[quadrant]) * Lookup->N;
		dif.Abs();

        
		//fprintf(stderr, "*** dif vector = ");	
		//PrintVector3(dif);
		Decomposition *Decomp = (Decomposition*)Lookup->UField[P->Level];
        for (int l = 0; l < P->Terms; l++)
			Lerp(Data.V[Data.NumFarPanels][l], Decomp->AB_xx[l], Decomp->AB_xy[l], dif);

		Data.NumFarPanels++;
	}
    else
	{		
        // Check to see if P has children
		if (!P->Childless && P->N > P->MinPointsToHaveChildren)
		{
			for (int i = 0; i < 8; i++)
				if (P->Child[i]->N > 0)
					Panel_PreEval_Point(P->Child[i], x, UField, Lookup, X, Data);
		}
        else
		{	
            // The panel is too close, do a direct summation
			int N = P->N;
			int Index;
			
			Data.Near[Data.NumNearPanels] = P;
			Data.A[Data.NumNearPanels] = new Tensor33[N];
			
			for (int i = 0; i < N; i++)
			{
				Index = P->Index[i];
				//fprintf(stderr, "i %d Index %d\n", i, Index);
				Vector3 y = (x - GetVector3(X, Index)) * Lookup->N;
				//PrintVector3(y);				
				Periodic_Lerp(Data.A[Data.NumNearPanels][i], UField, y);
				//Data.A[Data.NumNearPanels][i] = DirectLookup(UField, y);
				
				//fprintf(stderr, "lookup %f\n", Lerp(UField, y));
			}
			//fprintf(stderr, "%d / \n", Data.NumNearPanels);
			Data.NumNearPanels++;
		}
	}
}

void OrderPoints(Panel *P, int &i, int *Index)
{
	if (P->Childless)
	{
		for (int j = 0; j < P->N; j++)
			Index[i++] = P->Index[j];
	}
	else
	{
		for (int j = 0; j < 8; j++)
			OrderPoints(P->Child[j], i, Index);
	}
}

// Do the necessary precomputations to allow for efficient evaluations of the tree code
static PyObject *
Panel_PreEval(Panel* self, PyObject *args)
{
	PyArrayObject *UField; // The Green's function
	DecompGroup *Lookup;   // The DecompGroup storing all the lookup tables for the rank 1 decompositions
	PyArrayObject *X;      // The Nb x 3 array storing the immersed structure X
	if (!PyArg_ParseTuple(args, "OOO", &UField, &Lookup, &X)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function EvalPoint");

		return NULL;
	}

	int Nb = X->dimensions[0];
	self->PointDataList = new PointData[Nb];

	PointData DummyData;
	DummyData.Initialize(Lookup->Depth);

	//fprintf(stderr, "\n\n num of points in point data is %d\n", Nb);
	for (int i = 0; i < Nb; i++)
	{
		//fprintf(stderr, "\n\n%d / %d\n", i, Nb);
		//self->PointDataList[i].Initialize(Lookup->Depth);

		//fprintf(stderr, "Done initializing PointData\n", i, Nb);
		DummyData.Clean();
		Panel_PreEval_Point(self, GetVector3(X, i), UField, Lookup, X, DummyData);// self->PointDataList[i]);

		// Copy data from dummy
		self->PointDataList[i].InitializeFrom(DummyData);

		//fprintf(stderr, "\nDone\n");
		//Py_RETURN_NONE;
	}

	// Cleanup dummy
	// Use a shallow deallocation, as deeper objects are now held by the tree
	DummyData.shallow_dealloc();

	// Order points according to panels
	int i = 0;
	OrderPoints(self, i, self->Index);

	Py_RETURN_NONE;
}



// Evaluate the force at each point in X generated by a force distribution F on the immersed structure X
static PyObject *
Panel_FastEval(Panel* self, PyObject *args)
{
	PyArrayObject *X;      // The Nb x 3 array storing the immersed structure X
	PyArrayObject *dX;	   // A Nb x 3 array to hold the result M_n F, the motion induced on X by F
	PyArrayObject *F;      // The vector force distribution on the immersed structure X

	if (!PyArg_ParseTuple(args, "OOO", &X, &F, &dX)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function EvalPoint");

		return NULL;
	}

	int Nb = X->dimensions[0];
	int Terms;
	Tensor33 Sum;
	PointData *Data;
	Tensor33 Hold;
	for (int _i = 0; _i < Nb; _i++)
	{
		int i = self->Index[_i];

		Sum = 0;

		Data = &self->PointDataList[i];

		for (int j = 0; j < Data->NumFarPanels; j++)
		{
			//fprintf(stderr, "Far\n");
			
			Terms = self->Terms;

			Tensor33 *FarAddress = Data->Far[j];
			Tensor33 *VAddress = Data->V[j];
			for (int l = 0; l < Terms; l++)
			{
				// Calculate: Sum += P->UField_FarField[l] * Data->V[j][l];
				
/*				
				Hold = Data->Far[j][l];
				Hold *= Data->V[j][l];
				
				Sum += Hold;
				*/

				// Optomized version of: Sum += Data->Far[j][l] * Data->V[j][l];
				Sum.MultAdd(*(FarAddress++), *(VAddress++));
			}
		}

		for (int j = 0; j < Data->NumNearPanels; j++)
		{		
			int N = Data->Near[j]->N;
			
			int *IndexAddress = Data->Near[j]->Index;
			Tensor33 *AAddress = Data->A[j];
			for (int k = 0; k < N; k++)
			{
				// Calculate: Sum += GetDoubleF, Data->Near[j]->Index[k]) * Data->A[j][k];
				/*
				Hold = Data->A[j][k];
				
				// Optomized version of: Hold *= GetVector3(F, Data->Near[j]->Index[k]);
				Hold *= GetVector3(F, *IndexAddress);

				IndexAddress++;
				Sum += Hold;
				*/
				
				//Sum.MultAdd(Data->A[j][k], GetVector3(F, Data->Near[j]->Index[k]));
				Sum.MultAdd(*(AAddress++), GetVector3(F, *IndexAddress));
				IndexAddress++;
				//Sum.MultAdd(*(AAddress++), GetVector3(F, *(IndexAddress++)));
			}
		}

		GetDouble(dX, i, 0) = Sum.xx + Sum.yx + Sum.zx;
		GetDouble(dX, i, 1) = Sum.xy + Sum.yy + Sum.zy;
		GetDouble(dX, i, 2) = Sum.xz + Sum.yz + Sum.zz;
	}

	Py_RETURN_NONE;
}

static PyMethodDef Panel_methods[] = {
    {"Child", (PyCFunction)Panel_Child, METH_VARARGS,
     "Returns the child of the panel"
    },
    {"BL", (PyCFunction)Panel_BL, METH_NOARGS,
     "Returns the bottom left coordinate of the panel cube"
    },
    {"TR", (PyCFunction)Panel_TR, METH_NOARGS,
     "Returns the top right coordinate of the panel cube"
    },
    {"WellSeparated_BL", (PyCFunction)Panel_WellSeparated_BL, METH_NOARGS,
     "Returns the bottom left coordinate of the region beyond which points are well Separated from the panel"
    },
    {"WellSeparated_TR", (PyCFunction)Panel_WellSeparated_TR, METH_NOARGS,
     "Returns the top right coordinate of the region beyond which points are well Separated from the panel"
    },
    {"Center", (PyCFunction)Panel_Center, METH_NOARGS,
     "Returns the center coordinate of the panel cube"
    },
    {"Size", (PyCFunction)Panel_Size, METH_NOARGS,
     "Returns the size of the panel cube, as a 3-vector"
    },
    {"CalcSeriesTerms", (PyCFunction)Panel_CalcSeriesTerms, METH_VARARGS,
	 "(F) Calculates the far field series associated with each panel, given the force distribution F on the immersed structure X"
    },
    {"CalcUV", (PyCFunction)Panel_CalcUV, METH_VARARGS,
	 "(Terms, Lookup) Calculates the sources values of each source in each panel given the DecompGroup Lookup"
    },
    {"EvalPoint", (PyCFunction)Panel_EvalPoint, METH_VARARGS,
	 "(x, UField, Lookup, X, F) Evaluate the force at x generated by a force distribution F on the immersed structure X"
    },
    {"PreEval", (PyCFunction)Panel_PreEval, METH_VARARGS,
	 "(UField, Lookup, X) Do the necessary precomputations to allow for efficient evaluations of the tree code"
    },
    {"FastEval", (PyCFunction)Panel_FastEval, METH_VARARGS,
	 "(X, F, dF) Evaluate the force at each point in X generated by a force distribution F on the immersed structure X"
    },
	{NULL}  /* Sentinel */
};

static PyTypeObject PanelType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "Panel.Panel",             /*tp_name*/
    sizeof(Panel),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Panel_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_Level*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Panel objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Panel_methods,             /* tp_methods */
    Panel_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Panel_init,      /* tp_init */
    0,                         /* tp_alloc */
    Panel_new,                 /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initPanel(void) 
{
    PyObject* m;

    if (PyType_Ready(&PanelType) < 0)
        return;

    m = Py_InitModule3("Panel", module_methods,
                       "Implents the Panel class.");

    if (m == NULL)
      return;

    Py_INCREF(&PanelType);
    PyModule_AddObject(m, "Panel", (PyObject *)&PanelType);

	initDecomposition();
	initDecompGroup();
}