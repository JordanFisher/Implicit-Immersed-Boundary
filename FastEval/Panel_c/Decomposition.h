struct Decomposition {
	PyObject_HEAD
	int Terms;
	double U_Width;
	PyArrayObject **AB_xx, **AB_xy;
};

static void
Decomposition_dealloc(Decomposition* self)
{
	if (self->AB_xx != NULL)
	{
		for (int i = 0; i < self->Terms; i++)
			Py_XDECREF(self->AB_xx[i]);
		delete[] self->AB_xx;
	}

	if (self->AB_xy != NULL)
	{
		for (int i = 0; i < self->Terms; i++)
			Py_XDECREF(self->AB_xy[i]);
		delete[] self->AB_xy;
	}

	self->ob_type->tp_free((PyObject*)self);
	fprintf(stderr, "Decomposition Deallocated\n");
}

static PyObject *
Decomposition_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Decomposition *self;

    self = (Decomposition *)type->tp_alloc(type, 0);
    if (self != NULL)
	{    
		self->Terms = 0;
		self->U_Width = 0;
		self->AB_xx = self->AB_xy = NULL;
    }

    return (PyObject *)self;
}

static int
Decomposition_init(Decomposition *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|", kwlist))
        return -1; 

    return 0;
}
static PyMemberDef Decomposition_members[] = {
	{"Terms", T_INT, offsetof(Decomposition, Terms), 0,
     "Number of rank 1 terms in the decomposition"},
	{"U_Width", T_INT, offsetof(Decomposition, U_Width), 0,
     "The width of the cube in which sources must lie for the decomposition to be valid"},
	{NULL}  /* Sentinel */
};

static PyObject *
Decomposition_Initialize(Decomposition* self, PyObject *args)
{
	if (!PyArg_ParseTuple(args, "ii", &self->Terms, &self->U_Width)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function Initialization");

		return NULL;
	}

	self->AB_xx = new PyArrayObject*[self->Terms];
	self->AB_xy = new PyArrayObject*[self->Terms];

	for (int i = 0; i < self->Terms; i++)
		self->AB_xx[i] = self->AB_xy[i] = NULL;

	Py_RETURN_NONE;
}

static PyObject *
Decomposition_AB(Decomposition* self, PyArrayObject **AB, int Index)
{
	if (Index < 0 || Index >= self->Terms)
	{
		PyErr_SetString(PyExc_ValueError, "Index must satisfy 0 <= Index < Terms");

		return NULL;
	}

	PyObject *AB_Index = (PyObject*)AB[Index];

	if (AB_Index == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Rank 1 element doesn't exist");

		return NULL;
	}
	
	Py_INCREF(AB_Index);
	
	return AB_Index;
}
static PyObject *
Decomposition_AB_xx(Decomposition* self, PyObject *args)
{
	int Index;
	if (!PyArg_ParseTuple(args, "i", &Index)) {
		PyErr_SetString(PyExc_ValueError, "Function AB_xx expects an integer argument");

		return NULL;
	}

	if (self->AB_xx == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "AB_xx array hasn't been initialized");

		return NULL;
	}

	return Decomposition_AB(self, self->AB_xx, Index);
}
static PyObject *
Decomposition_AB_xy(Decomposition* self, PyObject *args)
{
	int Index;
	if (!PyArg_ParseTuple(args, "i", &Index)) {
		PyErr_SetString(PyExc_ValueError, "Function AB_xy expects an integer argument");

		return NULL;
	}

	if (self->AB_xy == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "AB_xy array hasn't been initialized");

		return NULL;
	}

	return Decomposition_AB(self, self->AB_xy, Index);
}




static PyObject *
Decomposition_Set_AB(Decomposition* self, PyArrayObject **AB, PyArrayObject *NewAB_Index, int Index)
{
	if (Index < 0 || Index >= self->Terms)
	{
		PyErr_SetString(PyExc_ValueError, "Index must satisfy 0 <= Index < Terms");

		return NULL;
	}

	Py_XDECREF(AB[Index]);
	
	AB[Index] = NewAB_Index;	
	Py_XINCREF(NewAB_Index);
	
	Py_RETURN_NONE;
}
static PyObject *
Decomposition_Set_AB_xx(Decomposition* self, PyObject *args)
{
	int Index;
	PyArrayObject *NewAB_Index_xx;
	if (!PyArg_ParseTuple(args, "iO", &Index, &NewAB_Index_xx)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function Set_U");

		return NULL;
	}

	if (self->AB_xx == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "AB_xx array hasn't been initialized");

		return NULL;
	}

	return Decomposition_Set_AB(self, self->AB_xx, NewAB_Index_xx, Index);
}
static PyObject *
Decomposition_Set_AB_xy(Decomposition* self, PyObject *args)
{
	int Index;
	PyArrayObject *NewAB_Index_xy;
	if (!PyArg_ParseTuple(args, "iO", &Index, &NewAB_Index_xy)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function Set_U");

		return NULL;
	}

	if (self->AB_xy == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "AB_xy array hasn't been initialized");

		return NULL;
	}

	return Decomposition_Set_AB(self, self->AB_xy, NewAB_Index_xy, Index);
}


static PyMethodDef Decomposition_methods[] = {
    {"Initialize", (PyCFunction)Decomposition_Initialize, METH_VARARGS,
     "(Terms, U_Width), Allocates space to store the decomposition"
    },
    {"AB_xx", (PyCFunction)Decomposition_AB_xx, METH_VARARGS,
	 "(i), Returns the xx component of the i-th tensor in the array \\{A_k + B_k\\}_{k=0}^{Terms-1}"
    },
    {"Set_AB_xx", (PyCFunction)Decomposition_Set_AB_xx, METH_VARARGS,
	 "(i, NewU), Sets the xx component of the i-th tensor in the array \\{A_k + B_k\\}_{k=0}^{Terms-1}"
    },
    {"AB_xy", (PyCFunction)Decomposition_AB_xy, METH_VARARGS,
	 "(i), Returns the xy component of the i-th tensor in the array \\{A_k + B_k\\}_{k=0}^{Terms-1}"
    },
    {"Set_AB_xy", (PyCFunction)Decomposition_Set_AB_xy, METH_VARARGS,
	 "(i, NewU), Sets the xy component of the i-th tensor in the array \\{A_k + B_k\\}_{k=0}^{Terms-1}"
    },
	{NULL}  /* Sentinel */
};

static PyTypeObject DecompositionType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "Decomposition.Decomposition",             /*tp_name*/
    sizeof(Decomposition),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Decomposition_dealloc, /*tp_dealloc*/
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
    "Decomposition objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Decomposition_methods,             /* tp_methods */
    Decomposition_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Decomposition_init,      /* tp_init */
    0,                         /* tp_alloc */
    Decomposition_new,                 /* tp_new */
};

static PyMethodDef decomposition_module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initDecomposition(void) 
{
    PyObject* m;

    if (PyType_Ready(&DecompositionType) < 0)
        return;

    m = Py_InitModule3("Decomposition", decomposition_module_methods,
                       "Implements the Decomposition class.");

    if (m == NULL)
      return;

    Py_INCREF(&DecompositionType);
    PyModule_AddObject(m, "Decomposition", (PyObject *)&DecompositionType);
}