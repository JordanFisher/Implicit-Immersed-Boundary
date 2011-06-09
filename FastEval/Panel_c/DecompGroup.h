struct DecompGroup {
	PyObject_HEAD
	int Depth, N;
	PyObject **UField;
};

static void
DecompGroup_dealloc(DecompGroup* self)
{
	if (self->UField != NULL)
	{
		for (int i = 0; i < self->Depth; i++)
			Py_XDECREF(self->UField[i]);
		delete[] self->UField;
	}

	self->ob_type->tp_free((PyObject*)self);	
	fprintf(stderr, "DecompGroup Deallocated\n");
}

static PyObject *
DecompGroup_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    DecompGroup *self;

    self = (DecompGroup *)type->tp_alloc(type, 0);
    if (self != NULL)
	{    
		self->Depth = 0;
		self->UField = NULL;
    }

    return (PyObject *)self;
}

static int
DecompGroup_init(DecompGroup *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|", kwlist))
        return -1; 

    return 0;
}
static PyMemberDef DecompGroup_members[] = {
	{"Depth", T_INT, offsetof(DecompGroup, Depth), 0,
     "The number of refinement levels, and the length of the array UField"},
	{"N", T_INT, offsetof(DecompGroup, N), 0,
     "The resolution of the underlying Eulerian grid"},
	{NULL}  /* Sentinel */
};

static PyObject *
DecompGroup_Initialize(DecompGroup* self, PyObject *args)
{
	if (!PyArg_ParseTuple(args, "ii", &self->Depth, &self->N)) {
		PyErr_SetString(PyExc_ValueError, "Function Initialize expects integer arguments");

		return NULL;
	}

	self->UField = new PyObject*[self->Depth];
	
	for (int i = 0; i < self->Depth; i++)
		self->UField[i] = NULL;

	Py_RETURN_NONE;
}

static PyObject *
DecompGroup_UField(DecompGroup* self, PyObject *args)
{
	int Index;
	if (!PyArg_ParseTuple(args, "i", &Index)) {
		PyErr_SetString(PyExc_ValueError, "Function UField expects an integer argument");

		return NULL;
	}

	if (self->UField == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "UField array hasn't been initialized");

		return NULL;
	}

	if (Index < 0 || Index >= self->Depth)
	{
		PyErr_SetString(PyExc_ValueError, "Index must satisfy 0 <= Index < Depth");

		return NULL;
	}

	PyObject *Decomp = (PyObject*)self->UField[Index];

	if (Decomp == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "Decomposition doesn't exist");

		return NULL;
	}
	
	Py_INCREF(Decomp);
	
	return Decomp;
}

static PyObject *
DecompGroup_Set_UField(DecompGroup* self, PyObject *args)
{
	int Index;
	PyObject* NewDecomp;

	if (!PyArg_ParseTuple(args, "iO", &Index, &NewDecomp)) {
		PyErr_SetString(PyExc_ValueError, "Type mismatch in parameter list for function Set_UField");

		return NULL;
	}

	if (self->UField == NULL)
	{
		PyErr_SetString(PyExc_ValueError, "UField array hasn't been initialized");

		return NULL;
	}

	if (Index < 0 || Index >= self->Depth)
	{
		PyErr_SetString(PyExc_ValueError, "Index must satisfy 0 <= Index < Depth");

		return NULL;
	}
	
	Py_XDECREF(self->UField[Index]);
	
	self->UField[Index] = NewDecomp;	
	Py_XINCREF(NewDecomp);
	
	Py_RETURN_NONE;
}

static PyMethodDef DecompGroup_methods[] = {
    {"Initialize", (PyCFunction)DecompGroup_Initialize, METH_VARARGS,
	 "(Depth, N), Allocates space to store the DecompGroup"
    },
    {"UField", (PyCFunction)DecompGroup_UField, METH_VARARGS,
	 "(i), Returns the i-th level UField decomposition"
    },
    {"Set_UField", (PyCFunction)DecompGroup_Set_UField, METH_VARARGS,
	 "(i, NewDecomp), Sets the i-th level UField decomposition to NewDecomp"
    },
	{NULL}  /* Sentinel */
};

static PyTypeObject DecompGroupType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "DecompGroup.DecompGroup",             /*tp_name*/
    sizeof(DecompGroup),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)DecompGroup_dealloc, /*tp_dealloc*/
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
    "DecompGroup objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    DecompGroup_methods,             /* tp_methods */
    DecompGroup_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)DecompGroup_init,      /* tp_init */
    0,                         /* tp_alloc */
    DecompGroup_new,                 /* tp_new */
};

static PyMethodDef DecompGroup_module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initDecompGroup(void) 
{
    PyObject* m;

    if (PyType_Ready(&DecompGroupType) < 0)
        return;

    m = Py_InitModule3("DecompGroup", DecompGroup_module_methods,
                       "Implents the DecompGroup class.");

    if (m == NULL)
      return;

    Py_INCREF(&DecompGroupType);
    PyModule_AddObject(m, "DecompGroup", (PyObject *)&DecompGroupType);
}