%{
    // Dirty trick which adds some Python object info to a plc_type.
    typedef struct plc_type_w {
        plCurve *p; // Honest pointer to plCurve C object
        PyObject **py_cmps; // Sequence which holds objects for components
        char no_own; // Whether or not this wrapper owns data to its plCurve
    } plCurve_w;

    // Create a new plCurve_w object from a plCurve.
    plCurve_w *plCurve_w_from_plCurve(plCurve *p) {
        int i;
        PyObject *o;
        plCurve_w *ret = malloc(sizeof(plCurve_w));
        ret->p = p;
        ret->py_cmps = PyMem_Malloc(p->nc * sizeof(PyObject*));
        for (i = 0; i < p->nc; i++) {
            o = SWIG_InternalNewPointerObj(p->cp+i,
                                           SWIGTYPE_p_plc_strand_type,
                                           0);
            ret->py_cmps[i] = o;
        }
        ret->no_own = 0;
        return ret;
    }
    plCurve_w *plCurve_w_from_plCurve_NOOWN(plCurve *p) {
        plCurve_w *ret = plCurve_w_from_plCurve(p);
        ret->no_own = false;
        return ret;
    }
%}
