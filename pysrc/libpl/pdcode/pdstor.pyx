from .pdstor cimport *
from .plctopology cimport *
from .pd_storage cimport *
from .diagram cimport PlanarDiagram, PlanarDiagram_wrap

cdef extern from "Python.h":
    cdef FILE* PyFile_AsFile(file_obj)

cdef class PDStorage:
    def __init__(self):
        self.p = pd_new_pdstor()

    def add(self, PlanarDiagram pd, pd_equivalence_t eq=DIAGRAM_ISOTOPY):
        pd_addto_pdstor(self.p, pd.p, eq)

    @classmethod
    def read(cls, f, pd_equivalence_t eq=DIAGRAM_ISOTOPY):
        if not isinstance(f, file):
            if isinstance(f, basestring):
                with open(f, "r") as real_file:
                    return cls.read(real_file)

        return PDStorage_wrap(pd_read_pdstor(PyFile_AsFile(f), eq))

    @staticmethod
    def start_incremental(f):
        if not isinstance(f, file):
            if isinstance(f, basestring):
                with open(f, "w") as real_file:
                    return PDStorage.start_incremental(real_file)
        pd_start_incremental_pdstor(PyFile_AsFile(f))

    def add_to_incremental(self, f, unsigned int nhashes, unsigned int nelts):
        if not isinstance(f, file):
            if isinstance(f, basestring):
                with open(f, "a") as real_file:
                    return self.add_to_incremental(real_file, nhashes, nelts)

        pd_addto_incremental_pdstor(PyFile_AsFile(f), self.p, &nhashes, &nelts)
        return nhashes, nelts

    @staticmethod
    def finish_incremental(f, unsigned int nhashes, unsigned int nelts):
        if not isinstance(f, file):
            if isinstance(f, basestring):
                with open(f, "a") as real_file:
                    return PDStorage.finish_incremental(real_file)

        pd_finish_incremental_pdstor(PyFile_AsFile(f), nhashes, nelts)
        return nhashes, nelts

    def first(self):
        return PlanarDiagram_wrap(pd_stor_firstelt(self.p))

    def next(self):
        return PlanarDiagram_wrap(pd_stor_nextelt(self.p))

    def __dealloc__(self):
        if self.p is not NULL:
            pd_free_pdstor(&self.p)
            self.p = NULL

cdef PDStorage PDStorage_wrap(pd_stor_t *stor):
    cdef PDStorage newobj = PDStorage.__new__(PDStorage)
    newobj.p = stor
    return newobj
