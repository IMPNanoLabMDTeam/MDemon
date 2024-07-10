cimport numpy as cnp
from scipy.sparse import csr_matrix


cdef class Source:
    cdef object _dtype
    cdef int ndim
    cdef cnp.npy_intp dim[2]
    cdef cnp.npy_intp tmp1[1]


    def __cinit__(self, str dtype, cnp.ndarray values, **kwargs):
        if dtype == 'bool':
            self._dtype = cnp.NPY_BOOL
        elif dtype == 'float':
            self._dtype = cnp.NPY_FLOAT32
        elif dtype == 'int':
            self._dtype = cnp.NPY_INT32
        self.dim = values.shape
        self.ndim = values.ndim

cdef class Source1D(Source):
    cdef cnp.ndarray _values
    cdef cnp.npy_intp tmp2[2]

    def __cinit__(self, str dtype, cnp.ndarray values):
        for i in range(self.ndim):
            self.tmp2[i] = self.dim[i]
        self.tmp1[0] = self.dim[0]
        if self.ndim == 2:
            self._values = cnp.PyArray_EMPTY(
                2, self.tmp2, self._dtype, 0)
        else:
            self._values = cnp.PyArray_EMPTY(
                1, self.tmp1, self._dtype, 0)
                
        cnp.PyArray_CopyInto(self._values, values)
    
    @property
    def values(self):
        return self._values

    @values.setter
    def values(self,valix):
        cdef cnp.ndarray[cnp.npy_int32, ndim=1] ix = valix[0]
        cdef cnp.ndarray values
        cdef cnp.npy_intp dim[1]
        dim[0] = valix[1].shape[0]
        values = cnp.PyArray_EMPTY(1, dim, self._dtype, 0)
        cnp.PyArray_CopyInto(values, valix[1])
        self._values[ix] = values

cdef class Source2D(Source):
    cdef cnp.ndarray row
    cdef cnp.ndarray col
    cdef cnp.ndarray data
    cdef object _values
    cdef int _N
    cdef int _M

    def __cinit__(self, str dtype, cnp.ndarray values, **kwargs):
        
        self._N = kwargs.get('N', 0)
        self._M = kwargs.get('M', 0)

        self.tmp1[0] = self.dim[1]

        self.row =  cnp.PyArray_EMPTY(1, self.tmp1, cnp.NPY_INT32, 0)
        cnp.PyArray_CopyInto(self.row, values[0])
        
        self.col =  cnp.PyArray_EMPTY(1, self.tmp1, cnp.NPY_INT32, 0)
        cnp.PyArray_CopyInto(self.col, values[1])       

        self.data =  cnp.PyArray_EMPTY(1, self.tmp1, self._dtype, 0)
        cnp.PyArray_CopyInto(self.data, values[2])


    def __init__(self, str dtype, cnp.ndarray values, **kwargs):
        self._values = csr_matrix((self.data, (self.row, self.col)), shape=(self._N, self._M))


    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, valmn):
        val = csr_matrix((valmn[2], (valmn[0], valmn[1])), 
                            shape=(self.N, self.M))
        self._values += val

    @property
    def N(self):
        return self._N

    @property
    def M(self):
        raise self._M