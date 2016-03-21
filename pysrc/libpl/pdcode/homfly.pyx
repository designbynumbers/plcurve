import re
from operator import itemgetter, mul
from .homfly cimport *

_latex_term_regex = re.compile("(-?[0-9]*)(a\\^\\{-?[0-9]+\\})?(z\\^\\{-?[0-9]+\\})?")
_latex_coef_regex = re.compile("[az]\\^\\{(-?[0-9]+)\\}")

cdef class HOMFLYTerm:
    def __cinit__(self, int C, int alpha, int zeta):
        self.C = C
        self.alpha = alpha
        self.zeta = zeta

    def __invert__(self):
        return HOMFLYTerm(
            self.C,
            -self.alpha,
            self.zeta)
    def __mul__(HOMFLYTerm x, HOMFLYTerm y):
        return HOMFLYTerm(
            x.C * y.C,
            x.alpha + y.alpha,
            x.zeta + y.zeta)

    def __add__(HOMFLYTerm x, HOMFLYTerm y):
        if x.alpha != y.alpha or x.zeta != y.zeta:
            raise Exception("Terms do not match for summation")
        return HOMFLYTerm(
            x.C + y.C,
            x.alpha, x.zeta)

    def __call__(self, a, z):
        return self.C * a**self.alpha * z**self.zeta

    def __nonzero__(self):
        if self.C == 0:
            return 0
        else:
            return 1
    def __cmp__(HOMFLYTerm x, HOMFLYTerm y):
        if x.alpha < y.alpha:
            return -1
        elif x.alpha == y.alpha:
            if x.zeta < y.zeta:
                return -1
            elif x.zeta == y.zeta:
                return 0
            else:
                return 1
        else:
            return 1
    cpdef bool equals(HOMFLYTerm x, HOMFLYTerm y):
        return x.alpha == y.alpha and x.zeta == y.zeta and x.C == y.C
    def __hash__(self):
        return hash((self.alpha, self.zeta))

    def __reduce__(self):
        """Required method to implement pickling of extension types.

        Returns a tuple: (Constructor-ish, Constructor args).
        """
        return (self.__class__, (self.C, self.alpha, self.zeta))


    def __str__(self):
        cdef str cterm = "-" if self.C == -1 else ("" if self.C == 1 else str(self.C))
        if self.alpha and self.zeta:
            return "%sa^{%s}z^{%s}"%(cterm, self.alpha, self.zeta)
        elif self.alpha:
            return "%sa^{%s}"%(cterm, self.alpha)
        elif self.zeta:
            return "%sz^{%s}"%(cterm, self.zeta)
        else:
            return "%s"%self.C
    def __repr__(self):
        cdef str cterm = "-" if self.C == -1 else ("" if self.C == 1 else str(self.C))
        if self.alpha and self.zeta:
            return "%sa^{%s}z^{%s}"%(cterm, self.alpha, self.zeta)
        elif self.alpha:
            return "%sa^{%s}"%(cterm, self.alpha)
        elif self.zeta:
            return "%sz^{%s}"%(cterm, self.zeta)
        else:
            return "%s"%self.C

cdef class HOMFLYPolynomial:
    def __cinit__(self):
        self.terms = None
    def __init__(self, data, bool sort=True):
        """Create a new HOMFLY polynomial from a LaTeX representation"""
        # let's do this pythonically now and then improve it later
        cdef basestring latex

        if isinstance(data, basestring):
            latex = data
            self._load_from_latex(latex, sort=sort)
        elif getattr(data, '__iter__', False):
            self._load_from_terms(data, sort=sort)

    def _load_from_latex(self, basestring latex, bool sort=True):
        cdef object coeff, a_part, z_part, term
        cdef list result = []
        # change '- x' into '+ -x' for consistency
        latex = latex.replace("- ", "+ -")
        for term in latex.split(" + "):
            coeff, a_part, z_part = _latex_term_regex.match(term).groups()
            if coeff == "-":
                coeff = -1
            elif coeff == "":
                coeff = 1
            if a_part:
                a_part = _latex_coef_regex.match(a_part).group(1)
            else:
                a_part = 0
            if z_part:
                z_part = _latex_coef_regex.match(z_part).group(1)
            else:
                z_part = 0

            term = HOMFLYTerm(int(coeff),
                              int(a_part),
                              int(z_part))
            if term:
                result.append(term)
        # TODO: Raise error on parse failure or inappropriate latex rep

        if sort:
            self.terms = tuple(sorted(result))
        else:
            self.terms = tuple(result)

    def _load_from_terms(self, terms, bool sort=True):
        homfly_terms = []
        for term in terms:
            if not isinstance(term, HOMFLYTerm):
                term = HOMFLYTerm(*term)
            try:
                homfly_terms[homfly_terms.index(term)] += term
            except ValueError:
                homfly_terms.append(term)

        if sort:
            self.terms = tuple(sorted(homfly_terms))
        else:
            self.terms = tuple(homfly_terms)

    def __getitem__(self, x):
        try:
            return self.terms[x]
        except IndexError:
            raise
    def __len__(self):
        return len(self.terms)

    def __call__(self, a, z):
        return sum(term(a,z) for term in self.terms)

    def __invert__(HOMFLYPolynomial self):
        cdef HOMFLYPolynomial newpoly = HOMFLYPolynomial.__new__(HOMFLYPolynomial)
        newpoly.terms = tuple(sorted((~A for A in self.terms)))
        return newpoly
    def __mul__(HOMFLYPolynomial x, HOMFLYPolynomial y):
        cdef list newterms = list()
        cdef HOMFLYPolynomial newpoly = HOMFLYPolynomial.__new__(HOMFLYPolynomial)
        cdef HOMFLYTerm A,B,C
        cdef dict found = dict()

        # There is room to replace this loop & sort with an O(n^2) solution
        for A in x.terms:
            for B in y.terms:
                C = A*B
                n = (C.alpha, C.zeta)
                if not C:
                    continue
                if n in found:
                    newterms[found[n]] = newterms[found[n]] + C
                else:
                    found[n] = len(newterms)
                    newterms.append(C)

        # remove any 0-coefficient terms... can possibly do this inline
        newpoly.terms = tuple(sorted(filter(None, newterms)))
        return newpoly
    def __pow__(HOMFLYPolynomial self, int k, z):
        if k == 0:
            return HOMFLYPolynomial('1')
        elif k < 0:
            raise Exception("Power must be positive")
        else:
            return reduce(mul, [self]*k)

    def __richcmp__(HOMFLYPolynomial x, HOMFLYPolynomial y, int op):
        cdef HOMFLYTerm A,B
        if op == 2:
            if not isinstance(x, HOMFLYPolynomial) or not isinstance(y, HOMFLYPolynomial):
                return False
            if len(x.terms) != len(y.terms):
                return False
            return not (False in (A.equals(B) for A,B in zip(x.terms, y.terms)))
        elif op == 3:
            if not isinstance(x, HOMFLYPolynomial) or not isinstance(y, HOMFLYPolynomial):
                return True
            if len(x.terms) != len(y.terms):
                return True
            return False in (A.equals(B) for A,B in zip(x.terms, y.terms))
        else:
            raise NotImplementedError(
                "Only equality checking implemented for HPs")

    def __reduce__(self):
        return (self.__class__, (self.terms,))

    def __hash__(self):
        return hash(self.terms)
    def __str__(self):
        return " + ".join(str(t) for t in self.terms)
    def __repr__(self):
        return " + ".join(repr(t) for t in self.terms)
