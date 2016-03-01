ctypedef bint bool

cdef class HOMFLYTerm:
    """HOMFLYTerm(int C, int alpha, int zeta)

    An [immutable] HOMFLY monomial object (i.e. a monomial of the form \\\\(Ca^\\\\alpha{}z^\\\\zeta\\\\))

    HOMFLYTerm support the following operations; given T and S two objects;

    - ``~T`` performs the substitution \\\\(a \mapsto a^{-1}\\\\)
    - ``T+S``, provided the exponents of a and z match, returns the monomial with coefficient the sum of coefficients of T and S
    - ``T*S`` produces the product monomial of the two monomials
    """
    cdef int C, alpha, zeta
    cpdef bool equals(HOMFLYTerm x, HOMFLYTerm y)

cdef class HOMFLYPolynomial:
    """HOMFLYPolynomial(latex_or_terms)

    An [immutable] HOMFLY polynomial object.

    HOMFLYPolynomials support the following operations; given T and S two objects;

    - ``~T`` performs the substitution \\\\(a \\\\mapsto a^{-1}\\\\), i.e. returns the HOMFLY of a mirrored knot
    - ``T+S`` produces the sum of the two polynomials
    - ``T*S`` produces the product of the two polynomials
    - ``T**k`` produces ``T*T*...*T``, provided k is an integer

    HOMFLYPolynomials can be checked for equality with ``==`` or ``!=`` as usual.

    :param latex_or_terms: Either a LaTeX-ish string from llmpoly, or an iterable of triples (C, alpha, zeta), representing monomials \\\\(Ca^\\\\alpha{}z^\\\\zeta\\\\)
    :type latex_or_terms: str or iterable
    """
    cdef readonly tuple terms

