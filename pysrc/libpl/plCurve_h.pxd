from gsl_rng cimport *
from libc.stdio cimport FILE
ctypedef bint bool

cdef extern from "plCurve.h":
    ctypedef struct plc_vector:
        double c[3]
    ctypedef struct plc_color:
        pass
    ctypedef struct plc_strand:
        int nv
        bool open
        int cc
        plc_vector* vt
        plc_color* clr
    ctypedef struct plc_constraint:
        pass
    ctypedef struct plc_cst_kind:
        pass
    ctypedef struct plc_vert_quant:
        pass
    ctypedef struct plc_symmetry_group:
        pass
    ctypedef struct plc_spline:
        pass

    ctypedef struct plc_matrix:
        pass
    ctypedef struct plc_symmetry:
        pass

    ctypedef struct plCurve:
        int nc
        plc_strand* cp
        plc_constraint* cst
        plc_vert_quant* quant
        plc_symmetry_group* G

    plc_vector plc_build_vect(const double x, const double y, const double z)

    # -------------------------------------------------------------------
    # plCurve Data Operations
    plCurve *plc_new(const int components,
                     const int * const nv,
                     const bool * const open,
                     const int * const cc)
    void plc_free(plCurve *L)

    void plc_add_component(plCurve *L, const int add_as, const int nv,
                           const bool open, const int cc,
                           const plc_vector * const vt,
                           const plc_color  * const clr)
    void plc_drop_component(plCurve *L, const int cmp)

    void plc_resize_colorbuf(plCurve *L, const int cp, const int cc)
    void plc_set_color(plCurve *L, const plc_color inColor)

    void plc_set_fixed(plCurve * const L,
                       const int cmp,
                       const int vert,
                       const plc_vector point)

    void plc_constrain_to_line(plCurve * const L,
                               const int cmp,
                               const int vert,
                               const int num_verts,
                               const plc_vector tangent,
                               const plc_vector point_on_line)

    void plc_constrain_to_plane(plCurve * const L,
                                const int cmp,
                                const int vert,
                                const int num_verts,
                                const plc_vector normal,
                                const double dist_from_origin)

    void plc_unconstrain(plCurve * const L, const int cmp,
                         const int vert, const int num_verts)
    int plc_remove_constraint(plCurve * const L,
                              const plc_cst_kind kind,
                              const plc_vector vect[])
    void plc_remove_all_constraints(plCurve * const L)

    bool plc_is_constrained(plCurve * const L,int cp, int vt,plc_constraint **constraint)

    plCurve *plc_read(FILE *file,
                      int *error_num,
                      char error_str[],
                      size_t error_str_len)
    void plc_write(FILE *outfile, plCurve * const L)

    void plc_fix_wrap(plCurve * const L)
    void plc_fix_cst(plCurve * const L)

    plCurve *plc_copy(const plCurve * const L)

    plCurve *plc_double_verts(plCurve * L)

    void plc_version( char *version, size_t strlen)

    # -------------------------------------------------------------------
    # plCurve spline package
    plc_spline *plc_spline_new(const int          components,
                               const int  * const ns,
                               const bool * const open,
                               const int  * const cc)
    void plc_spline_free(  plc_spline *L)

    plc_spline *plc_convert_to_spline(plCurve * const L, bool *ok)
    plCurve *plc_convert_from_spline(const plc_spline * const spL,
                                     const int * const nv)

    plc_vector plc_sample_spline(const plc_spline * const spL,
                                 const int cmp,
                                 double s)

    plc_vector plc_spline_tangent(const plc_spline * const spL,
                                  const int cmp,
                                  double s)

    # -------------------------------------------------------------------
    # plCurve Geometric Information
    int plc_num_edges(const plCurve * const L)
    int plc_edges(const plCurve * const L,
                  int *component_edges)
    int plc_num_verts(const plCurve * const L)

    int plc_vertex_num(const plCurve * const L, const int cp, const int vt)
    int plc_cp_num(const plCurve * const L, int wrapVt)
    int plc_vt_num(const plCurve * const L, int wrapVt)

    double plc_turning_angle(plCurve * const L, const int cmp, const int vert,
                             bool *ok)
    double plc_MR_curvature(plCurve * const L, const int cmp, const int vert)
    double plc_totalcurvature(const plCurve * const L,
                              double *component_tc)
    double plc_totaltorsion(const plCurve * const L,
                            double *component_tc)

    float *plc_dihedral_angles(const plCurve * const L,unsigned int *ndihedrals)
    plc_vector plc_mean_tangent(const plCurve * const L, const int cmp,
                                const int vert, bool *ok)

    double plc_arclength(const plCurve * const L,
                         double *component_lengths)
    double plc_subarc_length(const plCurve * const L, const int cmp,
                             const int vert1, const int vert2)
    double plc_s(const plCurve * const L, const int cmp, const int vert)

    void plc_edgelength_stats(const plCurve * const L,
                              double *longest, double *shortest,
                              double *mean, double *moment2)

    double plc_check_cst(const plCurve * const L)

    double plc_pointset_diameter(const plCurve * const L)
    plc_vector plc_center_of_mass(const plCurve * const L)
    double plc_gyradius(const plCurve * const L)
    double *plc_mean_squared_chordlengths( plCurve *L, int cp, int *skips,int nskips)

    cdef struct plc_nearest_neighbor_pc_data:
        plc_vector *check_buffer
        int search_dimension
        int *sorted_buffer

    cdef struct plc_nearest_vertex_pc_data:
        plCurve *check_curve
        plc_nearest_neighbor_pc_data **component_data

    plc_vector plc_nearest_vertex(const plc_vector pt,plCurve *L,int *cp, int *vt,
                                  plc_nearest_vertex_pc_data **pc_data, int *plc_error)

    void plc_nearest_vertex_pc_data_free(plc_nearest_vertex_pc_data **pc_data)

    int plc_nearest_neighbor(const plc_vector pt,const int n, plc_vector *buffer,
                             plc_nearest_neighbor_pc_data **pc_data, int *plc_error)
    void plc_nearest_neighbor_pc_data_free(plc_nearest_neighbor_pc_data **pc_data)

    # -------------------------------------------------------------------
    # plCurve Geometric Operations
    void plc_scale( plCurve *L, const double alpha)
    void plc_whitten(plCurve *L, int mirror, int *eps, int *perm)
    void plc_pfm( plCurve *L, int cp, int vt0, int vt1, double angle)

    void plc_rotate( plCurve *L, plc_vector axis, double angle)
    void plc_random_rotate(plCurve *link, plc_vector axis)
    void plc_translate(plCurve *link,plc_vector translation)
    void plc_perturb( plCurve *L, double radius)
    void plc_project(plCurve *L, plc_vector N)

    plCurve *plc_delete_arc(plCurve *L,int cp,int vt1, int vt2)

    # -------------------------------------------------------------------
    # plCurve Random Polygon Library
    plCurve *plc_random_closed_polygon(gsl_rng *r, int nEdges)
    plCurve *plc_random_open_polygon(gsl_rng *r,int nEdges)

    plCurve *plc_random_closed_plane_polygon(gsl_rng *r,int nEdges)
    plCurve *plc_random_open_plane_polygon(gsl_rng *r,int nEdges)

    plCurve *plc_random_equilateral_closed_polygon(gsl_rng *r,int nEdges)
    plCurve *plc_random_equilateral_open_polygon(gsl_rng *r,int nEdges)

    plCurve *plc_loop_closure(gsl_rng *r,int cp,plCurve *openL,int nEdges)

    # -------------------------------------------------------------------
    # plCurve Symmetry Functions
    void plc_identity_matrix(plc_matrix *A)
    void plc_rotation_matrix(plc_vector axis, double angle,plc_matrix *A)
    void plc_reflection_matrix(plc_vector axis,plc_matrix *A)

    plc_symmetry *plc_symmetry_new(plCurve *model)
    void plc_symmetry_free(plc_symmetry **A)

    plc_symmetry *plc_symmetry_copy(plc_symmetry *A)
    plc_symmetry *plc_build_symmetry(plc_matrix *A,plCurve *L)
    plc_symmetry *plc_compose_symmetries(plc_symmetry *A,plc_symmetry *B)

    plc_symmetry_group *plc_symmetry_group_new(int n)
    void plc_symmetry_group_free(plc_symmetry_group **G)
    plc_symmetry_group *plc_symmetry_group_copy(plc_symmetry_group *G)

    plc_symmetry_group *plc_rotation_group(plCurve *L,plc_vector axis, int n)
    plc_symmetry_group *plc_reflection_group(plCurve *L,plc_vector axis)
    plc_symmetry_group *plc_coordplanes_reflection_group(plCurve *L)

    void plc_symmetrize(plCurve *L)
    void plc_symmetrize_variation(plCurve *L,plc_vector *buffer)

    double plc_symmetry_check(plCurve *L,plc_symmetry *A)
    double plc_symmetry_variation_check(plCurve *L,plc_vector *buffer,plc_symmetry *A)

    double plc_symmetry_group_check(plCurve *L)
    double plc_symmetry_group_variation_check(plCurve *L,plc_vector *buffer)


cdef extern from "plcTopology.h":
    ctypedef struct pd_code_t:
        pass
    pd_code_t *pd_code_from_plCurve(gsl_rng *rng, plCurve *L)
