/*
     pd_permdata.h : This is an automatically generated header file
     produced by the Mathematica notebook /data/generate_permdata.nb
     It should not be hand-edited.
*/

#ifndef __PD_PERMDATA_H__
#define __PD_PERMDATA_H__ 1

#define PD_MAX_PC_PERM 6

typedef struct pd_pc_perm_struct {

    pd_idx_t p[6];

} pd_pc_perm_t; 

pd_idx_t pd_n_pc_perms[7] = {0, 1, 2, 6, 24, 120, 720};

pd_pc_perm_t pdint_pcperm_p1[1] = 
{
{{0}}
};

pd_pc_perm_t pdint_pcperm_p2[2] = 
{
{{0, 1}},
{{1, 0}}
};

pd_pc_perm_t pdint_pcperm_p3[6] = 
{
{{0, 1, 2}},
{{0, 2, 1}},
{{1, 0, 2}},
{{1, 2, 0}},
{{2, 0, 1}},
{{2, 1, 0}}
};

pd_pc_perm_t pdint_pcperm_p4[24] = 
{
{{0, 1, 2, 3}},
{{0, 1, 3, 2}},
{{0, 2, 1, 3}},
{{0, 2, 3, 1}},
{{0, 3, 1, 2}},
{{0, 3, 2, 1}},
{{1, 0, 2, 3}},
{{1, 0, 3, 2}},
{{1, 2, 0, 3}},
{{1, 2, 3, 0}},
{{1, 3, 0, 2}},
{{1, 3, 2, 0}},
{{2, 0, 1, 3}},
{{2, 0, 3, 1}},
{{2, 1, 0, 3}},
{{2, 1, 3, 0}},
{{2, 3, 0, 1}},
{{2, 3, 1, 0}},
{{3, 0, 1, 2}},
{{3, 0, 2, 1}},
{{3, 1, 0, 2}},
{{3, 1, 2, 0}},
{{3, 2, 0, 1}},
{{3, 2, 1, 0}}
};

pd_pc_perm_t pdint_pcperm_p5[120] = 
{
{{0, 1, 2, 3, 4}},
{{0, 1, 2, 4, 3}},
{{0, 1, 3, 2, 4}},
{{0, 1, 3, 4, 2}},
{{0, 1, 4, 2, 3}},
{{0, 1, 4, 3, 2}},
{{0, 2, 1, 3, 4}},
{{0, 2, 1, 4, 3}},
{{0, 2, 3, 1, 4}},
{{0, 2, 3, 4, 1}},
{{0, 2, 4, 1, 3}},
{{0, 2, 4, 3, 1}},
{{0, 3, 1, 2, 4}},
{{0, 3, 1, 4, 2}},
{{0, 3, 2, 1, 4}},
{{0, 3, 2, 4, 1}},
{{0, 3, 4, 1, 2}},
{{0, 3, 4, 2, 1}},
{{0, 4, 1, 2, 3}},
{{0, 4, 1, 3, 2}},
{{0, 4, 2, 1, 3}},
{{0, 4, 2, 3, 1}},
{{0, 4, 3, 1, 2}},
{{0, 4, 3, 2, 1}},
{{1, 0, 2, 3, 4}},
{{1, 0, 2, 4, 3}},
{{1, 0, 3, 2, 4}},
{{1, 0, 3, 4, 2}},
{{1, 0, 4, 2, 3}},
{{1, 0, 4, 3, 2}},
{{1, 2, 0, 3, 4}},
{{1, 2, 0, 4, 3}},
{{1, 2, 3, 0, 4}},
{{1, 2, 3, 4, 0}},
{{1, 2, 4, 0, 3}},
{{1, 2, 4, 3, 0}},
{{1, 3, 0, 2, 4}},
{{1, 3, 0, 4, 2}},
{{1, 3, 2, 0, 4}},
{{1, 3, 2, 4, 0}},
{{1, 3, 4, 0, 2}},
{{1, 3, 4, 2, 0}},
{{1, 4, 0, 2, 3}},
{{1, 4, 0, 3, 2}},
{{1, 4, 2, 0, 3}},
{{1, 4, 2, 3, 0}},
{{1, 4, 3, 0, 2}},
{{1, 4, 3, 2, 0}},
{{2, 0, 1, 3, 4}},
{{2, 0, 1, 4, 3}},
{{2, 0, 3, 1, 4}},
{{2, 0, 3, 4, 1}},
{{2, 0, 4, 1, 3}},
{{2, 0, 4, 3, 1}},
{{2, 1, 0, 3, 4}},
{{2, 1, 0, 4, 3}},
{{2, 1, 3, 0, 4}},
{{2, 1, 3, 4, 0}},
{{2, 1, 4, 0, 3}},
{{2, 1, 4, 3, 0}},
{{2, 3, 0, 1, 4}},
{{2, 3, 0, 4, 1}},
{{2, 3, 1, 0, 4}},
{{2, 3, 1, 4, 0}},
{{2, 3, 4, 0, 1}},
{{2, 3, 4, 1, 0}},
{{2, 4, 0, 1, 3}},
{{2, 4, 0, 3, 1}},
{{2, 4, 1, 0, 3}},
{{2, 4, 1, 3, 0}},
{{2, 4, 3, 0, 1}},
{{2, 4, 3, 1, 0}},
{{3, 0, 1, 2, 4}},
{{3, 0, 1, 4, 2}},
{{3, 0, 2, 1, 4}},
{{3, 0, 2, 4, 1}},
{{3, 0, 4, 1, 2}},
{{3, 0, 4, 2, 1}},
{{3, 1, 0, 2, 4}},
{{3, 1, 0, 4, 2}},
{{3, 1, 2, 0, 4}},
{{3, 1, 2, 4, 0}},
{{3, 1, 4, 0, 2}},
{{3, 1, 4, 2, 0}},
{{3, 2, 0, 1, 4}},
{{3, 2, 0, 4, 1}},
{{3, 2, 1, 0, 4}},
{{3, 2, 1, 4, 0}},
{{3, 2, 4, 0, 1}},
{{3, 2, 4, 1, 0}},
{{3, 4, 0, 1, 2}},
{{3, 4, 0, 2, 1}},
{{3, 4, 1, 0, 2}},
{{3, 4, 1, 2, 0}},
{{3, 4, 2, 0, 1}},
{{3, 4, 2, 1, 0}},
{{4, 0, 1, 2, 3}},
{{4, 0, 1, 3, 2}},
{{4, 0, 2, 1, 3}},
{{4, 0, 2, 3, 1}},
{{4, 0, 3, 1, 2}},
{{4, 0, 3, 2, 1}},
{{4, 1, 0, 2, 3}},
{{4, 1, 0, 3, 2}},
{{4, 1, 2, 0, 3}},
{{4, 1, 2, 3, 0}},
{{4, 1, 3, 0, 2}},
{{4, 1, 3, 2, 0}},
{{4, 2, 0, 1, 3}},
{{4, 2, 0, 3, 1}},
{{4, 2, 1, 0, 3}},
{{4, 2, 1, 3, 0}},
{{4, 2, 3, 0, 1}},
{{4, 2, 3, 1, 0}},
{{4, 3, 0, 1, 2}},
{{4, 3, 0, 2, 1}},
{{4, 3, 1, 0, 2}},
{{4, 3, 1, 2, 0}},
{{4, 3, 2, 0, 1}},
{{4, 3, 2, 1, 0}}
};

pd_pc_perm_t pdint_pcperm_p6[720] = 
{
{{0, 1, 2, 3, 4, 5}},
{{0, 1, 2, 3, 5, 4}},
{{0, 1, 2, 4, 3, 5}},
{{0, 1, 2, 4, 5, 3}},
{{0, 1, 2, 5, 3, 4}},
{{0, 1, 2, 5, 4, 3}},
{{0, 1, 3, 2, 4, 5}},
{{0, 1, 3, 2, 5, 4}},
{{0, 1, 3, 4, 2, 5}},
{{0, 1, 3, 4, 5, 2}},
{{0, 1, 3, 5, 2, 4}},
{{0, 1, 3, 5, 4, 2}},
{{0, 1, 4, 2, 3, 5}},
{{0, 1, 4, 2, 5, 3}},
{{0, 1, 4, 3, 2, 5}},
{{0, 1, 4, 3, 5, 2}},
{{0, 1, 4, 5, 2, 3}},
{{0, 1, 4, 5, 3, 2}},
{{0, 1, 5, 2, 3, 4}},
{{0, 1, 5, 2, 4, 3}},
{{0, 1, 5, 3, 2, 4}},
{{0, 1, 5, 3, 4, 2}},
{{0, 1, 5, 4, 2, 3}},
{{0, 1, 5, 4, 3, 2}},
{{0, 2, 1, 3, 4, 5}},
{{0, 2, 1, 3, 5, 4}},
{{0, 2, 1, 4, 3, 5}},
{{0, 2, 1, 4, 5, 3}},
{{0, 2, 1, 5, 3, 4}},
{{0, 2, 1, 5, 4, 3}},
{{0, 2, 3, 1, 4, 5}},
{{0, 2, 3, 1, 5, 4}},
{{0, 2, 3, 4, 1, 5}},
{{0, 2, 3, 4, 5, 1}},
{{0, 2, 3, 5, 1, 4}},
{{0, 2, 3, 5, 4, 1}},
{{0, 2, 4, 1, 3, 5}},
{{0, 2, 4, 1, 5, 3}},
{{0, 2, 4, 3, 1, 5}},
{{0, 2, 4, 3, 5, 1}},
{{0, 2, 4, 5, 1, 3}},
{{0, 2, 4, 5, 3, 1}},
{{0, 2, 5, 1, 3, 4}},
{{0, 2, 5, 1, 4, 3}},
{{0, 2, 5, 3, 1, 4}},
{{0, 2, 5, 3, 4, 1}},
{{0, 2, 5, 4, 1, 3}},
{{0, 2, 5, 4, 3, 1}},
{{0, 3, 1, 2, 4, 5}},
{{0, 3, 1, 2, 5, 4}},
{{0, 3, 1, 4, 2, 5}},
{{0, 3, 1, 4, 5, 2}},
{{0, 3, 1, 5, 2, 4}},
{{0, 3, 1, 5, 4, 2}},
{{0, 3, 2, 1, 4, 5}},
{{0, 3, 2, 1, 5, 4}},
{{0, 3, 2, 4, 1, 5}},
{{0, 3, 2, 4, 5, 1}},
{{0, 3, 2, 5, 1, 4}},
{{0, 3, 2, 5, 4, 1}},
{{0, 3, 4, 1, 2, 5}},
{{0, 3, 4, 1, 5, 2}},
{{0, 3, 4, 2, 1, 5}},
{{0, 3, 4, 2, 5, 1}},
{{0, 3, 4, 5, 1, 2}},
{{0, 3, 4, 5, 2, 1}},
{{0, 3, 5, 1, 2, 4}},
{{0, 3, 5, 1, 4, 2}},
{{0, 3, 5, 2, 1, 4}},
{{0, 3, 5, 2, 4, 1}},
{{0, 3, 5, 4, 1, 2}},
{{0, 3, 5, 4, 2, 1}},
{{0, 4, 1, 2, 3, 5}},
{{0, 4, 1, 2, 5, 3}},
{{0, 4, 1, 3, 2, 5}},
{{0, 4, 1, 3, 5, 2}},
{{0, 4, 1, 5, 2, 3}},
{{0, 4, 1, 5, 3, 2}},
{{0, 4, 2, 1, 3, 5}},
{{0, 4, 2, 1, 5, 3}},
{{0, 4, 2, 3, 1, 5}},
{{0, 4, 2, 3, 5, 1}},
{{0, 4, 2, 5, 1, 3}},
{{0, 4, 2, 5, 3, 1}},
{{0, 4, 3, 1, 2, 5}},
{{0, 4, 3, 1, 5, 2}},
{{0, 4, 3, 2, 1, 5}},
{{0, 4, 3, 2, 5, 1}},
{{0, 4, 3, 5, 1, 2}},
{{0, 4, 3, 5, 2, 1}},
{{0, 4, 5, 1, 2, 3}},
{{0, 4, 5, 1, 3, 2}},
{{0, 4, 5, 2, 1, 3}},
{{0, 4, 5, 2, 3, 1}},
{{0, 4, 5, 3, 1, 2}},
{{0, 4, 5, 3, 2, 1}},
{{0, 5, 1, 2, 3, 4}},
{{0, 5, 1, 2, 4, 3}},
{{0, 5, 1, 3, 2, 4}},
{{0, 5, 1, 3, 4, 2}},
{{0, 5, 1, 4, 2, 3}},
{{0, 5, 1, 4, 3, 2}},
{{0, 5, 2, 1, 3, 4}},
{{0, 5, 2, 1, 4, 3}},
{{0, 5, 2, 3, 1, 4}},
{{0, 5, 2, 3, 4, 1}},
{{0, 5, 2, 4, 1, 3}},
{{0, 5, 2, 4, 3, 1}},
{{0, 5, 3, 1, 2, 4}},
{{0, 5, 3, 1, 4, 2}},
{{0, 5, 3, 2, 1, 4}},
{{0, 5, 3, 2, 4, 1}},
{{0, 5, 3, 4, 1, 2}},
{{0, 5, 3, 4, 2, 1}},
{{0, 5, 4, 1, 2, 3}},
{{0, 5, 4, 1, 3, 2}},
{{0, 5, 4, 2, 1, 3}},
{{0, 5, 4, 2, 3, 1}},
{{0, 5, 4, 3, 1, 2}},
{{0, 5, 4, 3, 2, 1}},
{{1, 0, 2, 3, 4, 5}},
{{1, 0, 2, 3, 5, 4}},
{{1, 0, 2, 4, 3, 5}},
{{1, 0, 2, 4, 5, 3}},
{{1, 0, 2, 5, 3, 4}},
{{1, 0, 2, 5, 4, 3}},
{{1, 0, 3, 2, 4, 5}},
{{1, 0, 3, 2, 5, 4}},
{{1, 0, 3, 4, 2, 5}},
{{1, 0, 3, 4, 5, 2}},
{{1, 0, 3, 5, 2, 4}},
{{1, 0, 3, 5, 4, 2}},
{{1, 0, 4, 2, 3, 5}},
{{1, 0, 4, 2, 5, 3}},
{{1, 0, 4, 3, 2, 5}},
{{1, 0, 4, 3, 5, 2}},
{{1, 0, 4, 5, 2, 3}},
{{1, 0, 4, 5, 3, 2}},
{{1, 0, 5, 2, 3, 4}},
{{1, 0, 5, 2, 4, 3}},
{{1, 0, 5, 3, 2, 4}},
{{1, 0, 5, 3, 4, 2}},
{{1, 0, 5, 4, 2, 3}},
{{1, 0, 5, 4, 3, 2}},
{{1, 2, 0, 3, 4, 5}},
{{1, 2, 0, 3, 5, 4}},
{{1, 2, 0, 4, 3, 5}},
{{1, 2, 0, 4, 5, 3}},
{{1, 2, 0, 5, 3, 4}},
{{1, 2, 0, 5, 4, 3}},
{{1, 2, 3, 0, 4, 5}},
{{1, 2, 3, 0, 5, 4}},
{{1, 2, 3, 4, 0, 5}},
{{1, 2, 3, 4, 5, 0}},
{{1, 2, 3, 5, 0, 4}},
{{1, 2, 3, 5, 4, 0}},
{{1, 2, 4, 0, 3, 5}},
{{1, 2, 4, 0, 5, 3}},
{{1, 2, 4, 3, 0, 5}},
{{1, 2, 4, 3, 5, 0}},
{{1, 2, 4, 5, 0, 3}},
{{1, 2, 4, 5, 3, 0}},
{{1, 2, 5, 0, 3, 4}},
{{1, 2, 5, 0, 4, 3}},
{{1, 2, 5, 3, 0, 4}},
{{1, 2, 5, 3, 4, 0}},
{{1, 2, 5, 4, 0, 3}},
{{1, 2, 5, 4, 3, 0}},
{{1, 3, 0, 2, 4, 5}},
{{1, 3, 0, 2, 5, 4}},
{{1, 3, 0, 4, 2, 5}},
{{1, 3, 0, 4, 5, 2}},
{{1, 3, 0, 5, 2, 4}},
{{1, 3, 0, 5, 4, 2}},
{{1, 3, 2, 0, 4, 5}},
{{1, 3, 2, 0, 5, 4}},
{{1, 3, 2, 4, 0, 5}},
{{1, 3, 2, 4, 5, 0}},
{{1, 3, 2, 5, 0, 4}},
{{1, 3, 2, 5, 4, 0}},
{{1, 3, 4, 0, 2, 5}},
{{1, 3, 4, 0, 5, 2}},
{{1, 3, 4, 2, 0, 5}},
{{1, 3, 4, 2, 5, 0}},
{{1, 3, 4, 5, 0, 2}},
{{1, 3, 4, 5, 2, 0}},
{{1, 3, 5, 0, 2, 4}},
{{1, 3, 5, 0, 4, 2}},
{{1, 3, 5, 2, 0, 4}},
{{1, 3, 5, 2, 4, 0}},
{{1, 3, 5, 4, 0, 2}},
{{1, 3, 5, 4, 2, 0}},
{{1, 4, 0, 2, 3, 5}},
{{1, 4, 0, 2, 5, 3}},
{{1, 4, 0, 3, 2, 5}},
{{1, 4, 0, 3, 5, 2}},
{{1, 4, 0, 5, 2, 3}},
{{1, 4, 0, 5, 3, 2}},
{{1, 4, 2, 0, 3, 5}},
{{1, 4, 2, 0, 5, 3}},
{{1, 4, 2, 3, 0, 5}},
{{1, 4, 2, 3, 5, 0}},
{{1, 4, 2, 5, 0, 3}},
{{1, 4, 2, 5, 3, 0}},
{{1, 4, 3, 0, 2, 5}},
{{1, 4, 3, 0, 5, 2}},
{{1, 4, 3, 2, 0, 5}},
{{1, 4, 3, 2, 5, 0}},
{{1, 4, 3, 5, 0, 2}},
{{1, 4, 3, 5, 2, 0}},
{{1, 4, 5, 0, 2, 3}},
{{1, 4, 5, 0, 3, 2}},
{{1, 4, 5, 2, 0, 3}},
{{1, 4, 5, 2, 3, 0}},
{{1, 4, 5, 3, 0, 2}},
{{1, 4, 5, 3, 2, 0}},
{{1, 5, 0, 2, 3, 4}},
{{1, 5, 0, 2, 4, 3}},
{{1, 5, 0, 3, 2, 4}},
{{1, 5, 0, 3, 4, 2}},
{{1, 5, 0, 4, 2, 3}},
{{1, 5, 0, 4, 3, 2}},
{{1, 5, 2, 0, 3, 4}},
{{1, 5, 2, 0, 4, 3}},
{{1, 5, 2, 3, 0, 4}},
{{1, 5, 2, 3, 4, 0}},
{{1, 5, 2, 4, 0, 3}},
{{1, 5, 2, 4, 3, 0}},
{{1, 5, 3, 0, 2, 4}},
{{1, 5, 3, 0, 4, 2}},
{{1, 5, 3, 2, 0, 4}},
{{1, 5, 3, 2, 4, 0}},
{{1, 5, 3, 4, 0, 2}},
{{1, 5, 3, 4, 2, 0}},
{{1, 5, 4, 0, 2, 3}},
{{1, 5, 4, 0, 3, 2}},
{{1, 5, 4, 2, 0, 3}},
{{1, 5, 4, 2, 3, 0}},
{{1, 5, 4, 3, 0, 2}},
{{1, 5, 4, 3, 2, 0}},
{{2, 0, 1, 3, 4, 5}},
{{2, 0, 1, 3, 5, 4}},
{{2, 0, 1, 4, 3, 5}},
{{2, 0, 1, 4, 5, 3}},
{{2, 0, 1, 5, 3, 4}},
{{2, 0, 1, 5, 4, 3}},
{{2, 0, 3, 1, 4, 5}},
{{2, 0, 3, 1, 5, 4}},
{{2, 0, 3, 4, 1, 5}},
{{2, 0, 3, 4, 5, 1}},
{{2, 0, 3, 5, 1, 4}},
{{2, 0, 3, 5, 4, 1}},
{{2, 0, 4, 1, 3, 5}},
{{2, 0, 4, 1, 5, 3}},
{{2, 0, 4, 3, 1, 5}},
{{2, 0, 4, 3, 5, 1}},
{{2, 0, 4, 5, 1, 3}},
{{2, 0, 4, 5, 3, 1}},
{{2, 0, 5, 1, 3, 4}},
{{2, 0, 5, 1, 4, 3}},
{{2, 0, 5, 3, 1, 4}},
{{2, 0, 5, 3, 4, 1}},
{{2, 0, 5, 4, 1, 3}},
{{2, 0, 5, 4, 3, 1}},
{{2, 1, 0, 3, 4, 5}},
{{2, 1, 0, 3, 5, 4}},
{{2, 1, 0, 4, 3, 5}},
{{2, 1, 0, 4, 5, 3}},
{{2, 1, 0, 5, 3, 4}},
{{2, 1, 0, 5, 4, 3}},
{{2, 1, 3, 0, 4, 5}},
{{2, 1, 3, 0, 5, 4}},
{{2, 1, 3, 4, 0, 5}},
{{2, 1, 3, 4, 5, 0}},
{{2, 1, 3, 5, 0, 4}},
{{2, 1, 3, 5, 4, 0}},
{{2, 1, 4, 0, 3, 5}},
{{2, 1, 4, 0, 5, 3}},
{{2, 1, 4, 3, 0, 5}},
{{2, 1, 4, 3, 5, 0}},
{{2, 1, 4, 5, 0, 3}},
{{2, 1, 4, 5, 3, 0}},
{{2, 1, 5, 0, 3, 4}},
{{2, 1, 5, 0, 4, 3}},
{{2, 1, 5, 3, 0, 4}},
{{2, 1, 5, 3, 4, 0}},
{{2, 1, 5, 4, 0, 3}},
{{2, 1, 5, 4, 3, 0}},
{{2, 3, 0, 1, 4, 5}},
{{2, 3, 0, 1, 5, 4}},
{{2, 3, 0, 4, 1, 5}},
{{2, 3, 0, 4, 5, 1}},
{{2, 3, 0, 5, 1, 4}},
{{2, 3, 0, 5, 4, 1}},
{{2, 3, 1, 0, 4, 5}},
{{2, 3, 1, 0, 5, 4}},
{{2, 3, 1, 4, 0, 5}},
{{2, 3, 1, 4, 5, 0}},
{{2, 3, 1, 5, 0, 4}},
{{2, 3, 1, 5, 4, 0}},
{{2, 3, 4, 0, 1, 5}},
{{2, 3, 4, 0, 5, 1}},
{{2, 3, 4, 1, 0, 5}},
{{2, 3, 4, 1, 5, 0}},
{{2, 3, 4, 5, 0, 1}},
{{2, 3, 4, 5, 1, 0}},
{{2, 3, 5, 0, 1, 4}},
{{2, 3, 5, 0, 4, 1}},
{{2, 3, 5, 1, 0, 4}},
{{2, 3, 5, 1, 4, 0}},
{{2, 3, 5, 4, 0, 1}},
{{2, 3, 5, 4, 1, 0}},
{{2, 4, 0, 1, 3, 5}},
{{2, 4, 0, 1, 5, 3}},
{{2, 4, 0, 3, 1, 5}},
{{2, 4, 0, 3, 5, 1}},
{{2, 4, 0, 5, 1, 3}},
{{2, 4, 0, 5, 3, 1}},
{{2, 4, 1, 0, 3, 5}},
{{2, 4, 1, 0, 5, 3}},
{{2, 4, 1, 3, 0, 5}},
{{2, 4, 1, 3, 5, 0}},
{{2, 4, 1, 5, 0, 3}},
{{2, 4, 1, 5, 3, 0}},
{{2, 4, 3, 0, 1, 5}},
{{2, 4, 3, 0, 5, 1}},
{{2, 4, 3, 1, 0, 5}},
{{2, 4, 3, 1, 5, 0}},
{{2, 4, 3, 5, 0, 1}},
{{2, 4, 3, 5, 1, 0}},
{{2, 4, 5, 0, 1, 3}},
{{2, 4, 5, 0, 3, 1}},
{{2, 4, 5, 1, 0, 3}},
{{2, 4, 5, 1, 3, 0}},
{{2, 4, 5, 3, 0, 1}},
{{2, 4, 5, 3, 1, 0}},
{{2, 5, 0, 1, 3, 4}},
{{2, 5, 0, 1, 4, 3}},
{{2, 5, 0, 3, 1, 4}},
{{2, 5, 0, 3, 4, 1}},
{{2, 5, 0, 4, 1, 3}},
{{2, 5, 0, 4, 3, 1}},
{{2, 5, 1, 0, 3, 4}},
{{2, 5, 1, 0, 4, 3}},
{{2, 5, 1, 3, 0, 4}},
{{2, 5, 1, 3, 4, 0}},
{{2, 5, 1, 4, 0, 3}},
{{2, 5, 1, 4, 3, 0}},
{{2, 5, 3, 0, 1, 4}},
{{2, 5, 3, 0, 4, 1}},
{{2, 5, 3, 1, 0, 4}},
{{2, 5, 3, 1, 4, 0}},
{{2, 5, 3, 4, 0, 1}},
{{2, 5, 3, 4, 1, 0}},
{{2, 5, 4, 0, 1, 3}},
{{2, 5, 4, 0, 3, 1}},
{{2, 5, 4, 1, 0, 3}},
{{2, 5, 4, 1, 3, 0}},
{{2, 5, 4, 3, 0, 1}},
{{2, 5, 4, 3, 1, 0}},
{{3, 0, 1, 2, 4, 5}},
{{3, 0, 1, 2, 5, 4}},
{{3, 0, 1, 4, 2, 5}},
{{3, 0, 1, 4, 5, 2}},
{{3, 0, 1, 5, 2, 4}},
{{3, 0, 1, 5, 4, 2}},
{{3, 0, 2, 1, 4, 5}},
{{3, 0, 2, 1, 5, 4}},
{{3, 0, 2, 4, 1, 5}},
{{3, 0, 2, 4, 5, 1}},
{{3, 0, 2, 5, 1, 4}},
{{3, 0, 2, 5, 4, 1}},
{{3, 0, 4, 1, 2, 5}},
{{3, 0, 4, 1, 5, 2}},
{{3, 0, 4, 2, 1, 5}},
{{3, 0, 4, 2, 5, 1}},
{{3, 0, 4, 5, 1, 2}},
{{3, 0, 4, 5, 2, 1}},
{{3, 0, 5, 1, 2, 4}},
{{3, 0, 5, 1, 4, 2}},
{{3, 0, 5, 2, 1, 4}},
{{3, 0, 5, 2, 4, 1}},
{{3, 0, 5, 4, 1, 2}},
{{3, 0, 5, 4, 2, 1}},
{{3, 1, 0, 2, 4, 5}},
{{3, 1, 0, 2, 5, 4}},
{{3, 1, 0, 4, 2, 5}},
{{3, 1, 0, 4, 5, 2}},
{{3, 1, 0, 5, 2, 4}},
{{3, 1, 0, 5, 4, 2}},
{{3, 1, 2, 0, 4, 5}},
{{3, 1, 2, 0, 5, 4}},
{{3, 1, 2, 4, 0, 5}},
{{3, 1, 2, 4, 5, 0}},
{{3, 1, 2, 5, 0, 4}},
{{3, 1, 2, 5, 4, 0}},
{{3, 1, 4, 0, 2, 5}},
{{3, 1, 4, 0, 5, 2}},
{{3, 1, 4, 2, 0, 5}},
{{3, 1, 4, 2, 5, 0}},
{{3, 1, 4, 5, 0, 2}},
{{3, 1, 4, 5, 2, 0}},
{{3, 1, 5, 0, 2, 4}},
{{3, 1, 5, 0, 4, 2}},
{{3, 1, 5, 2, 0, 4}},
{{3, 1, 5, 2, 4, 0}},
{{3, 1, 5, 4, 0, 2}},
{{3, 1, 5, 4, 2, 0}},
{{3, 2, 0, 1, 4, 5}},
{{3, 2, 0, 1, 5, 4}},
{{3, 2, 0, 4, 1, 5}},
{{3, 2, 0, 4, 5, 1}},
{{3, 2, 0, 5, 1, 4}},
{{3, 2, 0, 5, 4, 1}},
{{3, 2, 1, 0, 4, 5}},
{{3, 2, 1, 0, 5, 4}},
{{3, 2, 1, 4, 0, 5}},
{{3, 2, 1, 4, 5, 0}},
{{3, 2, 1, 5, 0, 4}},
{{3, 2, 1, 5, 4, 0}},
{{3, 2, 4, 0, 1, 5}},
{{3, 2, 4, 0, 5, 1}},
{{3, 2, 4, 1, 0, 5}},
{{3, 2, 4, 1, 5, 0}},
{{3, 2, 4, 5, 0, 1}},
{{3, 2, 4, 5, 1, 0}},
{{3, 2, 5, 0, 1, 4}},
{{3, 2, 5, 0, 4, 1}},
{{3, 2, 5, 1, 0, 4}},
{{3, 2, 5, 1, 4, 0}},
{{3, 2, 5, 4, 0, 1}},
{{3, 2, 5, 4, 1, 0}},
{{3, 4, 0, 1, 2, 5}},
{{3, 4, 0, 1, 5, 2}},
{{3, 4, 0, 2, 1, 5}},
{{3, 4, 0, 2, 5, 1}},
{{3, 4, 0, 5, 1, 2}},
{{3, 4, 0, 5, 2, 1}},
{{3, 4, 1, 0, 2, 5}},
{{3, 4, 1, 0, 5, 2}},
{{3, 4, 1, 2, 0, 5}},
{{3, 4, 1, 2, 5, 0}},
{{3, 4, 1, 5, 0, 2}},
{{3, 4, 1, 5, 2, 0}},
{{3, 4, 2, 0, 1, 5}},
{{3, 4, 2, 0, 5, 1}},
{{3, 4, 2, 1, 0, 5}},
{{3, 4, 2, 1, 5, 0}},
{{3, 4, 2, 5, 0, 1}},
{{3, 4, 2, 5, 1, 0}},
{{3, 4, 5, 0, 1, 2}},
{{3, 4, 5, 0, 2, 1}},
{{3, 4, 5, 1, 0, 2}},
{{3, 4, 5, 1, 2, 0}},
{{3, 4, 5, 2, 0, 1}},
{{3, 4, 5, 2, 1, 0}},
{{3, 5, 0, 1, 2, 4}},
{{3, 5, 0, 1, 4, 2}},
{{3, 5, 0, 2, 1, 4}},
{{3, 5, 0, 2, 4, 1}},
{{3, 5, 0, 4, 1, 2}},
{{3, 5, 0, 4, 2, 1}},
{{3, 5, 1, 0, 2, 4}},
{{3, 5, 1, 0, 4, 2}},
{{3, 5, 1, 2, 0, 4}},
{{3, 5, 1, 2, 4, 0}},
{{3, 5, 1, 4, 0, 2}},
{{3, 5, 1, 4, 2, 0}},
{{3, 5, 2, 0, 1, 4}},
{{3, 5, 2, 0, 4, 1}},
{{3, 5, 2, 1, 0, 4}},
{{3, 5, 2, 1, 4, 0}},
{{3, 5, 2, 4, 0, 1}},
{{3, 5, 2, 4, 1, 0}},
{{3, 5, 4, 0, 1, 2}},
{{3, 5, 4, 0, 2, 1}},
{{3, 5, 4, 1, 0, 2}},
{{3, 5, 4, 1, 2, 0}},
{{3, 5, 4, 2, 0, 1}},
{{3, 5, 4, 2, 1, 0}},
{{4, 0, 1, 2, 3, 5}},
{{4, 0, 1, 2, 5, 3}},
{{4, 0, 1, 3, 2, 5}},
{{4, 0, 1, 3, 5, 2}},
{{4, 0, 1, 5, 2, 3}},
{{4, 0, 1, 5, 3, 2}},
{{4, 0, 2, 1, 3, 5}},
{{4, 0, 2, 1, 5, 3}},
{{4, 0, 2, 3, 1, 5}},
{{4, 0, 2, 3, 5, 1}},
{{4, 0, 2, 5, 1, 3}},
{{4, 0, 2, 5, 3, 1}},
{{4, 0, 3, 1, 2, 5}},
{{4, 0, 3, 1, 5, 2}},
{{4, 0, 3, 2, 1, 5}},
{{4, 0, 3, 2, 5, 1}},
{{4, 0, 3, 5, 1, 2}},
{{4, 0, 3, 5, 2, 1}},
{{4, 0, 5, 1, 2, 3}},
{{4, 0, 5, 1, 3, 2}},
{{4, 0, 5, 2, 1, 3}},
{{4, 0, 5, 2, 3, 1}},
{{4, 0, 5, 3, 1, 2}},
{{4, 0, 5, 3, 2, 1}},
{{4, 1, 0, 2, 3, 5}},
{{4, 1, 0, 2, 5, 3}},
{{4, 1, 0, 3, 2, 5}},
{{4, 1, 0, 3, 5, 2}},
{{4, 1, 0, 5, 2, 3}},
{{4, 1, 0, 5, 3, 2}},
{{4, 1, 2, 0, 3, 5}},
{{4, 1, 2, 0, 5, 3}},
{{4, 1, 2, 3, 0, 5}},
{{4, 1, 2, 3, 5, 0}},
{{4, 1, 2, 5, 0, 3}},
{{4, 1, 2, 5, 3, 0}},
{{4, 1, 3, 0, 2, 5}},
{{4, 1, 3, 0, 5, 2}},
{{4, 1, 3, 2, 0, 5}},
{{4, 1, 3, 2, 5, 0}},
{{4, 1, 3, 5, 0, 2}},
{{4, 1, 3, 5, 2, 0}},
{{4, 1, 5, 0, 2, 3}},
{{4, 1, 5, 0, 3, 2}},
{{4, 1, 5, 2, 0, 3}},
{{4, 1, 5, 2, 3, 0}},
{{4, 1, 5, 3, 0, 2}},
{{4, 1, 5, 3, 2, 0}},
{{4, 2, 0, 1, 3, 5}},
{{4, 2, 0, 1, 5, 3}},
{{4, 2, 0, 3, 1, 5}},
{{4, 2, 0, 3, 5, 1}},
{{4, 2, 0, 5, 1, 3}},
{{4, 2, 0, 5, 3, 1}},
{{4, 2, 1, 0, 3, 5}},
{{4, 2, 1, 0, 5, 3}},
{{4, 2, 1, 3, 0, 5}},
{{4, 2, 1, 3, 5, 0}},
{{4, 2, 1, 5, 0, 3}},
{{4, 2, 1, 5, 3, 0}},
{{4, 2, 3, 0, 1, 5}},
{{4, 2, 3, 0, 5, 1}},
{{4, 2, 3, 1, 0, 5}},
{{4, 2, 3, 1, 5, 0}},
{{4, 2, 3, 5, 0, 1}},
{{4, 2, 3, 5, 1, 0}},
{{4, 2, 5, 0, 1, 3}},
{{4, 2, 5, 0, 3, 1}},
{{4, 2, 5, 1, 0, 3}},
{{4, 2, 5, 1, 3, 0}},
{{4, 2, 5, 3, 0, 1}},
{{4, 2, 5, 3, 1, 0}},
{{4, 3, 0, 1, 2, 5}},
{{4, 3, 0, 1, 5, 2}},
{{4, 3, 0, 2, 1, 5}},
{{4, 3, 0, 2, 5, 1}},
{{4, 3, 0, 5, 1, 2}},
{{4, 3, 0, 5, 2, 1}},
{{4, 3, 1, 0, 2, 5}},
{{4, 3, 1, 0, 5, 2}},
{{4, 3, 1, 2, 0, 5}},
{{4, 3, 1, 2, 5, 0}},
{{4, 3, 1, 5, 0, 2}},
{{4, 3, 1, 5, 2, 0}},
{{4, 3, 2, 0, 1, 5}},
{{4, 3, 2, 0, 5, 1}},
{{4, 3, 2, 1, 0, 5}},
{{4, 3, 2, 1, 5, 0}},
{{4, 3, 2, 5, 0, 1}},
{{4, 3, 2, 5, 1, 0}},
{{4, 3, 5, 0, 1, 2}},
{{4, 3, 5, 0, 2, 1}},
{{4, 3, 5, 1, 0, 2}},
{{4, 3, 5, 1, 2, 0}},
{{4, 3, 5, 2, 0, 1}},
{{4, 3, 5, 2, 1, 0}},
{{4, 5, 0, 1, 2, 3}},
{{4, 5, 0, 1, 3, 2}},
{{4, 5, 0, 2, 1, 3}},
{{4, 5, 0, 2, 3, 1}},
{{4, 5, 0, 3, 1, 2}},
{{4, 5, 0, 3, 2, 1}},
{{4, 5, 1, 0, 2, 3}},
{{4, 5, 1, 0, 3, 2}},
{{4, 5, 1, 2, 0, 3}},
{{4, 5, 1, 2, 3, 0}},
{{4, 5, 1, 3, 0, 2}},
{{4, 5, 1, 3, 2, 0}},
{{4, 5, 2, 0, 1, 3}},
{{4, 5, 2, 0, 3, 1}},
{{4, 5, 2, 1, 0, 3}},
{{4, 5, 2, 1, 3, 0}},
{{4, 5, 2, 3, 0, 1}},
{{4, 5, 2, 3, 1, 0}},
{{4, 5, 3, 0, 1, 2}},
{{4, 5, 3, 0, 2, 1}},
{{4, 5, 3, 1, 0, 2}},
{{4, 5, 3, 1, 2, 0}},
{{4, 5, 3, 2, 0, 1}},
{{4, 5, 3, 2, 1, 0}},
{{5, 0, 1, 2, 3, 4}},
{{5, 0, 1, 2, 4, 3}},
{{5, 0, 1, 3, 2, 4}},
{{5, 0, 1, 3, 4, 2}},
{{5, 0, 1, 4, 2, 3}},
{{5, 0, 1, 4, 3, 2}},
{{5, 0, 2, 1, 3, 4}},
{{5, 0, 2, 1, 4, 3}},
{{5, 0, 2, 3, 1, 4}},
{{5, 0, 2, 3, 4, 1}},
{{5, 0, 2, 4, 1, 3}},
{{5, 0, 2, 4, 3, 1}},
{{5, 0, 3, 1, 2, 4}},
{{5, 0, 3, 1, 4, 2}},
{{5, 0, 3, 2, 1, 4}},
{{5, 0, 3, 2, 4, 1}},
{{5, 0, 3, 4, 1, 2}},
{{5, 0, 3, 4, 2, 1}},
{{5, 0, 4, 1, 2, 3}},
{{5, 0, 4, 1, 3, 2}},
{{5, 0, 4, 2, 1, 3}},
{{5, 0, 4, 2, 3, 1}},
{{5, 0, 4, 3, 1, 2}},
{{5, 0, 4, 3, 2, 1}},
{{5, 1, 0, 2, 3, 4}},
{{5, 1, 0, 2, 4, 3}},
{{5, 1, 0, 3, 2, 4}},
{{5, 1, 0, 3, 4, 2}},
{{5, 1, 0, 4, 2, 3}},
{{5, 1, 0, 4, 3, 2}},
{{5, 1, 2, 0, 3, 4}},
{{5, 1, 2, 0, 4, 3}},
{{5, 1, 2, 3, 0, 4}},
{{5, 1, 2, 3, 4, 0}},
{{5, 1, 2, 4, 0, 3}},
{{5, 1, 2, 4, 3, 0}},
{{5, 1, 3, 0, 2, 4}},
{{5, 1, 3, 0, 4, 2}},
{{5, 1, 3, 2, 0, 4}},
{{5, 1, 3, 2, 4, 0}},
{{5, 1, 3, 4, 0, 2}},
{{5, 1, 3, 4, 2, 0}},
{{5, 1, 4, 0, 2, 3}},
{{5, 1, 4, 0, 3, 2}},
{{5, 1, 4, 2, 0, 3}},
{{5, 1, 4, 2, 3, 0}},
{{5, 1, 4, 3, 0, 2}},
{{5, 1, 4, 3, 2, 0}},
{{5, 2, 0, 1, 3, 4}},
{{5, 2, 0, 1, 4, 3}},
{{5, 2, 0, 3, 1, 4}},
{{5, 2, 0, 3, 4, 1}},
{{5, 2, 0, 4, 1, 3}},
{{5, 2, 0, 4, 3, 1}},
{{5, 2, 1, 0, 3, 4}},
{{5, 2, 1, 0, 4, 3}},
{{5, 2, 1, 3, 0, 4}},
{{5, 2, 1, 3, 4, 0}},
{{5, 2, 1, 4, 0, 3}},
{{5, 2, 1, 4, 3, 0}},
{{5, 2, 3, 0, 1, 4}},
{{5, 2, 3, 0, 4, 1}},
{{5, 2, 3, 1, 0, 4}},
{{5, 2, 3, 1, 4, 0}},
{{5, 2, 3, 4, 0, 1}},
{{5, 2, 3, 4, 1, 0}},
{{5, 2, 4, 0, 1, 3}},
{{5, 2, 4, 0, 3, 1}},
{{5, 2, 4, 1, 0, 3}},
{{5, 2, 4, 1, 3, 0}},
{{5, 2, 4, 3, 0, 1}},
{{5, 2, 4, 3, 1, 0}},
{{5, 3, 0, 1, 2, 4}},
{{5, 3, 0, 1, 4, 2}},
{{5, 3, 0, 2, 1, 4}},
{{5, 3, 0, 2, 4, 1}},
{{5, 3, 0, 4, 1, 2}},
{{5, 3, 0, 4, 2, 1}},
{{5, 3, 1, 0, 2, 4}},
{{5, 3, 1, 0, 4, 2}},
{{5, 3, 1, 2, 0, 4}},
{{5, 3, 1, 2, 4, 0}},
{{5, 3, 1, 4, 0, 2}},
{{5, 3, 1, 4, 2, 0}},
{{5, 3, 2, 0, 1, 4}},
{{5, 3, 2, 0, 4, 1}},
{{5, 3, 2, 1, 0, 4}},
{{5, 3, 2, 1, 4, 0}},
{{5, 3, 2, 4, 0, 1}},
{{5, 3, 2, 4, 1, 0}},
{{5, 3, 4, 0, 1, 2}},
{{5, 3, 4, 0, 2, 1}},
{{5, 3, 4, 1, 0, 2}},
{{5, 3, 4, 1, 2, 0}},
{{5, 3, 4, 2, 0, 1}},
{{5, 3, 4, 2, 1, 0}},
{{5, 4, 0, 1, 2, 3}},
{{5, 4, 0, 1, 3, 2}},
{{5, 4, 0, 2, 1, 3}},
{{5, 4, 0, 2, 3, 1}},
{{5, 4, 0, 3, 1, 2}},
{{5, 4, 0, 3, 2, 1}},
{{5, 4, 1, 0, 2, 3}},
{{5, 4, 1, 0, 3, 2}},
{{5, 4, 1, 2, 0, 3}},
{{5, 4, 1, 2, 3, 0}},
{{5, 4, 1, 3, 0, 2}},
{{5, 4, 1, 3, 2, 0}},
{{5, 4, 2, 0, 1, 3}},
{{5, 4, 2, 0, 3, 1}},
{{5, 4, 2, 1, 0, 3}},
{{5, 4, 2, 1, 3, 0}},
{{5, 4, 2, 3, 0, 1}},
{{5, 4, 2, 3, 1, 0}},
{{5, 4, 3, 0, 1, 2}},
{{5, 4, 3, 0, 2, 1}},
{{5, 4, 3, 1, 0, 2}},
{{5, 4, 3, 1, 2, 0}},
{{5, 4, 3, 2, 0, 1}},
{{5, 4, 3, 2, 1, 0}}
};

pd_pc_perm_t* pd_pc_perm[7] = {NULL, &(pdint_pcperm_p1[0]), &(pdint_pcperm_p2[0]), &(pdint_pcperm_p3[0]), &(pdint_pcperm_p4[0]), &(pdint_pcperm_p5[0]), &(pdint_pcperm_p6[0])};

typedef char pd_perm_pp_t[18];

pd_perm_pp_t pdint_pcperm_print1[1] = 
{
{"e"}
};

pd_perm_pp_t pdint_pcperm_print2[2] = 
{
{"e"},
{"(0 1)"}
};

pd_perm_pp_t pdint_pcperm_print3[6] = 
{
{"e"},
{"(1 2)"},
{"(0 1)"},
{"(0 2 1)"},
{"(0 1 2)"},
{"(0 2)"}
};

pd_perm_pp_t pdint_pcperm_print4[24] = 
{
{"e"},
{"(2 3)"},
{"(1 2)"},
{"(1 3 2)"},
{"(1 2 3)"},
{"(1 3)"},
{"(0 1)"},
{"(0 1)(2 3)"},
{"(0 2 1)"},
{"(0 3 2 1)"},
{"(0 2 3 1)"},
{"(0 3 1)"},
{"(0 1 2)"},
{"(0 1 3 2)"},
{"(0 2)"},
{"(0 3 2)"},
{"(0 2)(1 3)"},
{"(0 3 1 2)"},
{"(0 1 2 3)"},
{"(0 1 3)"},
{"(0 2 3)"},
{"(0 3)"},
{"(0 2 1 3)"},
{"(0 3)(1 2)"}
};

pd_perm_pp_t pdint_pcperm_print5[120] = 
{
{"e"},
{"(3 4)"},
{"(2 3)"},
{"(2 4 3)"},
{"(2 3 4)"},
{"(2 4)"},
{"(1 2)"},
{"(1 2)(3 4)"},
{"(1 3 2)"},
{"(1 4 3 2)"},
{"(1 3 4 2)"},
{"(1 4 2)"},
{"(1 2 3)"},
{"(1 2 4 3)"},
{"(1 3)"},
{"(1 4 3)"},
{"(1 3)(2 4)"},
{"(1 4 2 3)"},
{"(1 2 3 4)"},
{"(1 2 4)"},
{"(1 3 4)"},
{"(1 4)"},
{"(1 3 2 4)"},
{"(1 4)(2 3)"},
{"(0 1)"},
{"(0 1)(3 4)"},
{"(0 1)(2 3)"},
{"(0 1)(2 4 3)"},
{"(0 1)(2 3 4)"},
{"(0 1)(2 4)"},
{"(0 2 1)"},
{"(0 2 1)(3 4)"},
{"(0 3 2 1)"},
{"(0 4 3 2 1)"},
{"(0 3 4 2 1)"},
{"(0 4 2 1)"},
{"(0 2 3 1)"},
{"(0 2 4 3 1)"},
{"(0 3 1)"},
{"(0 4 3 1)"},
{"(0 3 1)(2 4)"},
{"(0 4 2 3 1)"},
{"(0 2 3 4 1)"},
{"(0 2 4 1)"},
{"(0 3 4 1)"},
{"(0 4 1)"},
{"(0 3 2 4 1)"},
{"(0 4 1)(2 3)"},
{"(0 1 2)"},
{"(0 1 2)(3 4)"},
{"(0 1 3 2)"},
{"(0 1 4 3 2)"},
{"(0 1 3 4 2)"},
{"(0 1 4 2)"},
{"(0 2)"},
{"(0 2)(3 4)"},
{"(0 3 2)"},
{"(0 4 3 2)"},
{"(0 3 4 2)"},
{"(0 4 2)"},
{"(0 2)(1 3)"},
{"(0 2)(1 4 3)"},
{"(0 3 1 2)"},
{"(0 4 3 1 2)"},
{"(0 3 1 4 2)"},
{"(0 4 2)(1 3)"},
{"(0 2)(1 3 4)"},
{"(0 2)(1 4)"},
{"(0 3 4 1 2)"},
{"(0 4 1 2)"},
{"(0 3 2)(1 4)"},
{"(0 4 1 3 2)"},
{"(0 1 2 3)"},
{"(0 1 2 4 3)"},
{"(0 1 3)"},
{"(0 1 4 3)"},
{"(0 1 3)(2 4)"},
{"(0 1 4 2 3)"},
{"(0 2 3)"},
{"(0 2 4 3)"},
{"(0 3)"},
{"(0 4 3)"},
{"(0 3)(2 4)"},
{"(0 4 2 3)"},
{"(0 2 1 3)"},
{"(0 2 1 4 3)"},
{"(0 3)(1 2)"},
{"(0 4 3)(1 2)"},
{"(0 3)(1 4 2)"},
{"(0 4 2 1 3)"},
{"(0 2 4 1 3)"},
{"(0 2 3)(1 4)"},
{"(0 3)(1 2 4)"},
{"(0 4 1 2 3)"},
{"(0 3)(1 4)"},
{"(0 4 1 3)"},
{"(0 1 2 3 4)"},
{"(0 1 2 4)"},
{"(0 1 3 4)"},
{"(0 1 4)"},
{"(0 1 3 2 4)"},
{"(0 1 4)(2 3)"},
{"(0 2 3 4)"},
{"(0 2 4)"},
{"(0 3 4)"},
{"(0 4)"},
{"(0 3 2 4)"},
{"(0 4)(2 3)"},
{"(0 2 1 3 4)"},
{"(0 2 1 4)"},
{"(0 3 4)(1 2)"},
{"(0 4)(1 2)"},
{"(0 3 2 1 4)"},
{"(0 4)(1 3 2)"},
{"(0 2 4)(1 3)"},
{"(0 2 3 1 4)"},
{"(0 3 1 2 4)"},
{"(0 4)(1 2 3)"},
{"(0 3 1 4)"},
{"(0 4)(1 3)"}
};

pd_perm_pp_t pdint_pcperm_print6[720] = 
{
{"e"},
{"(4 5)"},
{"(3 4)"},
{"(3 5 4)"},
{"(3 4 5)"},
{"(3 5)"},
{"(2 3)"},
{"(2 3)(4 5)"},
{"(2 4 3)"},
{"(2 5 4 3)"},
{"(2 4 5 3)"},
{"(2 5 3)"},
{"(2 3 4)"},
{"(2 3 5 4)"},
{"(2 4)"},
{"(2 5 4)"},
{"(2 4)(3 5)"},
{"(2 5 3 4)"},
{"(2 3 4 5)"},
{"(2 3 5)"},
{"(2 4 5)"},
{"(2 5)"},
{"(2 4 3 5)"},
{"(2 5)(3 4)"},
{"(1 2)"},
{"(1 2)(4 5)"},
{"(1 2)(3 4)"},
{"(1 2)(3 5 4)"},
{"(1 2)(3 4 5)"},
{"(1 2)(3 5)"},
{"(1 3 2)"},
{"(1 3 2)(4 5)"},
{"(1 4 3 2)"},
{"(1 5 4 3 2)"},
{"(1 4 5 3 2)"},
{"(1 5 3 2)"},
{"(1 3 4 2)"},
{"(1 3 5 4 2)"},
{"(1 4 2)"},
{"(1 5 4 2)"},
{"(1 4 2)(3 5)"},
{"(1 5 3 4 2)"},
{"(1 3 4 5 2)"},
{"(1 3 5 2)"},
{"(1 4 5 2)"},
{"(1 5 2)"},
{"(1 4 3 5 2)"},
{"(1 5 2)(3 4)"},
{"(1 2 3)"},
{"(1 2 3)(4 5)"},
{"(1 2 4 3)"},
{"(1 2 5 4 3)"},
{"(1 2 4 5 3)"},
{"(1 2 5 3)"},
{"(1 3)"},
{"(1 3)(4 5)"},
{"(1 4 3)"},
{"(1 5 4 3)"},
{"(1 4 5 3)"},
{"(1 5 3)"},
{"(1 3)(2 4)"},
{"(1 3)(2 5 4)"},
{"(1 4 2 3)"},
{"(1 5 4 2 3)"},
{"(1 4 2 5 3)"},
{"(1 5 3)(2 4)"},
{"(1 3)(2 4 5)"},
{"(1 3)(2 5)"},
{"(1 4 5 2 3)"},
{"(1 5 2 3)"},
{"(1 4 3)(2 5)"},
{"(1 5 2 4 3)"},
{"(1 2 3 4)"},
{"(1 2 3 5 4)"},
{"(1 2 4)"},
{"(1 2 5 4)"},
{"(1 2 4)(3 5)"},
{"(1 2 5 3 4)"},
{"(1 3 4)"},
{"(1 3 5 4)"},
{"(1 4)"},
{"(1 5 4)"},
{"(1 4)(3 5)"},
{"(1 5 3 4)"},
{"(1 3 2 4)"},
{"(1 3 2 5 4)"},
{"(1 4)(2 3)"},
{"(1 5 4)(2 3)"},
{"(1 4)(2 5 3)"},
{"(1 5 3 2 4)"},
{"(1 3 5 2 4)"},
{"(1 3 4)(2 5)"},
{"(1 4)(2 3 5)"},
{"(1 5 2 3 4)"},
{"(1 4)(2 5)"},
{"(1 5 2 4)"},
{"(1 2 3 4 5)"},
{"(1 2 3 5)"},
{"(1 2 4 5)"},
{"(1 2 5)"},
{"(1 2 4 3 5)"},
{"(1 2 5)(3 4)"},
{"(1 3 4 5)"},
{"(1 3 5)"},
{"(1 4 5)"},
{"(1 5)"},
{"(1 4 3 5)"},
{"(1 5)(3 4)"},
{"(1 3 2 4 5)"},
{"(1 3 2 5)"},
{"(1 4 5)(2 3)"},
{"(1 5)(2 3)"},
{"(1 4 3 2 5)"},
{"(1 5)(2 4 3)"},
{"(1 3 5)(2 4)"},
{"(1 3 4 2 5)"},
{"(1 4 2 3 5)"},
{"(1 5)(2 3 4)"},
{"(1 4 2 5)"},
{"(1 5)(2 4)"},
{"(0 1)"},
{"(0 1)(4 5)"},
{"(0 1)(3 4)"},
{"(0 1)(3 5 4)"},
{"(0 1)(3 4 5)"},
{"(0 1)(3 5)"},
{"(0 1)(2 3)"},
{"(0 1)(2 3)(4 5)"},
{"(0 1)(2 4 3)"},
{"(0 1)(2 5 4 3)"},
{"(0 1)(2 4 5 3)"},
{"(0 1)(2 5 3)"},
{"(0 1)(2 3 4)"},
{"(0 1)(2 3 5 4)"},
{"(0 1)(2 4)"},
{"(0 1)(2 5 4)"},
{"(0 1)(2 4)(3 5)"},
{"(0 1)(2 5 3 4)"},
{"(0 1)(2 3 4 5)"},
{"(0 1)(2 3 5)"},
{"(0 1)(2 4 5)"},
{"(0 1)(2 5)"},
{"(0 1)(2 4 3 5)"},
{"(0 1)(2 5)(3 4)"},
{"(0 2 1)"},
{"(0 2 1)(4 5)"},
{"(0 2 1)(3 4)"},
{"(0 2 1)(3 5 4)"},
{"(0 2 1)(3 4 5)"},
{"(0 2 1)(3 5)"},
{"(0 3 2 1)"},
{"(0 3 2 1)(4 5)"},
{"(0 4 3 2 1)"},
{"(0 5 4 3 2 1)"},
{"(0 4 5 3 2 1)"},
{"(0 5 3 2 1)"},
{"(0 3 4 2 1)"},
{"(0 3 5 4 2 1)"},
{"(0 4 2 1)"},
{"(0 5 4 2 1)"},
{"(0 4 2 1)(3 5)"},
{"(0 5 3 4 2 1)"},
{"(0 3 4 5 2 1)"},
{"(0 3 5 2 1)"},
{"(0 4 5 2 1)"},
{"(0 5 2 1)"},
{"(0 4 3 5 2 1)"},
{"(0 5 2 1)(3 4)"},
{"(0 2 3 1)"},
{"(0 2 3 1)(4 5)"},
{"(0 2 4 3 1)"},
{"(0 2 5 4 3 1)"},
{"(0 2 4 5 3 1)"},
{"(0 2 5 3 1)"},
{"(0 3 1)"},
{"(0 3 1)(4 5)"},
{"(0 4 3 1)"},
{"(0 5 4 3 1)"},
{"(0 4 5 3 1)"},
{"(0 5 3 1)"},
{"(0 3 1)(2 4)"},
{"(0 3 1)(2 5 4)"},
{"(0 4 2 3 1)"},
{"(0 5 4 2 3 1)"},
{"(0 4 2 5 3 1)"},
{"(0 5 3 1)(2 4)"},
{"(0 3 1)(2 4 5)"},
{"(0 3 1)(2 5)"},
{"(0 4 5 2 3 1)"},
{"(0 5 2 3 1)"},
{"(0 4 3 1)(2 5)"},
{"(0 5 2 4 3 1)"},
{"(0 2 3 4 1)"},
{"(0 2 3 5 4 1)"},
{"(0 2 4 1)"},
{"(0 2 5 4 1)"},
{"(0 2 4 1)(3 5)"},
{"(0 2 5 3 4 1)"},
{"(0 3 4 1)"},
{"(0 3 5 4 1)"},
{"(0 4 1)"},
{"(0 5 4 1)"},
{"(0 4 1)(3 5)"},
{"(0 5 3 4 1)"},
{"(0 3 2 4 1)"},
{"(0 3 2 5 4 1)"},
{"(0 4 1)(2 3)"},
{"(0 5 4 1)(2 3)"},
{"(0 4 1)(2 5 3)"},
{"(0 5 3 2 4 1)"},
{"(0 3 5 2 4 1)"},
{"(0 3 4 1)(2 5)"},
{"(0 4 1)(2 3 5)"},
{"(0 5 2 3 4 1)"},
{"(0 4 1)(2 5)"},
{"(0 5 2 4 1)"},
{"(0 2 3 4 5 1)"},
{"(0 2 3 5 1)"},
{"(0 2 4 5 1)"},
{"(0 2 5 1)"},
{"(0 2 4 3 5 1)"},
{"(0 2 5 1)(3 4)"},
{"(0 3 4 5 1)"},
{"(0 3 5 1)"},
{"(0 4 5 1)"},
{"(0 5 1)"},
{"(0 4 3 5 1)"},
{"(0 5 1)(3 4)"},
{"(0 3 2 4 5 1)"},
{"(0 3 2 5 1)"},
{"(0 4 5 1)(2 3)"},
{"(0 5 1)(2 3)"},
{"(0 4 3 2 5 1)"},
{"(0 5 1)(2 4 3)"},
{"(0 3 5 1)(2 4)"},
{"(0 3 4 2 5 1)"},
{"(0 4 2 3 5 1)"},
{"(0 5 1)(2 3 4)"},
{"(0 4 2 5 1)"},
{"(0 5 1)(2 4)"},
{"(0 1 2)"},
{"(0 1 2)(4 5)"},
{"(0 1 2)(3 4)"},
{"(0 1 2)(3 5 4)"},
{"(0 1 2)(3 4 5)"},
{"(0 1 2)(3 5)"},
{"(0 1 3 2)"},
{"(0 1 3 2)(4 5)"},
{"(0 1 4 3 2)"},
{"(0 1 5 4 3 2)"},
{"(0 1 4 5 3 2)"},
{"(0 1 5 3 2)"},
{"(0 1 3 4 2)"},
{"(0 1 3 5 4 2)"},
{"(0 1 4 2)"},
{"(0 1 5 4 2)"},
{"(0 1 4 2)(3 5)"},
{"(0 1 5 3 4 2)"},
{"(0 1 3 4 5 2)"},
{"(0 1 3 5 2)"},
{"(0 1 4 5 2)"},
{"(0 1 5 2)"},
{"(0 1 4 3 5 2)"},
{"(0 1 5 2)(3 4)"},
{"(0 2)"},
{"(0 2)(4 5)"},
{"(0 2)(3 4)"},
{"(0 2)(3 5 4)"},
{"(0 2)(3 4 5)"},
{"(0 2)(3 5)"},
{"(0 3 2)"},
{"(0 3 2)(4 5)"},
{"(0 4 3 2)"},
{"(0 5 4 3 2)"},
{"(0 4 5 3 2)"},
{"(0 5 3 2)"},
{"(0 3 4 2)"},
{"(0 3 5 4 2)"},
{"(0 4 2)"},
{"(0 5 4 2)"},
{"(0 4 2)(3 5)"},
{"(0 5 3 4 2)"},
{"(0 3 4 5 2)"},
{"(0 3 5 2)"},
{"(0 4 5 2)"},
{"(0 5 2)"},
{"(0 4 3 5 2)"},
{"(0 5 2)(3 4)"},
{"(0 2)(1 3)"},
{"(0 2)(1 3)(4 5)"},
{"(0 2)(1 4 3)"},
{"(0 2)(1 5 4 3)"},
{"(0 2)(1 4 5 3)"},
{"(0 2)(1 5 3)"},
{"(0 3 1 2)"},
{"(0 3 1 2)(4 5)"},
{"(0 4 3 1 2)"},
{"(0 5 4 3 1 2)"},
{"(0 4 5 3 1 2)"},
{"(0 5 3 1 2)"},
{"(0 3 1 4 2)"},
{"(0 3 1 5 4 2)"},
{"(0 4 2)(1 3)"},
{"(0 5 4 2)(1 3)"},
{"(0 4 2)(1 5 3)"},
{"(0 5 3 1 4 2)"},
{"(0 3 1 4 5 2)"},
{"(0 3 1 5 2)"},
{"(0 4 5 2)(1 3)"},
{"(0 5 2)(1 3)"},
{"(0 4 3 1 5 2)"},
{"(0 5 2)(1 4 3)"},
{"(0 2)(1 3 4)"},
{"(0 2)(1 3 5 4)"},
{"(0 2)(1 4)"},
{"(0 2)(1 5 4)"},
{"(0 2)(1 4)(3 5)"},
{"(0 2)(1 5 3 4)"},
{"(0 3 4 1 2)"},
{"(0 3 5 4 1 2)"},
{"(0 4 1 2)"},
{"(0 5 4 1 2)"},
{"(0 4 1 2)(3 5)"},
{"(0 5 3 4 1 2)"},
{"(0 3 2)(1 4)"},
{"(0 3 2)(1 5 4)"},
{"(0 4 1 3 2)"},
{"(0 5 4 1 3 2)"},
{"(0 4 1 5 3 2)"},
{"(0 5 3 2)(1 4)"},
{"(0 3 5 2)(1 4)"},
{"(0 3 4 1 5 2)"},
{"(0 4 1 3 5 2)"},
{"(0 5 2)(1 3 4)"},
{"(0 4 1 5 2)"},
{"(0 5 2)(1 4)"},
{"(0 2)(1 3 4 5)"},
{"(0 2)(1 3 5)"},
{"(0 2)(1 4 5)"},
{"(0 2)(1 5)"},
{"(0 2)(1 4 3 5)"},
{"(0 2)(1 5)(3 4)"},
{"(0 3 4 5 1 2)"},
{"(0 3 5 1 2)"},
{"(0 4 5 1 2)"},
{"(0 5 1 2)"},
{"(0 4 3 5 1 2)"},
{"(0 5 1 2)(3 4)"},
{"(0 3 2)(1 4 5)"},
{"(0 3 2)(1 5)"},
{"(0 4 5 1 3 2)"},
{"(0 5 1 3 2)"},
{"(0 4 3 2)(1 5)"},
{"(0 5 1 4 3 2)"},
{"(0 3 5 1 4 2)"},
{"(0 3 4 2)(1 5)"},
{"(0 4 2)(1 3 5)"},
{"(0 5 1 3 4 2)"},
{"(0 4 2)(1 5)"},
{"(0 5 1 4 2)"},
{"(0 1 2 3)"},
{"(0 1 2 3)(4 5)"},
{"(0 1 2 4 3)"},
{"(0 1 2 5 4 3)"},
{"(0 1 2 4 5 3)"},
{"(0 1 2 5 3)"},
{"(0 1 3)"},
{"(0 1 3)(4 5)"},
{"(0 1 4 3)"},
{"(0 1 5 4 3)"},
{"(0 1 4 5 3)"},
{"(0 1 5 3)"},
{"(0 1 3)(2 4)"},
{"(0 1 3)(2 5 4)"},
{"(0 1 4 2 3)"},
{"(0 1 5 4 2 3)"},
{"(0 1 4 2 5 3)"},
{"(0 1 5 3)(2 4)"},
{"(0 1 3)(2 4 5)"},
{"(0 1 3)(2 5)"},
{"(0 1 4 5 2 3)"},
{"(0 1 5 2 3)"},
{"(0 1 4 3)(2 5)"},
{"(0 1 5 2 4 3)"},
{"(0 2 3)"},
{"(0 2 3)(4 5)"},
{"(0 2 4 3)"},
{"(0 2 5 4 3)"},
{"(0 2 4 5 3)"},
{"(0 2 5 3)"},
{"(0 3)"},
{"(0 3)(4 5)"},
{"(0 4 3)"},
{"(0 5 4 3)"},
{"(0 4 5 3)"},
{"(0 5 3)"},
{"(0 3)(2 4)"},
{"(0 3)(2 5 4)"},
{"(0 4 2 3)"},
{"(0 5 4 2 3)"},
{"(0 4 2 5 3)"},
{"(0 5 3)(2 4)"},
{"(0 3)(2 4 5)"},
{"(0 3)(2 5)"},
{"(0 4 5 2 3)"},
{"(0 5 2 3)"},
{"(0 4 3)(2 5)"},
{"(0 5 2 4 3)"},
{"(0 2 1 3)"},
{"(0 2 1 3)(4 5)"},
{"(0 2 1 4 3)"},
{"(0 2 1 5 4 3)"},
{"(0 2 1 4 5 3)"},
{"(0 2 1 5 3)"},
{"(0 3)(1 2)"},
{"(0 3)(1 2)(4 5)"},
{"(0 4 3)(1 2)"},
{"(0 5 4 3)(1 2)"},
{"(0 4 5 3)(1 2)"},
{"(0 5 3)(1 2)"},
{"(0 3)(1 4 2)"},
{"(0 3)(1 5 4 2)"},
{"(0 4 2 1 3)"},
{"(0 5 4 2 1 3)"},
{"(0 4 2 1 5 3)"},
{"(0 5 3)(1 4 2)"},
{"(0 3)(1 4 5 2)"},
{"(0 3)(1 5 2)"},
{"(0 4 5 2 1 3)"},
{"(0 5 2 1 3)"},
{"(0 4 3)(1 5 2)"},
{"(0 5 2 1 4 3)"},
{"(0 2 4 1 3)"},
{"(0 2 5 4 1 3)"},
{"(0 2 3)(1 4)"},
{"(0 2 3)(1 5 4)"},
{"(0 2 5 3)(1 4)"},
{"(0 2 4 1 5 3)"},
{"(0 3)(1 2 4)"},
{"(0 3)(1 2 5 4)"},
{"(0 4 1 2 3)"},
{"(0 5 4 1 2 3)"},
{"(0 4 1 2 5 3)"},
{"(0 5 3)(1 2 4)"},
{"(0 3)(1 4)"},
{"(0 3)(1 5 4)"},
{"(0 4 1 3)"},
{"(0 5 4 1 3)"},
{"(0 4 1 5 3)"},
{"(0 5 3)(1 4)"},
{"(0 3)(1 4)(2 5)"},
{"(0 3)(1 5 2 4)"},
{"(0 4 1 3)(2 5)"},
{"(0 5 2 4 1 3)"},
{"(0 4 1 5 2 3)"},
{"(0 5 2 3)(1 4)"},
{"(0 2 4 5 1 3)"},
{"(0 2 5 1 3)"},
{"(0 2 3)(1 4 5)"},
{"(0 2 3)(1 5)"},
{"(0 2 5 1 4 3)"},
{"(0 2 4 3)(1 5)"},
{"(0 3)(1 2 4 5)"},
{"(0 3)(1 2 5)"},
{"(0 4 5 1 2 3)"},
{"(0 5 1 2 3)"},
{"(0 4 3)(1 2 5)"},
{"(0 5 1 2 4 3)"},
{"(0 3)(1 4 5)"},
{"(0 3)(1 5)"},
{"(0 4 5 1 3)"},
{"(0 5 1 3)"},
{"(0 4 3)(1 5)"},
{"(0 5 1 4 3)"},
{"(0 3)(1 4 2 5)"},
{"(0 3)(1 5)(2 4)"},
{"(0 4 2 5 1 3)"},
{"(0 5 1 3)(2 4)"},
{"(0 4 2 3)(1 5)"},
{"(0 5 1 4 2 3)"},
{"(0 1 2 3 4)"},
{"(0 1 2 3 5 4)"},
{"(0 1 2 4)"},
{"(0 1 2 5 4)"},
{"(0 1 2 4)(3 5)"},
{"(0 1 2 5 3 4)"},
{"(0 1 3 4)"},
{"(0 1 3 5 4)"},
{"(0 1 4)"},
{"(0 1 5 4)"},
{"(0 1 4)(3 5)"},
{"(0 1 5 3 4)"},
{"(0 1 3 2 4)"},
{"(0 1 3 2 5 4)"},
{"(0 1 4)(2 3)"},
{"(0 1 5 4)(2 3)"},
{"(0 1 4)(2 5 3)"},
{"(0 1 5 3 2 4)"},
{"(0 1 3 5 2 4)"},
{"(0 1 3 4)(2 5)"},
{"(0 1 4)(2 3 5)"},
{"(0 1 5 2 3 4)"},
{"(0 1 4)(2 5)"},
{"(0 1 5 2 4)"},
{"(0 2 3 4)"},
{"(0 2 3 5 4)"},
{"(0 2 4)"},
{"(0 2 5 4)"},
{"(0 2 4)(3 5)"},
{"(0 2 5 3 4)"},
{"(0 3 4)"},
{"(0 3 5 4)"},
{"(0 4)"},
{"(0 5 4)"},
{"(0 4)(3 5)"},
{"(0 5 3 4)"},
{"(0 3 2 4)"},
{"(0 3 2 5 4)"},
{"(0 4)(2 3)"},
{"(0 5 4)(2 3)"},
{"(0 4)(2 5 3)"},
{"(0 5 3 2 4)"},
{"(0 3 5 2 4)"},
{"(0 3 4)(2 5)"},
{"(0 4)(2 3 5)"},
{"(0 5 2 3 4)"},
{"(0 4)(2 5)"},
{"(0 5 2 4)"},
{"(0 2 1 3 4)"},
{"(0 2 1 3 5 4)"},
{"(0 2 1 4)"},
{"(0 2 1 5 4)"},
{"(0 2 1 4)(3 5)"},
{"(0 2 1 5 3 4)"},
{"(0 3 4)(1 2)"},
{"(0 3 5 4)(1 2)"},
{"(0 4)(1 2)"},
{"(0 5 4)(1 2)"},
{"(0 4)(1 2)(3 5)"},
{"(0 5 3 4)(1 2)"},
{"(0 3 2 1 4)"},
{"(0 3 2 1 5 4)"},
{"(0 4)(1 3 2)"},
{"(0 5 4)(1 3 2)"},
{"(0 4)(1 5 3 2)"},
{"(0 5 3 2 1 4)"},
{"(0 3 5 2 1 4)"},
{"(0 3 4)(1 5 2)"},
{"(0 4)(1 3 5 2)"},
{"(0 5 2 1 3 4)"},
{"(0 4)(1 5 2)"},
{"(0 5 2 1 4)"},
{"(0 2 4)(1 3)"},
{"(0 2 5 4)(1 3)"},
{"(0 2 3 1 4)"},
{"(0 2 3 1 5 4)"},
{"(0 2 5 3 1 4)"},
{"(0 2 4)(1 5 3)"},
{"(0 3 1 2 4)"},
{"(0 3 1 2 5 4)"},
{"(0 4)(1 2 3)"},
{"(0 5 4)(1 2 3)"},
{"(0 4)(1 2 5 3)"},
{"(0 5 3 1 2 4)"},
{"(0 3 1 4)"},
{"(0 3 1 5 4)"},
{"(0 4)(1 3)"},
{"(0 5 4)(1 3)"},
{"(0 4)(1 5 3)"},
{"(0 5 3 1 4)"},
{"(0 3 1 4)(2 5)"},
{"(0 3 1 5 2 4)"},
{"(0 4)(1 3)(2 5)"},
{"(0 5 2 4)(1 3)"},
{"(0 4)(1 5 2 3)"},
{"(0 5 2 3 1 4)"},
{"(0 2 4)(1 3 5)"},
{"(0 2 5 1 3 4)"},
{"(0 2 3 5 1 4)"},
{"(0 2 3 4)(1 5)"},
{"(0 2 5 1 4)"},
{"(0 2 4)(1 5)"},
{"(0 3 5 1 2 4)"},
{"(0 3 4)(1 2 5)"},
{"(0 4)(1 2 3 5)"},
{"(0 5 1 2 3 4)"},
{"(0 4)(1 2 5)"},
{"(0 5 1 2 4)"},
{"(0 3 5 1 4)"},
{"(0 3 4)(1 5)"},
{"(0 4)(1 3 5)"},
{"(0 5 1 3 4)"},
{"(0 4)(1 5)"},
{"(0 5 1 4)"},
{"(0 3 2 5 1 4)"},
{"(0 3 2 4)(1 5)"},
{"(0 4)(1 3 2 5)"},
{"(0 5 1 3 2 4)"},
{"(0 4)(1 5)(2 3)"},
{"(0 5 1 4)(2 3)"},
{"(0 1 2 3 4 5)"},
{"(0 1 2 3 5)"},
{"(0 1 2 4 5)"},
{"(0 1 2 5)"},
{"(0 1 2 4 3 5)"},
{"(0 1 2 5)(3 4)"},
{"(0 1 3 4 5)"},
{"(0 1 3 5)"},
{"(0 1 4 5)"},
{"(0 1 5)"},
{"(0 1 4 3 5)"},
{"(0 1 5)(3 4)"},
{"(0 1 3 2 4 5)"},
{"(0 1 3 2 5)"},
{"(0 1 4 5)(2 3)"},
{"(0 1 5)(2 3)"},
{"(0 1 4 3 2 5)"},
{"(0 1 5)(2 4 3)"},
{"(0 1 3 5)(2 4)"},
{"(0 1 3 4 2 5)"},
{"(0 1 4 2 3 5)"},
{"(0 1 5)(2 3 4)"},
{"(0 1 4 2 5)"},
{"(0 1 5)(2 4)"},
{"(0 2 3 4 5)"},
{"(0 2 3 5)"},
{"(0 2 4 5)"},
{"(0 2 5)"},
{"(0 2 4 3 5)"},
{"(0 2 5)(3 4)"},
{"(0 3 4 5)"},
{"(0 3 5)"},
{"(0 4 5)"},
{"(0 5)"},
{"(0 4 3 5)"},
{"(0 5)(3 4)"},
{"(0 3 2 4 5)"},
{"(0 3 2 5)"},
{"(0 4 5)(2 3)"},
{"(0 5)(2 3)"},
{"(0 4 3 2 5)"},
{"(0 5)(2 4 3)"},
{"(0 3 5)(2 4)"},
{"(0 3 4 2 5)"},
{"(0 4 2 3 5)"},
{"(0 5)(2 3 4)"},
{"(0 4 2 5)"},
{"(0 5)(2 4)"},
{"(0 2 1 3 4 5)"},
{"(0 2 1 3 5)"},
{"(0 2 1 4 5)"},
{"(0 2 1 5)"},
{"(0 2 1 4 3 5)"},
{"(0 2 1 5)(3 4)"},
{"(0 3 4 5)(1 2)"},
{"(0 3 5)(1 2)"},
{"(0 4 5)(1 2)"},
{"(0 5)(1 2)"},
{"(0 4 3 5)(1 2)"},
{"(0 5)(1 2)(3 4)"},
{"(0 3 2 1 4 5)"},
{"(0 3 2 1 5)"},
{"(0 4 5)(1 3 2)"},
{"(0 5)(1 3 2)"},
{"(0 4 3 2 1 5)"},
{"(0 5)(1 4 3 2)"},
{"(0 3 5)(1 4 2)"},
{"(0 3 4 2 1 5)"},
{"(0 4 2 1 3 5)"},
{"(0 5)(1 3 4 2)"},
{"(0 4 2 1 5)"},
{"(0 5)(1 4 2)"},
{"(0 2 4 5)(1 3)"},
{"(0 2 5)(1 3)"},
{"(0 2 3 1 4 5)"},
{"(0 2 3 1 5)"},
{"(0 2 5)(1 4 3)"},
{"(0 2 4 3 1 5)"},
{"(0 3 1 2 4 5)"},
{"(0 3 1 2 5)"},
{"(0 4 5)(1 2 3)"},
{"(0 5)(1 2 3)"},
{"(0 4 3 1 2 5)"},
{"(0 5)(1 2 4 3)"},
{"(0 3 1 4 5)"},
{"(0 3 1 5)"},
{"(0 4 5)(1 3)"},
{"(0 5)(1 3)"},
{"(0 4 3 1 5)"},
{"(0 5)(1 4 3)"},
{"(0 3 1 4 2 5)"},
{"(0 3 1 5)(2 4)"},
{"(0 4 2 5)(1 3)"},
{"(0 5)(1 3)(2 4)"},
{"(0 4 2 3 1 5)"},
{"(0 5)(1 4 2 3)"},
{"(0 2 4 1 3 5)"},
{"(0 2 5)(1 3 4)"},
{"(0 2 3 5)(1 4)"},
{"(0 2 3 4 1 5)"},
{"(0 2 5)(1 4)"},
{"(0 2 4 1 5)"},
{"(0 3 5)(1 2 4)"},
{"(0 3 4 1 2 5)"},
{"(0 4 1 2 3 5)"},
{"(0 5)(1 2 3 4)"},
{"(0 4 1 2 5)"},
{"(0 5)(1 2 4)"},
{"(0 3 5)(1 4)"},
{"(0 3 4 1 5)"},
{"(0 4 1 3 5)"},
{"(0 5)(1 3 4)"},
{"(0 4 1 5)"},
{"(0 5)(1 4)"},
{"(0 3 2 5)(1 4)"},
{"(0 3 2 4 1 5)"},
{"(0 4 1 3 2 5)"},
{"(0 5)(1 3 2 4)"},
{"(0 4 1 5)(2 3)"},
{"(0 5)(1 4)(2 3)"}
};

pd_perm_pp_t* pd_pc_perm_print[7] = {NULL, pdint_pcperm_print1, pdint_pcperm_print2, pdint_pcperm_print3, pdint_pcperm_print4, pdint_pcperm_print5, pdint_pcperm_print6};

#endif
