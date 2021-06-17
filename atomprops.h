/* name: atomprops.h                                     */
/* author: J. Michael Word     date written: 6/12/97     */
/* purpose: define atom properties                       */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef ATOMPROPS_H
#define ATOMPROPS_H 1

#define HE_RESNAMES \
 "currently there are none RMI 070711"

#define HF_RESNAMES \
 ":PHF:HF3:HF5:"

#define HG_RESNAMES \
 ": HG:HG2:HGB:HGC:HGI:MAC:MBO:MMC:PHG:PMB:AAS:AMS:BE7:CMH:EMC:EMT:"

#define HO_RESNAMES \
 ": HO:HO3:"

#define HS_RESNAMES \
 "currently there are none RMI 070711"

typedef struct atomProp_t {
   int    type;     /* atom identifier                 */
   int    atno;     /* atomic number                   */
   char*  name;     /* atom name                       */
   float  eRad;     /* (explicit H) ecloud VDW radius  */
   float  eRad_nuc; /* (explicit H) nuclear VDW radius */
   float  iRad;     /* (implicit H) VDW radius         */
   float  covRad;   /* covalent radius                 */
   char* color;     /* dot color                       */
   int    flags;    /* element features                */
} atomProp;

void initalizeAtomTbl(void);
int fixAtomName(const char* atomname, char resname[], int position); /*rmi070719 add resname info to identifyAtom*/
int identifyAtom(char* name, char resname[], int Verbose); /*dcr041007 allow warning choice*/

int   getAtno(int a);
char* getAtomName(int a);
float getExplRad(int a);
float getImplRad(int a);
float getCovRad(int a);
char* getColor(int a);
float getMaxRadius(int isImplicit);
int   atomHasProp(int a, int f);

/* Each atom identifier (except endAtomTypes) must be in the AtomTbl      */
/* A new atom must also be added to the function identifyAtomName() and   */
/* may need to be entered in setProperties()                              */
enum atomIdentifiers {
   noAtom=0, ignoreAtom,
   atomH, atomHarom, atomHpolar, atomHOd,
   atomC, atomCarom, atomN, atomO, atomF, atomS, atomP,

   atomAc, atomAg, atomAl, atomAm, atomAr, atomAs, atomAt, atomAu,
   atomB, atomBa, atomBe, atomBi, atomBk, atomBr, atomCa, atomCd,
   atomCe, atomCf, atomCl, atomCm, atomCo, atomCr, atomCs, atomCu,
   atomDy, atomEr, atomEs, atomEu, atomFe, atomFm, atomFr, atomGa,
   atomGd, atomGe, atomHe, atomHf, atomHg, atomHo, atomI, atomIn,
   atomIr, atomK, atomKr, atomLa, atomLi, atomLu, atomMd, atomMg,
   atomMn, atomMo, atomNa, atomNb, atomNd, atomNe, atomNi, atomNo,
   atomNp, atomOs, atomPa, atomPb, atomPd, atomPm, atomPo, atomPr,
   atomPt, atomPu, atomRa, atomRb, atomRe, atomRh, atomRn, atomRu,
   atomSb, atomSc, atomSe, atomSi, atomSm, atomSn, atomSr, atomTa,
   atomTb, atomTc, atomTe, atomTh, atomTi, atomTl, atomTm, atomU,
   atomV, atomW, atomXe, atomY, atomYb, atomZn, atomZr,

   /* base identifiers are used in a real hack to paint nucleic acid  */
   /* contacts with colors for bases as an alternative to atom colors */
   /* (dots are catagorized in a table NUMATOMTYPES wide)             */
   baseA, baseTU, baseC, baseG, baseOther, nonBase,

   endAtomTypes
};

#define isHatom(atype) (getAtno(atype) == 1)
#define isCatom(atype) (getAtno(atype) == 6)
#define NUMATOMTYPES endAtomTypes

#define COVRADFUDGE 0.2

#define METALIC_ATOM_FLAG (1 <<  0) /*atomProp AtomTbl flags: element features*/
#define   IONIC_ATOM_FLAG (1 <<  1) /*dcr041007 partial, now on halides*/

#ifdef INIT_ATOM_TABLE
/* For non-metals, explicit VDW radii from                      */
/* Gavezzotti, J. Am. Chem. Soc. (1983) 105, 5220-5225.         */
/* or, if unavailable,                                          */
/* Bondi, J. Phys. Chem. (1964), V68, N3, 441-451.              */
/* Covalent and ionic radii from                                */
/* Advanced Inorganic Chemistry, Cotton & Wilkinson, 1962, p93. */

/* from above: struct atomProp, NUMATOMTYPES from enum of atomIdentifiers*/
   /*type;   atom identifier                   */
   /*atno;   atomic number                     */
   /*name;   atom name                         */
   /*eRad;   (explicit H) ecloud VDW radius    */
   /*eRad_nuc; (explicit H) nuclear VDW radius */ //jjh 130307
   /*iRad;   (implicit H) VDW radius           */
   /*covRad; covalent radius                   */
   /*color;  dot color                         */
   /*flags;  element features                  */

atomProp AtomTbl[NUMATOMTYPES] = { /* noAtom must be first */
   {noAtom,      0, "?",    1.05, 1.05, 0.00, 0.00, "magenta", 0},
   {ignoreAtom,  0, ".",    0.00, 0.00, 0.00, 0.00, "grey",   0},

   {atomH,       1, "H",    1.22, 1.17, 0.00, 0.30, "grey",   0},
   {atomHarom,   1, "Har",  1.05, 1.00, 0.00, 0.30, "grey",   0},
   {atomHpolar,  1, "Hpol", 1.05, 1.00, 0.00, 0.30, "grey",   0},
   {atomHOd,     1, "HOd",  1.05, 1.00, 0.00, 0.30, "grey",   0},/*hb-only-dummy*/
   {atomC,       6, "C",    1.70, 1.70, 1.90, 0.77, "white",  0},
   {atomCarom,   6, "Car",  1.75, 1.75, 1.90, 0.77, "white",  0},
   {atomN,       7, "N",    1.55, 1.55, 1.70, 0.70, "sky",    0},
   {atomO,       8, "O",    1.40, 1.40, 1.50, 0.66, "red",    0},
   {atomP,      15, "P",    1.80, 1.80, 1.80, 1.10, "pink",   0},
   {atomS,      16, "S",    1.80, 1.80, 1.90, 1.04, "yellow", 0},
   {atomAs,     33, "As",   2.00, 2.00, 2.10, 1.21, "grey",   0},
   {atomSe,     34, "Se",   1.90, 1.90, 2.00, 1.17, "green",  0},
/*dcr041007 halides given IONIC_ATOM_FLAG*/
   {atomF,       9, "F",    1.30, 1.30, 1.30, 0.58, "green", IONIC_ATOM_FLAG},
   {atomCl,     17, "Cl",   1.77, 1.77, 1.77, 0.99, "green", IONIC_ATOM_FLAG},
   {atomBr,     35, "Br",   1.95, 1.95, 1.95, 1.14, "brown", IONIC_ATOM_FLAG},
   {atomI,      53, "I",    2.10, 2.10, 2.10, 1.33, "brown", IONIC_ATOM_FLAG},

 /* for most common metals we use Pauling's ionic radii                       */
 /* "covalent radii" = ionic + 0.74 (i.e., oxygenVDW(1.4) - oxygenCov(0.66))  */
 /* because the ionic radii are usually calculated from Oxygen-Metal distance */

   {atomLi,      3, "Li", 1.34, 1.34, 1.34,   0.60, "grey", METALIC_ATOM_FLAG},
   {atomNa,     11, "Na", 1.69, 1.69, 1.69,   0.95, "grey", METALIC_ATOM_FLAG},
   {atomAl,     13, "Al", 1.24, 1.24, 1.24,   0.50, "grey", METALIC_ATOM_FLAG},
   {atomK,      19, "K", 2.07, 2.07, 2.07,    1.33, "grey", METALIC_ATOM_FLAG},
   {atomMg,     12, "Mg", 1.39, 1.39, 1.39,   0.65, "grey", METALIC_ATOM_FLAG},
   {atomCa,     20, "Ca", 1.73, 1.73, 1.73,   0.99, "grey", METALIC_ATOM_FLAG},
   {atomMn,     25, "Mn", 1.54, 1.54, 1.54,   0.80, "grey", METALIC_ATOM_FLAG},
   {atomFe,     26, "Fe", 1.48, 1.48, 1.48,   0.74, "grey", METALIC_ATOM_FLAG},
   {atomCo,     27, "Co", 1.44, 1.44, 1.44,   0.70, "blue", METALIC_ATOM_FLAG},
   {atomNi,     28, "Ni", 1.40, 1.40, 1.40,   0.66, "grey", METALIC_ATOM_FLAG},
   {atomCu,     29, "Cu", 1.46, 1.46, 1.46,   0.72,"orange",METALIC_ATOM_FLAG},
   {atomZn,     30, "Zn", 1.45, 1.45, 1.45,   0.71, "grey", METALIC_ATOM_FLAG},
   {atomRb,     37, "Rb", 2.22, 2.22, 2.22,   1.48, "grey", METALIC_ATOM_FLAG},
   {atomSr,     38, "Sr", 1.84, 1.84, 1.84,   1.10, "grey", METALIC_ATOM_FLAG},
   {atomMo,     42, "Mo", 1.67, 1.67, 1.67,   0.93, "grey", METALIC_ATOM_FLAG},
   {atomAg,     47, "Ag", 2.00, 2.00, 2.00,   1.26, "white",METALIC_ATOM_FLAG},
   {atomCd,     48, "Cd", 1.65, 1.65, 1.65,   0.91, "grey", METALIC_ATOM_FLAG},
   {atomIn,     49, "In", 1.55, 1.55, 1.55,   0.81, "grey", METALIC_ATOM_FLAG},
   {atomCs,     55, "Cs", 2.43, 2.43, 2.43,   1.69, "grey", METALIC_ATOM_FLAG},
   {atomBa,     56, "Ba", 2.03, 2.03, 2.03,   1.29, "grey", METALIC_ATOM_FLAG},
   {atomAu,     79, "Au", 1.84, 1.84, 1.84,   1.10, "gold", METALIC_ATOM_FLAG},
   {atomHg,     80, "Hg", 1.74, 1.74, 1.74,   1.00, "grey", METALIC_ATOM_FLAG},
   {atomTl,     81, "Tl", 2.18, 2.18, 2.18,   1.44, "grey", METALIC_ATOM_FLAG},
   {atomPb,     82, "Pb", 1.58, 1.58, 1.58,   0.84, "grey", METALIC_ATOM_FLAG},

/* for other metals we use Shannon's ionic radii */
/* Acta Crystallogr. (1975) A32, pg751.          */
   {atomV,      23, "V", 1.53, 1.53, 1.53,    0.79, "grey", METALIC_ATOM_FLAG},
   {atomCr,     24, "Cr", 1.47, 1.47, 1.47,   0.73, "grey", METALIC_ATOM_FLAG},
   {atomTe,     52, "Te", 1.71, 1.71, 1.71,   0.97, "grey", METALIC_ATOM_FLAG},
   {atomSm,     62, "Sm", 1.82, 1.82, 1.82,   1.08, "grey", METALIC_ATOM_FLAG},
   {atomGd,     64, "Gd", 1.79, 1.79, 1.79,   1.05, "grey", METALIC_ATOM_FLAG},
   {atomYb,     70, "Yb", 1.88, 1.88, 1.88,   1.14, "grey", METALIC_ATOM_FLAG},
   {atomW,      74, "W", 1.40, 1.40, 1.40,    0.66, "grey", METALIC_ATOM_FLAG},
   {atomPt,     78, "Pt", 1.37, 1.37, 1.37,   0.63, "grey", METALIC_ATOM_FLAG},
   {atomU,      92, "U", 1.77, 1.77, 1.77,    1.03, "grey", METALIC_ATOM_FLAG},

/* Cotton & Wilkinson and also- */
/* L.E. Sutton (ed.) in Table of interatomic distances and configuration in molecules */
/* and ions, Supplement 1956-1959, Special publication No. 18, Chemical Society,      */
/* London, UK, 1965 (as listed in web-elements by Mark Winter)                        */
/*                   http://www.shef.ac.uk/chemistry/web-elements                     */

   {atomHe,      2, "He",    1.60, 1.60, 1.60, 0.00, "sky",             0},
   {atomBe,      4, "Be", 0.90, 0.90, 0.90,    0.31, "grey", METALIC_ATOM_FLAG},
   {atomB,       5, "B",     0.20, 0.20, 0.20, 0.86, "grey",            0},
   {atomNe,     10, "Ne",    1.60, 1.60, 1.60, 0.00, "pink",            0},
   {atomSi,     14, "Si", 1.17, 1.17, 1.17,    2.10, "grey", METALIC_ATOM_FLAG},
   {atomAr,     18, "Ar",    1.89, 1.89, 1.89, 0.00, "orange",          0},
   {atomSc,     21, "Sc", 0.44, 0.44, 0.44,    0.68, "grey", METALIC_ATOM_FLAG},
   {atomTi,     22, "Ti", 1.49, 1.49, 1.49,    0.75, "grey", METALIC_ATOM_FLAG},
   {atomGa,     31, "Ga", 1.27, 1.27, 1.27,    0.53, "grey", METALIC_ATOM_FLAG},
   {atomGe,     32, "Ge", 1.34, 1.34, 1.34,    0.60, "grey", METALIC_ATOM_FLAG},
   {atomKr,     36, "Kr",    2.01, 2.01, 2.01, 1.15, "greentint",       0},
   {atomY,      39, "Y", 1.64, 1.64, 1.64,     0.90, "grey", METALIC_ATOM_FLAG},
   {atomZr,     40, "Zr", 1.51, 1.51, 1.51,    0.77, "grey", METALIC_ATOM_FLAG},
   {atomSn,     50, "Sn", 1.45, 1.45, 1.45,    0.71, "grey", METALIC_ATOM_FLAG},
   {atomSb,     51, "Sb", 1.41, 1.41, 1.41,    2.20, "grey", METALIC_ATOM_FLAG},
   {atomXe,     54, "Xe",    2.18, 2.18, 2.18, 1.28, "magenta",         0},
   {atomLa,     57, "La", 1.77, 1.77, 1.77,    1.03, "grey", METALIC_ATOM_FLAG},
   {atomCe,     58, "Ce", 1.61, 1.61, 1.61,    0.87, "grey", METALIC_ATOM_FLAG},
   {atomFr,     87, "Fr", 2.68, 2.68, 2.68,    1.94, "grey", METALIC_ATOM_FLAG},
   {atomRa,     88, "Ra", 2.36, 2.36, 2.36,    1.62, "grey", METALIC_ATOM_FLAG},
   {atomTh,     90, "Th", 1.82, 1.82, 1.82,    1.08, "grey", METALIC_ATOM_FLAG},

/* finally, we have a set of elements where the radii are unknown    */
/* so we use estimates and extrapolations based on web-elements data */
   {atomNb,     41, "Nb", 1.40, 1.40, 1.40,    0.86, "grey", METALIC_ATOM_FLAG},
   {atomTc,     43, "Tc", 1.25, 1.25, 1.25,    0.71, "grey", METALIC_ATOM_FLAG},
   {atomRu,     44, "Ru", 1.36, 1.36, 1.36,    0.82, "grey", METALIC_ATOM_FLAG},
   {atomRh,     45, "Rh", 1.30, 1.30, 1.30,    0.76, "grey", METALIC_ATOM_FLAG},
   {atomPd,     46, "Pd", 1.59, 1.59, 1.59,    1.05, "grey", METALIC_ATOM_FLAG},
   {atomPr,     59, "Pr", 1.65, 1.65, 1.65,    1.11, "grey", METALIC_ATOM_FLAG},
   {atomNd,     60, "Nd", 1.64, 1.64, 1.64,    1.10, "grey", METALIC_ATOM_FLAG},
   {atomPm,     61, "Pm", 1.89, 1.89, 1.89,    1.15, "grey", METALIC_ATOM_FLAG},
   {atomEu,     63, "Eu", 1.85, 1.85, 1.85,    1.31, "grey", METALIC_ATOM_FLAG},
   {atomTb,     65, "Tb", 1.59, 1.59, 1.59,    1.05, "grey", METALIC_ATOM_FLAG},
   {atomDy,     66, "Dy", 1.59, 1.59, 1.59,    1.05, "grey", METALIC_ATOM_FLAG},
   {atomHo,     67, "Ho", 1.58, 1.58, 1.58,    1.04, "grey", METALIC_ATOM_FLAG},
   {atomEr,     68, "Er", 1.57, 1.57, 1.57,    1.03, "grey", METALIC_ATOM_FLAG},
   {atomTm,     69, "Tm", 1.56, 1.56, 1.56,    1.02, "grey", METALIC_ATOM_FLAG},
   {atomLu,     71, "Lu", 1.56, 1.56, 1.56,    1.02, "grey", METALIC_ATOM_FLAG},
   {atomHf,     72, "Hf", 1.46, 1.46, 1.46,    0.85, "grey", METALIC_ATOM_FLAG},
   {atomTa,     73, "Ta", 1.40, 1.40, 1.40,    0.86, "grey", METALIC_ATOM_FLAG},
   {atomRe,     75, "Re", 1.31, 1.31, 1.31,    0.77, "grey", METALIC_ATOM_FLAG},
   {atomOs,     76, "Os", 1.32, 1.32, 1.32,    0.78, "grey", METALIC_ATOM_FLAG},
   {atomIr,     77, "Ir", 1.34, 1.34, 1.34,    0.80, "grey", METALIC_ATOM_FLAG},
   {atomBi,     83, "Bi", 1.71, 1.71, 1.71,    1.17, "grey", METALIC_ATOM_FLAG},
   {atomPo,     84, "Po", 1.53, 1.53, 1.53,    0.99, "grey", METALIC_ATOM_FLAG},
   {atomAt,     85, "At", 1.45, 1.45, 1.45,    0.91, "grey", METALIC_ATOM_FLAG},
   {atomRn,     86, "Rn",    2.50, 2.50, 2.50, 1.25, "pinktint", 0},
   {atomAc,     89, "Ac", 2.00, 2.00, 2.00,    1.30, "grey", METALIC_ATOM_FLAG},
   {atomPa,     91, "Pa", 1.85, 1.85, 1.85,    1.10, "grey", METALIC_ATOM_FLAG},
   {atomNp,     93, "Np", 1.72, 1.72, 1.72,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomPu,     94, "Pu", 1.67, 1.67, 1.67,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomAm,     95, "Am", 1.63, 1.63, 1.63,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomCm,     96, "Cm", 1.60, 1.60, 1.60,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomBk,     97, "Bk", 1.58, 1.58, 1.58,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomCf,     98, "Cf", 1.57, 1.57, 1.57,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomEs,     99, "Es", 1.56, 1.56, 1.56,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomFm,     100,"Fm", 1.55, 1.55, 1.55,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomMd,     101,"Md", 1.55, 1.55, 1.55,    1.00, "grey", METALIC_ATOM_FLAG},
   {atomNo,     102,"No", 1.55, 1.55, 1.55,    1.00, "grey", METALIC_ATOM_FLAG},

   /* base identifiers are used in color hack (see note in enum above) */
   {baseA,      0, "a",      0.00, 0.00, 0.00, 0.00, "pink",         0},
   {baseC,      0, "c",      0.00, 0.00, 0.00, 0.00, "yellow",       0},
   {baseTU,     0, "t/u",    0.00, 0.00, 0.00, 0.00, "sky",          0},
   {baseG,      0, "g",      0.00, 0.00, 0.00, 0.00, "sea",          0},
   {baseOther,  0,"other na",0.00, 0.00, 0.00, 0.00, "white",        0},
   {nonBase,    0, "nonbase",0.00, 0.00, 0.00, 0.00, "grey",         0}
};
#endif

#endif
