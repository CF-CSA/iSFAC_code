/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   defines.h
 * Author: tg
 *
 * Created on October 11, 2024, 1:59 PM
 */

#ifndef DEFINES_H
#define DEFINES_H

namespace Physics {
    const double a0 = 0.5291772; // Bohr radius in A A/Bohr
}

namespace PSE {
    /**
    Van der Waals radii in [A] taken from:
    A cartography of the van der Waals territories
    S. Alvarez, Dalton Trans., 2013, 42, 8617-8636
    DOI: 10.1039/C3DT50599E
    */
    const double vdw_radii[] = {
        -1.0, // X
        1.2, // H
        1.43, // He [larger uncertainty]
        2.12, // Li
        1.98, // Be
        1.91, // B
        1.77, // C
        1.66, // N
        1.5, // O
        1.46, // F
        1.58, // Ne [larger uncertainty]
        2.5, // Na
        2.51, // Mg
        2.25, // Al
        2.19, // Si
        1.9, // P
        1.89, // S
        1.82, // Cl
        1.83, // Ar
        2.73, // K
        2.62, // Ca
        2.58, // Sc
        2.46, // Ti
        2.42, // V
        2.45, // Cr
        2.45, // Mn
        2.44, // Fe
        2.4, // Co
        2.4, // Ni
        2.38, // Cu
        2.39, // Zn
        2.32, // Ga
        2.29, // Ge
        1.88, // As
        1.82, // Se
        1.86, // Br
        2.25, // Kr
        3.21, // Rb
        2.84, // Sr
        2.75, // Y
        2.52, // Zr
        2.56, // Nb
        2.45, // Mo
        2.44, // Tc
        2.46, // Ru
        2.44, // Rh
        2.15, // Pd
        2.53, // Ag
        2.49, // Cd
        2.43, // In
        2.42, // Sn
        2.47, // Sb
        1.99, // Te
        2.04, // I
        2.06, // Xe
        3.48, // Cs
        3.03, // Ba
        2.98, // La
        2.88, // Ce
        2.92, // Pr
        2.95, // Nd
        -1.00, // Pm
        2.9, // Sm
        2.87, // Eu
        2.83, // Gd
        2.79, // Tb
        2.87, // Dy
        2.81, // Ho
        2.83, // Er
        2.79, // Tm
        2.8, // Yb
        2.74, // Lu
        2.63, // Hf
        2.53, // Ta
        2.57, // W
        2.49, // Re
        2.48, // Os
        2.41, // Ir
        2.29, // Pt
        2.32, // Au
        2.45, // Hg
        2.47, // Tl
        2.6, // Pb
        2.54, // Bi
        -1.00, // Po
        -1.00, // At
        -1.00, // Rn
        -1.00, // Fr
        -1.00, // Ra
        2.8, // Ac [larger uncertainty]
        2.93, // Th
        2.88, // Pa [larger uncertainty]
        2.71, // U
        2.82, // Np
        2.81, // Pu
        2.83, // Am
        3.05, // Cm [larger uncertainty]
        3.4, // Bk [larger uncertainty]
        3.05, // Cf [larger uncertainty]
        2.7, // Es [larger uncertainty]
        -1.00, // Fm
        -1.00, // Md
        -1.00, // No
        -1.00, // Lr
    };


    const char* const Elements[] = {
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
        "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti",
        "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
        "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
        "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
        "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
        "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Uun", "Uuu", "Uub"
    };
    
    /**
        from  ase.data.atomic_mases
     *  print (f'', *atomic_masses, sep=", ")
     */
    const double AtomicWeight[] = { -1.0, 1.008, 4.002602, 6.94, 9.0121831, 10.81, 
        12.011, 14.007, 15.999, 18.998403163, 20.1797, 22.98976928, 24.305, 
        26.9815385, 28.085, 30.973761998, 32.06, 35.45, 39.948, 39.0983, 
        40.078, 44.955908, 47.867, 50.9415, 51.9961, 54.938044, 55.845, 
        58.933194, 58.6934, 63.546, 65.38, 69.723, 72.63, 74.921595, 78.971, 
        79.904, 83.798, 85.4678, 87.62, 88.90584, 91.224, 92.90637, 95.95, 
        97.90721, 101.07, 102.9055, 106.42, 107.8682, 112.414, 114.818, 
        118.71, 121.76, 127.6, 126.90447, 131.293, 132.90545196, 137.327, 
        138.90547, 140.116, 140.90766, 144.242, 144.91276, 150.36, 151.964, 
        157.25, 158.92535, 162.5, 164.93033, 167.259, 168.93422, 173.054, 
        174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 
        195.084, 196.966569, 200.592, 204.38, 207.2, 208.9804, 208.98243, 
        209.98715, 222.01758, 223.01974, 226.02541, 227.02775, 232.0377, 
        231.03588, 238.02891, 237.04817, 244.06421, 243.06138, 247.07035, 
        247.07031, 251.07959, 252.083, 257.09511, 258.09843, 259.101, 262.11, 
        267.122, 268.126, 271.134, 270.133, 269.1338, 278.156, 281.165, 
        281.166, 285.177, 286.182, 289.19, 289.194, 293.204, 293.208, 294.214        
    };

    // H=1 so that a value of -1 denotes an invalid type
    // NONE is used for empty atoms (e.g. to mark end of a list,
    // ZERO is used for dummy atoms with no scattering power
    //

    enum Ord_Num {
        NONE = -2, DUMMY = -1, ZERO = 0,
        H = 1, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl,
        Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br,
        Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,
        Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb,
        Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr,
        Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf,
        Db, Sg, Bh, Hs, Mt, Uun, Uuu, Uub, NUM_ELEMENTS
    };
}


#endif /* DEFINES_H */

