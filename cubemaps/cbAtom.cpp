/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cbAtom.cpp
 * Author: tg
 * 
 * Created on October 8, 2024, 3:35 PM
 */

#include "cbAtom.h"
#include "defines.h"

cbAtom::cbAtom(int Z, int q, double x, double y, double z):
Z_(Z), q_(q), x_(x), y_(y), z_(z) {
}

/**
 convert this cbAtom into a full atom
 */
Atom cbAtom::atom() const {
    // get vdw
    double vdw = PSE::vdw_radii[Z_];
    std::string el = PSE::Elements[Z_];
    double wght = PSE::AtomicWeight[Z_];
    //const std::string& el, const Vec3& xyz, int Z, double weight, double vdw_radius
    Atom atom (el, Vec3(x_, y_, z_), Z_, wght, vdw);
    return atom;
}

