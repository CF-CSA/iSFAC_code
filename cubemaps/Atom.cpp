/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Atom.cpp
 * Author: tg
 * 
 * Created on October 11, 2024, 1:58 PM
 */

#include "Atom.h"

Atom::Atom(const std::string& el, const Vec3& xyz, int Z, double weight, double vdw_radius):
element_(el), xyz_(xyz), Z_(Z), weight_(weight), vdw_radius_(vdw_radius) {
    
}


/*
 * return true if point xyz is within vdW sphere of this atom
 */
bool Atom::insphere(const Vec3 xyz) const {
    double d2 = (xyz_ - xyz).lengthsq();
    return (d2 < vdw_radius_*vdw_radius_);
}