/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Atom.h
 * Author: tg
 *
 * Created on October 11, 2024, 1:58 PM
 */

#ifndef ATOM_H
#define ATOM_H

#include "Vec3.h"

#include <string>

class Atom {
private:
    //! name as in PSE
    std::string element_;
    //! position
    Vec3 xyz_;
    //! ordinal number
    int Z_;
    //! atomic weight
    double weight_;
    //! van-der-Waals radius in A
    double vdw_radius_;
public:
    Atom(const std::string& el, const Vec3& xyz, int Z, double weight, double vdw_radius);
    ~Atom() = default;
    
    std::string element() const { return element_; }
    Vec3 xyz() const { return xyz_; }
    double x() const { return xyz_.x();}
    double y() const { return xyz_.y();}
    double z() const { return xyz_.z();}
    int Z() const { return Z_; }
    double weight () const { return weight_; }
    double vdw_radius() const { return vdw_radius_; }
    
    //! check whether point xyz is within van der waals radius
    bool insphere(const Vec3 xyz) const;

};

#endif /* ATOM_H */

