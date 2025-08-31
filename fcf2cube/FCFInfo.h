/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FCFInfo.h
 * Author: tg
 *
 * Created on July 30, 2023, 9:51 PM
 */

#ifndef FCFINFO_H
#define FCFINFO_H

#include "Int3x3.h"

#include <string>
#include <vector>


class FCFInfo {
    double a_,b_,c_;
    double alpha_, beta_, gamma_;
    double dhighres_;
    double F000_;
    
    bool centrosymmetric_;

    int verbosity_;
    
    std::vector<Int3x3> symops_;

public:
    FCFInfo(double a, double b, double c, double alpha_, double beta_, double gamma_,
            double dhighres, double F000, const std::vector<Int3x3> symops, int verbosity);
    FCFInfo() = default;
    ~FCFInfo() = default;
    
    double a() const { return a_; }
    double b() const { return b_; }
    double c() const { return c_; }
    double alpha() const { return alpha_; }
    double beta() const { return beta_; }
    double gamma() const { return gamma_; }
    double dhighres() const { return dhighres_; }
    
    double fftscale() const;
    
    std::vector<Int3x3> symops() const { return symops_; }
    // eliminate lattice and inversion operators
    std::vector<Int3x3> symops_no_inv();
    int nsymops() const { return symops_.size();}
    
    bool centrosymmetric() const { return centrosymmetric_; }

};

#endif /* FCFINFO_H */

