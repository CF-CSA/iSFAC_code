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
    float a_,b_,c_;
    float alpha_, beta_, gamma_;
    double dhighres_;
    double F000_;
    
    bool centrosymmetric_;

    int verbosity_;
    
    std::vector<Int3x3> symops_;
    bool check_inversion() const;

public:
    FCFInfo(float a, float b, float c, float alpha_, float beta_, float gamma_,
            float dhighres, float F000, const std::vector<Int3x3> symops, int verbosity);
    FCFInfo() = default;
    ~FCFInfo() = default;
    
    float a() const { return a_; }
    float b() const { return b_; }
    float c() const { return c_; }
    float alpha() const { return alpha_; }
    float beta() const { return beta_; }
    float gamma() const { return gamma_; }
    float dhighres() const { return dhighres_; }
    
    float fftscale() const;
    
    std::vector<Int3x3> symops() const { return symops_; }
    // eliminate lattice and inversion operators
    std::vector<Int3x3> symops_no_inv() const;
    int nsymops() const { return symops_.size();}
    
    bool centrosymmetric() const { return centrosymmetric_; }

};

#endif /* FCFINFO_H */

