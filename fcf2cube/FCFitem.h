/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FCFitem.h
 * Author: tg
 *
 * Created on July 30, 2023, 9:17 PM
 */

#ifndef FCFITEM_H
#define FCFITEM_H

#include "HKL.h"
class FCFitem {
        HKL hkl_;
        double Imeas_;
        double sigImeas_;
        double Fcalc_;
        double phicalc_;

public:
    FCFitem() = default;
    FCFitem(const HKL& hkl, const double& I, const double& sigI, const double& Fc, const double& phic);
    ~FCFitem() = default ;
    
    HKL& hkl(const HKL& hkl) { hkl_ = hkl; return hkl_;}
    double& Imeas(const double I) { Imeas_ = I; return Imeas_; }
    double& sigImeas(const double sigI) { sigImeas_ = sigI; return sigImeas_; }
    double& Fcalc(const double Fc) { Fcalc_ = Fc; return Fcalc_; }
    double& phicalc(const double phi) { phicalc_ = phi; return phicalc_;}

    HKL hkl() const { return hkl_; }
    double Imeas() const { return Imeas_; }
    double sigImeas() const{ return sigImeas_; }
    double Fcalc() const { return Fcalc_; }
    double phicalc() const { return phicalc_;}
    
    friend bool operator<(const FCFitem& f1, const FCFitem& f2);
    friend std::ostream& operator<< (std::ostream& outp, const FCFitem& fcfitem);

};

#endif /* FCFITEM_H */

