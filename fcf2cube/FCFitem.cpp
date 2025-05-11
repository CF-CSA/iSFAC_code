/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FCFitem.cpp
 * Author: tg
 * 
 * Created on July 30, 2023, 9:17 PM
 */

#include "FCFitem.h"

#include <iomanip>
#include <iostream>

FCFitem::FCFitem(const HKL& hkl, const double& I, const double& sigI, const double& Fc, const double& phic):
hkl_(hkl), Imeas_(I), sigImeas_(sigI), Fcalc_(Fc), phicalc_(phic){
}

/**
 comparison operator for sorting; passes through operator< for HKL
 */
bool operator< (const FCFitem& f1, const FCFitem& f2) {
    return f1.hkl_ < f2.hkl_;
}


