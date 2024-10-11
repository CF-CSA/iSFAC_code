/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Utils.h
 * Author: tg
 *
 * Created on October 11, 2024, 2:12 PM
 */

#ifndef UTILS_H
#define UTILS_H

#include "Vec3.h"
#include "Atom.h"


std::vector<Vec3> surfacegrid (const std::vector<Atom>& atoms, double gridspacing, int verbosity);


#endif /* UTILS_H */

