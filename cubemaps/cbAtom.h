/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cbAtom.h
 * Author: tg
 *
 * Created on October 8, 2024, 3:35 PM
 */

#ifndef CBATOM_H
#define CBATOM_H

/**
 * Atom description from a cube file, containing Z, q, x,y,z
 */
class cbAtom {
private:
    int Z_;
    double q_;
    double x_, y_, z_;
public:
    cbAtom() = default;
    cbAtom(int Z, int q, double x, double y, double z);
    ~cbAtom() = default;

};

#endif /* CBATOM_H */

