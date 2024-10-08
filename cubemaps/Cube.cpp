/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Cube.cpp
 * Author: tg
 * 
 * Created on October 8, 2024, 3:35 PM
 */

#include "Cube.h"
#include "myExceptions.h"
#include <fstream>
#include <iostream>

Cube::Cube(const std::string& filename, short verbosity) :
verbosity_(verbosity) {
    readMap(filename);
}

void Cube::readMap(const std::string& fname) {

    std::ifstream inp(fname);

    try {
        std::getline(inp, h1_);
        std::getline(inp, h2_);
        double x, y, z;
        inp >> numAtoms_ >> orgx_ >> orgy_ >> orgz_;
        inp >> Vx_ >> x >> y >> z;
        ex_ = Vec3(x, y, z);

        inp >> Vy_ >> x >> y >> z;
        ey_ = Vec3(x, y, z);

        inp >> Vz_ >> x >> y >> z;
        ez_ = Vec3(x, y, z);

        for (int i = 0; i < numAtoms_; ++i) {
            int Z;
            double q;
            double x, y, z;
            inp >> Z >> q >> x >> y >> z;
            if (inp.fail()) {
                std::cout << "---> Failure reading atom.\n";
            } else if (verbosity_ > 3) {
                std::cout << "---> Read atom " << Z << ' ' << q << ' '
                        << x << ' ' << y << ' ' << z << '\n';
            }
            cbatoms_.push_back(cbAtom(Z, q, x, y, z));
        }

        if (inp.eof()) {
            std::cout << "---> Premature end of file\n";
        }
        if (inp.fail()) {

        }
        while ((!inp.eof()) && (!inp.fail())) {
            double g;
            inp >> g;
            if (verbosity_ > 3)
                std::cout << g << '\n';
            gridvalues_.push_back(g);
        }
    }    catch (std::ifstream::failure e) {
        std::cout << "*** Error reading Cube file. Check format\n";
        throw myExcepts::FileIO("Format Cube file");
    }
    inp.close();
    if (verbosity_ > 1) {
        std::cout << "---> Info reading Cube file from " << fname << "\n"
                << "    Number of atoms: " << numAtoms_ << '\n'
                << "    Number of grid points: " << gridvalues_.size() << '\n'
                << "    along x y z: " << Vx_ << ' ' << Vy_ << ' ' << Vz_ 
                << " = " << Vx_ * Vy_ * Vz_ << '\n';
    }
}
