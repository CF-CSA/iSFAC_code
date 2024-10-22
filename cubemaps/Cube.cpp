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
#include "Utils.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <numeric>

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
        inp >> numAtoms_ >> x >> y >> z;
        origin_ = Vec3(x,y,z);
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

size_t Cube::gridIndex(int ix, int iy, int iz) const {
    size_t idx = iz + Vz_*(iy + Vy_*ix);
    return idx;
}

double Cube::mapValue(int ix, int iy, int iz) const{
    const size_t idx = gridIndex(ix, iy, iz);
    assert (idx >= 0);
    assert (idx < Vx_*Vy_*Vz_);
    double val = gridvalues_[idx];
    return val;
}

/**
 * simple interpolation of grid to get a value
 * coordinates outside grid throws exception (runtime_error)
For interpolating the value at a point \((x, y, z)\) within a 3D cube defined by the eight corner points \((0, 0, 0)\), \((0, 0, 1)\), \((0, 1, 0)\), \((0, 1, 1)\), \((1, 0, 0)\), \((1, 0, 1)\), \((1, 1, 0)\), and \((1, 1, 1)\), a suitable approach is **trilinear interpolation**.

### Trilinear Interpolation Algorithm:

1. **Define the values at the eight corners**:
   Let \( V_{000}, V_{001}, V_{010}, V_{011}, V_{100}, V_{101}, V_{110}, V_{111} \) be the values at the eight corner points:
   - \( V_{000} \) at \( (0, 0, 0) \)
   - \( V_{001} \) at \( (0, 0, 1) \)
   - \( V_{010} \) at \( (0, 1, 0) \)
   - \( V_{011} \) at \( (0, 1, 1) \)
   - \( V_{100} \) at \( (1, 0, 0) \)
   - \( V_{101} \) at \( (1, 0, 1) \)
   - \( V_{110} \) at \( (1, 1, 0) \)
   - \( V_{111} \) at \( (1, 1, 1) \)

2. **Step-by-step interpolation**:
   First, interpolate along the x-axis at the four pairs of points at \((y,z)\) locations:
   V_x(y,z) = (1 - x) \cdot V_{0yz} + x \cdot V_{1yz}
   - For \((y, z) = (0, 0)\), interpolate between \(V_{000}\) and \(V_{100}\):
     V_x(0, 0) = (1 - x) \cdot V_{000} + x \cdot V_{100}
   - For \((y, z) = (0, 1)\), interpolate between \(V_{001}\) and \(V_{101}\):
     V_x(0, 1) = (1 - x) \cdot V_{001} + x \cdot V_{101}
   - For \((y, z) = (1, 0)\), interpolate between \(V_{010}\) and \(V_{110}\):
     V_x(1, 0) = (1 - x) \cdot V_{010} + x \cdot V_{110}
   - For \((y, z) = (1, 1)\), interpolate between \(V_{011}\) and \(V_{111}\):
     V_x(1, 1) = (1 - x) \cdot V_{011} + x \cdot V_{111}

3. **Interpolate along the y-axis**:
   Now, interpolate along the y-axis between the results from step 2:
   V_{xy}(z) = (1 - y) \cdot V_x(0, z) + y \cdot V_x(1, z)
   - For \(z = 0\):
     V_{xy}(0) = (1 - y) \cdot V_x(0, 0) + y \cdot V_x(1, 0)
   - For \(z = 1\):
     V_{xy}(1) = (1 - y) \cdot V_x(0, 1) + y \cdot V_x(1, 1)

4. **Interpolate along the z-axis**:
   Finally, interpolate along the z-axis between the results from step 3:
   V_{xyz} = (1 - z) \cdot V_{xy}(0) + z \cdot V_{xy}(1)

### Result:
The final value \(V_{xyz}\) at the point \((x, y, z)\) is the result of trilinear interpolation.

Would you like me to provide an example calculation or help implement this in code?

 * @param 
 * @return 
 */
double Cube::mapValue(const Vec3& pos) const{
    
    //! get coordinages of pos and ensure its inside grid
    double pos_x = ex_*(pos - origin_);
    double pos_y = ey_*(pos - origin_);
    double pos_z = ez_*(pos - origin_);
    
    // fractional coordinates 
    const double fx = pos_x / std::sqrt(ex_.lengthsq());
    const double fy = pos_y / std::sqrt(ey_.lengthsq());
    const double fz = pos_x / std::sqrt(ez_.lengthsq());
    
    // low left index
    int ix = fx;
    int iy = fy;
    int iz = fz;
    
    if (ix < 0 || ix >= Vx_ || iy < 0 || iy >= Vy_ || iz < 0 || iz > Vz_) {
        if (verbosity_ > 1) {
            std::cout << "*** Error: coordinate " << " out of range\n";
        }
        throw std::logic_error("Index out of range");
    }
    
    const double Vx00 = (1-fx)*mapValue(ix, iy, iz)     + fx*mapValue(ix+1, iy, iz);
    const double Vx01 = (1-fx)*mapValue(ix, iy, iz+1)   + fx*mapValue(ix+1, iy, iz+1);
    const double Vx10 = (1-fx)*mapValue(ix, iy+1, iz)   + fx*mapValue(ix+1,iy+1, iz);
    const double Vx11 = (1-fx)*mapValue(ix, iy+1, iz+1) + fx*mapValue(ix+1, iy+1, iz+1);
    
    const double Vxy0 = (1-fy)*Vx00 + fy*Vx10;
    const double Vxy1 = (1-fy)*Vx01 + fy*Vx11;
    
    const double Vxyz = (1-fz)*Vxy0 + fz*Vxy1;
    
    return Vxyz;
}

/**
 create a list of atoms from cbAtoms
 */
std::vector<Atom> Cube::atoms() const {
    std::vector<Atom> a;
    for (auto c: cbatoms_) {
        a.push_back(c.atom());
    }
    return a;
}

std::vector<double> Cube::distances_sq(const Vec3& pos) const {
    std::vector<double> d2s(0);
    for (auto atom: cbatoms_) {
        const double d = (pos - atom.pos()).lengthsq();
        d2s.push_back(d);
    }
    return d2s;
}

Vec3 Cube::pos(unsigned short idx) const {
    if (idx > numAtoms_) {
        if (verbosity_ > 1) {
            std::cout << "*** Error: atom index out of range\n";
        }
        throw std::logic_error("Atom index out of range");
    }
    return cbatoms_[idx].pos();
}

/**
 * compare coordinates with another map
 * throws an error if different number of atoms
 */
double Cube::deltaTrace(const Cube& cubemap) const {
    if (numAtoms_ != cubemap.numAtoms()) {
        if (verbosity_ > 1) {
            std::cout << "*** Error: comparing Cube map of " << numAtoms_
                    << " atoms with Cube map of " << cubemap.numAtoms() << " atoms.\n";
            throw myExcepts::Format("Unequal number of atoms in maps");
        }
    }
    double trsq = 0.0;
    // create distance matrix
    std::vector<double> distMatrix(numAtoms_*numAtoms_, 0.0);
    for (int i = 0; i < numAtoms_; ++i) {
        double d2 = (cbatoms_[i].pos() - cubemap.pos(i)).lengthsq();
        trsq += d2;
        // std::vector<double> d_row = cubemap.distances_sq(myatom.pos());
        // distMatrix.insert(distMatrix.end(), d_row.begin(), d_row.end());
    }
    return trsq;
}


/**
 * Compute Pearson Correlation coefficient with another cube. Grid values
 * take from this cube and interpolated for the second cube
 * @param cube
 * @return 
 */
double Cube::CC(const Cube& cube) const {
    std::vector<double> g1, g2; 
    g1.reserve(gridvalues_.size());
    g2.reserve(gridvalues_.size());
    // get values from second cube
    for (int ix = 0; ix < Vx_; ++ix) {
        for (int iy = 0; iy < Vy_; ++iy) {
            for (int iz = 0; iz < Vz_; ++iz) {
                const Vec3 pos = origin_ + ix*ex_ + iy*ey_ + iz* ez_;
                try {
                    const double val = cube.mapValue(pos);
                    g1.push_back(mapValue(ix, iy, iz));
                    g2.push_back(val);
                }
                // non-overlapping position
                catch (std::logic_error& e) {
                    continue;
                }
            }
        }
    }
    const double cc =Utils::CC<double>(g1, g2);
    if (verbosity_ > 1) {
        std::cout << "---> Number of overlapping gridpoints: " << g1.size() << '\n'
                << "    Pearson CC for maps: " << cc << std::endl;
    }
    return cc;
}

/**
 Procedure to determine CC between two grid Cfix, Cmove:
 * - cenre of mass CoM for either grid, ctrf, ctrm
 * - move both sets of coordinates to origin
 * - perform Kabsch algorithm to retrieve rotation matrix R that rotates Atoms_mv to
 *   atoms_fix
 * - run through  grid Cmove
 *   -- get val_i and gridpoint g_i
 *   -- compute respective grid point in Cfix
 *   --  R(g_i - ctrm) + ctrf -> coordinate in Cfix
 *   -- get respective value of Cfix
 * - with two lists, compute CC
 */


std::vector<Vec3> Cube::coords() const {
    std::vector<Vec3>X;
    for (auto x: cbatoms_) {
        Vec3 pos = x.pos();
        X.push_back(pos);
    }
    return X;
}

/**
 * Compute centre of geometry for @c N top atoms.
 * If N < 0, use all atoms
 * @param N
 * @return 
 */
Vec3 Cube::CoM(int N) {
    std::vector<Vec3> coords;
    if (N < 0) {
        for (auto x: cbatoms_) {
            Vec3 p = x.pos();
            coords.push_back(p);
        }
    }
    else {
        for (int i = 0; i < N; ++i) {
            Vec3 p = cbatoms_[i].pos();
            coords.push_back(p);
        }
    }
    com_ = 1./coords.size()*std::accumulate(coords.begin(), coords.end(), Vec3(0,0,0));
    return com_;
}

double Cube::prepCC(const Cube& cube) const {
    std::vector<Vec3> coords1 (this->coords());
    std::vector<Vec3> coords2 (cube.coords());
    
    coords1 = Utils::centroid(coords1);
    coords2 = Utils::centroid(coords2);
    
    Mat33 R2_onto1 =Utils::KabschR(coords1, coords2);
    
    
}