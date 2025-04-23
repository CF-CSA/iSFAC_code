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
#include "kabsch.h"
#include "defines.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <numeric>
#include <iomanip>

Cube::Cube(const std::string& filename, short verbosity) :
verbosity_(verbosity) {
    readMap(filename);
    // when successful, compute centroid
    centroid();
}

/**
 * According to http://gaussian.com/cubegen/
 * Vx is the slow direction, and Vz is the fast direction
 * in the first line, there could be a number for NVal, number of points per grid
 * this might be useful for vector fields. My example cube files omit this entry, 
 * and the below code doesn't capture both formats.
 * @param fname
 */
void Cube::readMap(const std::string& fname) {

    std::ifstream inp(fname);
    // scale for Bohr radius
    double a0toA = Physics::a0;

    try {
        std::getline(inp, h1_);
        std::getline(inp, h2_);
        double x, y, z;
        double orgx, orgy, orgz;
        inp >> numAtoms_ >> orgx >> orgy >> orgz;
        
        inp >> Vx_ >> x >> y >> z;
        if (Vx_ < 0) {
            ex_ = Vec3(x, y, z);
            
        } else {
            ex_ = Vec3(a0toA*x, a0toA*y, a0toA * z);
            orgx *= a0toA;
        }

        inp >> Vy_ >> x >> y >> z;
        if (Vy_ < 0) { // unit already in A
            ey_ = Vec3(x, y, z);
        } else {
            ey_ = Vec3(a0toA*x, a0toA*y, a0toA * z);
            orgy *= a0toA;
        }

        inp >> Vz_ >> x >> y >> z;
        if (Vz_ < 0) {// unit already in A
            ez_ = Vec3(x, y, z);
        } else {
            ez_ = Vec3(a0toA*x, a0toA*y, a0toA * z);
            orgz *= a0toA;
        }
        origin_ = Vec3(orgx, orgy, orgz);


        if (Vx_ < 0 && Vy_ < 0 && Vz_ < 0) {
            a0toA = 1.0;
        }

        for (int i = 0; i < numAtoms_; ++i) {
            int Z;
            double q;
            double x, y, z;
            inp >> Z >> q >> x >> y >> z;
            if (Vx_ > 0) x *= a0toA;
            if (Vy_ > 0) y *= a0toA;
            if (Vz_ > 0) z *= a0toA;
            if (inp.fail()) {
                std::cout << "---> Failure reading atom.\n";
                throw myExcepts::Format("Premature end of Atom list in Cube file");
            } else if (verbosity_ > 2) {
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
            if (inp.fail() || inp.eof()) {
                break;
            }
            if (verbosity_ > 3)
                //std::cout << g << '\n';
                gridvalues_.push_back(g);
        }
    } catch (std::ifstream::failure e) {
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
    // define ``reciprocal vectors'' for coordinate retrieval
    ux_ = cross(ey_, ez_);
    uy_ = cross(ez_, ex_);
    uz_ = cross(ex_, ey_);
}

size_t Cube::gridIndex(int ix, int iy, int iz) const {
    size_t idx = iz + Vz_ * (iy + Vy_ * ix);
    return idx;
}

double Cube::mapValue(int ix, int iy, int iz) const {
    const size_t idx = gridIndex(ix, iy, iz);
    assert(idx >= 0);
    assert(idx < Vx_ * Vy_ * Vz_);
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
double Cube::mapValue(const Vec3& pos) const {

    //! get coordinages of pos and ensure its inside grid
    double fx = (pos - origin_) * ux_ / (ex_ * ux_);
    double fy = (pos - origin_) * uy_ / (ey_ * uy_);
    double fz = (pos - origin_) * uz_ / (ez_ * uz_);

    // low left index
    int ix = std::round(fx);
    int iy = std::round(fy);
    int iz = std::round(fz);


    if (ix < 0 || ix + 1 >= Vx_ || iy < 0 || iy + 1 >= Vy_ || iz < 0 || iz + 1 >= Vz_) {
        throw std::logic_error("Index out of range");
    }

    const double val = mapValue(ix, iy, iz);
    if (verbosity_ > 2) {
        std::cout << "-----> Reading map value at position " << pos
                << ", converted to indices "
                << std::setw(4) << ix
                << std::setw(4) << iy
                << std::setw(4) << iz
                << " = " << val << '\n';
    }
    // debugging (2024/11/07) just use the value at ix, iy, iz
    return val;

    const double Vx00 = (1 - fx) * mapValue(ix, iy, iz) + fx * mapValue(ix + 1, iy, iz);
    const double Vx01 = (1 - fx) * mapValue(ix, iy, iz + 1) + fx * mapValue(ix + 1, iy, iz + 1);
    const double Vx10 = (1 - fx) * mapValue(ix, iy + 1, iz) + fx * mapValue(ix + 1, iy + 1, iz);
    const double Vx11 = (1 - fx) * mapValue(ix, iy + 1, iz + 1) + fx * mapValue(ix + 1, iy + 1, iz + 1);

    const double Vxy0 = (1 - fy) * Vx00 + fy*Vx10;
    const double Vxy1 = (1 - fy) * Vx01 + fy*Vx11;

    const double Vxyz = (1 - fz) * Vxy0 + fz*Vxy1;

    return Vxyz;
}

/**
 create a list of atoms from cbAtoms
 */
std::vector<Atom> Cube::atoms() const {
    std::vector<Atom> a;
    for (auto c : cbatoms_) {
        a.push_back(c.atom());
    }
    return a;
}

std::vector<double> Cube::distances_sq(const Vec3& pos) const {
    std::vector<double> d2s(0);
    for (auto atom : cbatoms_) {
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
 * KabschTrafo must be computed for this map
 * @param cube
 * @return 
 */
double Cube::CC(const Cube& other, const std::pair<Mat33, Vec3>& KabschTrafo) const {
    const Mat33 U = KabschTrafo.first;
    const Vec3 T = KabschTrafo.second;
    std::vector<double> g1, g2;
    //! ensure centroid is properly computed
    //
    g1.reserve(gridvalues_.size());
    g2.reserve(gridvalues_.size());

    // iterate through the reference cube convert point to respective position
    // in other cube and get value from there
    for (int ix = 0; ix < Vx_; ++ix) {
        for (int iy = 0; iy < Vy_; ++iy) {
            for (int iz = 0; iz < Vz_; ++iz) {
                // compute coordinate of current grid point
                Vec3 pos = origin_ + ix * ex_ + iy * ey_ + iz* ez_;
                if (verbosity_ > 2) {
                    std::cout << "-----> Indices "
                            << std::setw(3) << ix
                            << std::setw(3) << iy
                            << std::setw(3) << iz
                            << " correspond to position " << pos;
                }
                pos = pos - centroid_;
                // compute position of grid index corresponding to position
                // in other map
                pos = U * pos + T + other.centroid();
                if (verbosity_ > 1) {
                    std::cout << " after trafo: " << pos << '\n';
                }
                try {
                    double val = other.mapValue(pos);
                    g2.push_back(val);
                    if (verbosity_ > 1) {
                        std::cout << "----> value from moving map at "
                                << pos << " = " << val << '\n';
                    }
                    val = mapValue(ix, iy, iz);
                    if (verbosity_ > 1) {
                        std::cout << "---> value from reference map at "
                                << ix << ' ' << iy << ' ' << iz << " = " << val << '\n';
                    }
                    // only push first value if transform pos is inside second grid
                    g1.push_back(val);
                }// non-overlapping position
                catch (std::logic_error& e) {
                    if (verbosity_ > 2) {
                        std::cout << "----> no matching grid point\n";
                    }
                    continue;
                }
            }
        }
    }
    const double cc = Utils::CC<double>(g1, g2);
    if (verbosity_ > 1) {
        std::cout << "---> Number of overlapping gridpoints: " << g1.size() << '\n'
                << "    Pearson CC for maps: " << cc << std::endl;
    }
    const double cc_gsl = Utils::CC_gsl(g1, g2);
    if (verbosity_ > 1) {
        std::cout << "---> Number of overlapping gridpoints: " << g1.size() << '\n'
                << "    Pearson CC for maps from GSL CC: " << cc_gsl << std::endl;
    }
    return cc;
}

/**
 * Compute Pearson Correlation coefficient with another cube. Grid values
 * take from this cube and interpolated for the second cube
 * KabschTrafo must be computed for this map
 * @param cube
 * @return 
 */
double Cube::CC_VdW(const Cube& other, const std::pair<Mat33, Vec3>& KabschTrafo) const {
    const Mat33 U = KabschTrafo.first;
    const Vec3 T = KabschTrafo.second;
    std::vector<double> g1, g2;
    
    // generate VdW surface for cube
    std::vector<Atom> atoms1;
    for (auto x: cbatoms_) {
        atoms1.push_back(x.atom());
    }
    std::vector<Vec3> vdw1 = Utils::surfacegrid(atoms1, 1.1, verbosity_);
    for (auto x: vdw1) {
        Vec3 pos (x.x(), x.y(), x.z());
        try {
        double v1 = mapValue(pos);
        double v2 = other.mapValue(pos);
        g1.push_back(v1);
        g2.push_back(v2);
        }
        catch(std::logic_error& e) {
            continue;
        }
    }

    const double cc = Utils::CC<double>(g1, g2);
    if (verbosity_ > 1) {
        std::cout << "---> Number of overlapping gridpoints: " << g1.size() << '\n'
                << "    Pearson CC for maps: " << cc << std::endl;
    }
    const double cc_gsl = Utils::CC_gsl(g1, g2);
    if (verbosity_ > 1) {
        std::cout << "---> Number of overlapping gridpoints: " << g1.size() << '\n'
                << "    Pearson CC for maps from GSL CC: " << cc_gsl << std::endl;
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
    for (auto x : cbatoms_) {
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
Vec3 Cube::centroid(int N) {
    std::vector<Vec3> coords;
    if (N < 0) {
        for (auto x : cbatoms_) {
            Vec3 p = x.pos();
            coords.push_back(p);
        }
    } else {
        for (int i = 0; i < N; ++i) {
            Vec3 p = cbatoms_[i].pos();
            coords.push_back(p);
        }
    }
    centroid_ = 1. / coords.size() * std::accumulate(coords.begin(), coords.end(), Vec3(0, 0, 0));
    if (verbosity_ > 1) {
        std::cout << "---> Cntroid: " << centroid_ <<'\n';
    }
    return centroid_;
}

/**
 * Determine how to move @c other to this:
 * move to centroid, get R from Kabsch, move to centroid to this
 * @param cube
 * @return 
 */
std::pair<Mat33, Vec3> Cube::makeKabsch(const Cube& other) const {

    // prepare calling kabsch function
    gsl_matrix* fixed;
    gsl_matrix* moved;
    // compute Kabsch transform where there is fixed, and this is moved
    // that way, iterating through the map of fixed can be transformed to 
    // moved

    fixed = gsl_matrix_alloc(cbatoms_.size(), 3);
    size_t idx(0);
    for (auto x : cbatoms_) {
        if (x.Z() == 1) continue;
        gsl_matrix_set(fixed, idx, 0, x.pos().x());
        gsl_matrix_set(fixed, idx, 1, x.pos().y());
        gsl_matrix_set(fixed, idx, 2, x.pos().z());
        ++idx;
    }

    moved = gsl_matrix_alloc(cbatoms_.size(), 3);
    idx = 0;
    for (auto x : other.cbatoms_) {
        if (x.Z() == 1) continue;
        gsl_matrix_set(moved, idx, 0, x.pos().x());
        gsl_matrix_set(moved, idx, 1, x.pos().y());
        gsl_matrix_set(moved, idx, 2, x.pos().z());
        ++idx;
    }

    gsl_matrix* U = gsl_matrix_alloc(3, 3);
    gsl_vector* t = gsl_vector_alloc(3);

    kabsch(cbatoms_.size(), moved, fixed, U, t, NULL);
    // Mat33 kabschR =Utils::KabschR(coords1, coords2);
    Mat33 kabschR;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            kabschR.operator()(r, c) = gsl_matrix_get(U, r, c);
        }
    }

    Vec3 kabschT(gsl_vector_get(t, 0), gsl_vector_get(t, 1), gsl_vector_get(t, 2));


    std::cout << "---> Det(U) = " << Utils::determinant(U) << std::endl;
    std::cout << "---> KR = \n" << kabschR << std::endl;
    std::cout << "---> Translation: " << gsl_vector_get(t, 0) << ' '
            << gsl_vector_get(t, 1) << ' '
            << gsl_vector_get(t, 2) << ' '
            << std::endl;

    std::pair<Mat33, Vec3> K(kabschR, kabschT);

    gsl_matrix_free(fixed);
    gsl_matrix_free(moved);
    gsl_matrix_free(U);
    gsl_vector_free(t);

    //return is Kabsch matrix and inverse of centre
    return K;
}

void Cube::transform_coords(const std::pair<Mat33, Vec3>& kabschTrafo, const Vec3& ctr_target) {
    const Mat33 M = kabschTrafo.first;
    const Vec3 T = kabschTrafo.second;
    for (const auto& a : cbatoms_) {
        Vec3 x = a.pos() - centroid_;
        x = M * x + T + centroid_;
        std::cout << a.Z() << " before: " << a.pos() << " after: " << x << '\n';
    }
}

/**
 * - compare number of atoms, keep only lower number
 * - compute both distance matrices, and compare them element wise
 * - recompute centroid
 * @param other
 * @return 
 */
bool consistency_checks(Cube& one, Cube& two, unsigned char verbosity) {
    int N = std::min(one.numAtoms(), two.numAtoms());

    if (verbosity > 1) {
        std::cout << "---> Consistency check: limiting number of atoms in both maps to " << N << '\n';
    }
    one.cbatoms_.resize(N);
    two.cbatoms_.resize(N);

    one.centroid(-1);
    two.centroid(-1);

    one.info();
    two.info();
    if (verbosity > 1) {
        std::cout << "---> Producing distance matrix for first map with " << one.coords().size() << " coordinates.\n";
    }
    const std::vector<double> d2one(Utils::distance_matrix(one.coords()));
    if (verbosity > 1) {
        std::cout << "---> Producing distance matrix for second map with " << two.coords().size() << " coordinates.\n";
    }
    const std::vector<double> d2two(Utils::distance_matrix(two.coords()));
    // compute sum pairwise distances
    std::vector<double>::const_iterator it1, it2;
    it2 = d2two.begin();
    double sum(0.0);
    for (const auto x1 : d2one) {
        sum += std::abs(x1 - (*it2));
        ++it2;
    }
    if (verbosity > 1) {
        std::cout << "---> Sum of absolute differences between distance matrices: " << sum << '\n';
    }

    return true;

}

void Cube::info() const {
    std::vector<double> dists = Utils::distance_matrix(coords());
    std::cout << "---> Information about Cube map\n"
            << "   Number of atoms: " << cbatoms_.size() << '\n';
    for (int idx = 0; idx < dists.size();) {
        const int tabs = idx / (cbatoms_.size() - 1) + 1;
        for (int i = 0; i < tabs; ++i) {
            std::cout << std::setw(7) << ' ';
        }
        for (int i = tabs; i < cbatoms_.size(); ++i) {
            std::cout << std::fixed << std::setw(7) << std::setprecision(2) << dists[idx];
            ++idx;
        }
        std::cout << '\n';
    }
}
