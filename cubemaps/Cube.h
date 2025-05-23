/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Cube.h
 * Author: tg
 *
 * Created on October 8, 2024, 3:35 PM
 */

#ifndef CUBE_H
#define CUBE_H

#include "cbAtom.h"
#include "Vec3.h"
#include "Mat33.h"
#include <vector>
#include <string>

/**
 *Cube map format description from 
 * https://paulbourke.net/dataformats/cube/
 * -> Z is the fast changing direction, X is slowest, i.e.
 * idx = iz + Vz_*(iy + Vy_*ix)
 * positioning:
 * xyz = origin + ix*ex_ + iy*ey_ + iz* ez_;
 */
class Cube {
private:
    //! two header lines
    std::string h1_, h2_;
    //! number of atoms
    int numAtoms_;
    //! origin from cube file
    Vec3 origin_;
    //! centroid from coordinates computed in kabsch algorithm
    Vec3 centroid_;
    //! number of voxels in each direction
    int Vx_, Vy_, Vz_;
    //! directions of grid
    Vec3 ex_, ey_, ez_;
    //! cross product of grid directions: ux_ = ey_ x ez_, normalised to ex_*ux_=1 etc.
    Vec3 ux_, uy_, uz_;
    //! mininum and maximum coordinates of atoms:
    double minx_, maxx_, miny_, maxy_, minz_, maxz_;
    
    std::vector<cbAtom> cbatoms_;
    std::vector<double> gridvalues_;
    short verbosity_;

    void readMap(const std::string& fname);
    //! compute geometric centroid ('CoM') for N coordinates
    Vec3 calc_centroid(int N);

public:
    Cube() = default;
    Cube(const std::string& filename, short verbosity_);
    ~Cube() = default;
    
    //! compute index from three grid coordinates
    size_t gridIndex(int ix, int iy, int iz) const;

    //! getter functions
    int numAtoms() const { return cbatoms_.size(); }
    
    //! get grid value at coordinate
    double mapValue(const Vec3&) const;
    //! get grid value by indices
    double mapValue(int ix, int iy, int iz) const;
    
    // return position of atom idx
    Vec3 pos(unsigned short idx) const;
    
    Vec3 centroid() const { return centroid_; }
    // set centroid with three coordinates
    void centroid(double x, double y, double z);

    // list of atoms
    std::vector<Atom> atoms() const;
    // only list of coordinates
    std::vector<Vec3> coords() const;
    
    // return list of distances to position
    std::vector<double> distances_sq(const Vec3& pos) const;
    
    // compute the trace of moduli of distances^2
    double deltaTrace(const Cube& cubemap) const;
    
    //! get Kabsch rotation to rotate other cube onto this one
    Mat33 makeKabsch(const Cube& cube, Vec3& ctrd_this, Vec3& ctrd_other) const;
    
    //! print coordinates before and after moving with Kabsch transform
    void transform_coords(const Mat33& kabschTrafo, const Vec3& ctr_target);
    
    //! compute pearson coefficient with a second grid
    double CC(const Cube& other, const Mat33& KabschTrafo) const;
    
    //! compute Pearson coefficient with second cube on VdW surface
    double CC_VdW(const Cube& other, const Mat33& KabschTrafo, const double& vdw_grid_spacing) const;
    
    //! print some information about map
    void info() const;
    
    //! make some checks, e.g. compare number of atoms, etc.
    friend bool consistency_checks (Cube& one, Cube& other, unsigned char verbosity);

};

#endif /* CUBE_H */

