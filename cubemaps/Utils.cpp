/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <regex>

#include "Utils.h"
#include <iostream>
/**
 * algorithm to create a surface grid at vdW radii of atoms, with a grid spacing 
 * og @c gridspacing (in Angstrom)
 * algortihm from ChatGPT
 * @param atoms
 * @param gridspacing
 * @return 
 */
std::vector<Vec3> surfacegrid (const std::vector<Atom>& atoms, double gridspacing, int verbosity) {
    
    Vec3 bbox_min, bbox_max;
    std::vector<Vec3> grid, surface;
    double minx(9999), miny(9999), minz(9999), maxx(-9999), maxy(-9999), maxz(-9999);

    // find min and max of coordinates
    for (auto atom : atoms) {
        if ((atom.xyz()).x()-atom.vdw_radius() < minx) minx = (atom.xyz()).x()-atom.vdw_radius() ;
        if ((atom.xyz()).x()+atom.vdw_radius() > maxx) maxx = (atom.xyz()).x()+atom.vdw_radius() ;

        if ((atom.xyz()).y()-atom.vdw_radius() < miny) miny = (atom.xyz()).y()-atom.vdw_radius();
        if ((atom.xyz()).y()+atom.vdw_radius() > maxy) maxy = (atom.xyz()).y()+atom.vdw_radius();
        
        if ((atom.xyz()).z()-atom.vdw_radius() < minz) minz = (atom.xyz()).z()-atom.vdw_radius();
        if ((atom.xyz()).z()+atom.vdw_radius() > maxz) maxz = (atom.xyz()).z()+atom.vdw_radius();
    }
    bbox_min = Vec3(minx, miny, minz);
    bbox_max = Vec3(maxx, maxy, maxz);
    
    int stepsx = (maxx-minx)/gridspacing;
    int stepsy = (maxy-miny)/gridspacing;
    int stepsz = (maxz-minz)/gridspacing;
    
    //! create list of grid points
    for (int ix = 0; ix < stepsx; ++ix){
        for (int iy = 0; iy < stepsy; ++iy) {
            for (int iz = 0; iz < stepsz; ++iz) {
                Vec3 xyz (minx + ix*gridspacing, miny + iy*gridspacing, minz + iz*gridspacing);
                grid.push_back(xyz);
                for (auto atom : atoms) {
                    if (atom.insphere(xyz)) {
                        surface.push_back(xyz);
                        break;
                    }
                }
            }
        }
    }
    if (verbosity > 1) {
        std::cout << "---> Generated bounding box with " << grid.size() << " and surface with " 
                << surface.size() << " grid points\n";
    }
    
    return surface;
    
}
