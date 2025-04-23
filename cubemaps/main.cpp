/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: tg
 *
 * Created on October 8, 2024, 3:35 PM
 */

#include "Cube.h"
#include "Parser.h"
#include "Usage.h"
#include "myExceptions.h"

/*
 * 
 */
int main(int argc, char** argv) {

    hello(std::cout);

    try {
        Parser parser(argc, argv);

        Cube reference = Cube(parser.cubeRef(), parser.verbosity());
        Cube moving = Cube(parser.cubeMoving(), parser.verbosity());

        //! do some sanity checks before preparing for CC
        consistency_checks(reference, moving, parser.verbosity());
        if (parser.verbosity() > 2) {
            std::cout << "---> Information about maps after consistence checks:\n";
            reference.info();
            moving.info();
        }

        std::pair<Mat33, Vec3> kabschTrafo = moving.makeKabsch(reference);
        moving.transform_coords(kabschTrafo, reference.centroid());
        
        // double cc = reference.CC(moving, kabschTrafo);
        double cc_vdw = reference.CC_VdW(moving, kabschTrafo, parser.vdw_grid_spacing());

    } catch (myExcepts::Usage& e) {
        usage();
    } catch (std::logic_error& e) {
        // just finish
    }
    return 0;
}

