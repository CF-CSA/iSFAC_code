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

        if (parser.verbosity() > 2) {
            std::cout << "---> Information about maps before consistence checks:\n";
            reference.info();
            moving.info();
        }
        //! do some sanity checks before preparing for CC
        consistency_checks(reference, moving, parser.verbosity());
        if (parser.verbosity() > 2) {
            std::cout << "---> Information about maps after consistence checks:\n";
            reference.info();
            moving.info();
        }
        
        Mat33 kabschR = reference.makeKabsch(moving);
        
        double cc = reference.CC(moving, kabschR);
        
        std::cout << "---> CC between both maps: " << cc << '\n'
                << "    with Kabsch Matrix \n" << kabschR << '\n';
        
        //! check for ID
        kabschR = reference.makeKabsch(reference);
        cc = reference.CC(reference, kabschR);
        std::cout << "---> CC between itself: " << cc << '\n'
                << "    with Kabsch Matrix \n" << kabschR << '\n'
                << "    det(K) = " << kabschR.determinant() << '\n';
        
    }    catch (myExcepts::Usage& e) {
        usage();
    }    catch (std::logic_error& e) {
        // just finish
    }
    return 0;
}

