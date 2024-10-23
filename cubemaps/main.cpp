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
    Cube moving    = Cube(parser.cubeMoving(), parser.verbosity());
    }
    catch (myExcepts::Usage& e){
        usage();
    }
    catch (std::logic_error& e) {
        // just finish
    }
    return 0;
}

