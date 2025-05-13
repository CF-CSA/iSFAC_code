/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: tg
 *
 * Created on May 9, 2025, 9:11 PM
 */

#include <iostream>
#include "Parser.h"
#include "FCFfile.h"
#include "Usage.h"
#include "ResFile.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    
    Parser parser(argc, argv);
    if (parser.verbosity() > 0) {
        hello(std::cout);
    }
    FCFfile fcffile(parser.fcffile(), parser.f000(), parser.verbosity());
    ResFile resfile(parser.resfile(), parser.read_qs(), parser.verbosity());

    return 0;
}

