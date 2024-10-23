/* 
 * File:   Parser.h
 * Author: tg
 *
 * Created on December 1, 2015, 11:45 AM
 */

#ifndef PARSER_H
#define	PARSER_H

#include <string>
#include <vector>

/**
 * Simplified command line parser
 * @param orig
 */
class Parser {
private:
    std::string cubeRef_, cubeMoving_;
    unsigned char verbosity_;
  
    template <typename T> bool getoption(const std::string& option,
            const std::string& opt, T& optval, int& idx, int argc,
            const char* const argv[]);
    
    // prepare output files if more than 1 input file given
    void make_outfiles();

public:
    Parser(int argc, const char* const argv[]);
    ~Parser(){};
    
    std::string cubeRef() const { return cubeRef_; }
    std::string cubeMoving() const { return cubeMoving_; }
    
    
    unsigned char verbosity() const { return verbosity_; }
};

#endif	/* PARSER_H */

