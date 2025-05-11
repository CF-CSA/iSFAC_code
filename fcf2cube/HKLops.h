/* 
 * File:   hklops.h
 * Author: tg
 *
 * Created on July 3, 2012, 6:49 PM
 */

#ifndef HKLOPS_H
#define	HKLOPS_H

class HKL;
class Int3x3;
class FCFitem;

/* forward operators of basic operations for other files e.g. for sorting */
bool operator< (const HKL& h1, const HKL& h2);
bool operator== (const HKL& h1, const HKL& h2);
bool operator> (const HKL& h1, const HKL& h2);
HKL operator+(const HKL& h1, const HKL& h2);
HKL operator-(const HKL& h1, const HKL& h2);
HKL operator*(const Int3x3& R, const HKL& h);
std::ostream& operator<<(std::ostream& outp, const FCFitem& fcfitem);
std::ostream& operator<<(std::ostream& outp, const HKL& hkl);


#endif	/* HKLOPS_H */

