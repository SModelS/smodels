* types.h
* real-based type declarations
* this file is part of LoopTools
* last modified 9 Jul 12 th


#ifndef TYPES_H
#define TYPES_H

#define RealType double precision
#define ComplexType double complex
#define Re DBLE
#define Im DIMAG
#define Conjugate DCONJG
#define ToComplex DCMPLX

#define Sq(c) Re((c)*Conjugate(c))
#define Sqrtc(c) sqrt(ToComplex(c))

#endif

