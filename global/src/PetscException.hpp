#ifndef PETSCEXCEPTION_HPP_
#define PETSCEXCEPTION_HPP_

#include <petsc.h>
#include <petscksp.h>

#include "Exception.hpp"

extern void PetscException(int petscError, int line,
                           const char* funct, const char* file);

extern void KspException(int kspError, int line,
                         const char* funct, const char* file);

//Positive codes mean that there's an error
//Zero means success
//Negative codes should never happen, but we'll throw anyway
#define PETSCEXCEPT(n) if (n) PetscException(n, __LINE__, __FUNCT__,__FILE__)

//Positive codes mean that the KSP converged
//Negative codes mean that the KSP diverged i.e. there's a problem
#define KSPEXCEPT(n) if (n<0) KspException(n, __LINE__, __FUNCT__,__FILE__)

#endif /*PETSCEXCEPTION_HPP_*/
