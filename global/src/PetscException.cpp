#include "PetscException.hpp"

       
//Positive codes mean that there's an error
//Zero means success
//Negative codes should never happen, but we'll throw anyway
void PetscException(PetscInt petscError,
                    unsigned line,
                    const char* funct,
                    const char* file)
{
    if (petscError != 0)
    {
        const char*  p_text;
        char default_message[30]="Unknown PETSc error code";
        
        //PetscErrorMessage will swing p_text to point to the error code's message
        //...but only if it's a valid code
        PetscErrorMessage(petscError,  &p_text, NULL); 
        if (p_text == 0)
        {
            p_text=default_message;
        }
        
        std::stringstream err_string;
        err_string << p_text;
        err_string << " in function '";
        err_string << funct;
        err_string << "' on line " ;
        err_string << line;
        err_string << " of file ";
        err_string << file;
        
        EXCEPTION(err_string.str());
    }
}

//Positive codes mean that the KSP converged
//Negative codes mean that the KSP diverged i.e. there's a problem
void KspException(PetscInt kspError,
                  unsigned line,
                  const char* funct,
                  const char* file)
{
    if (kspError < 0)
    {
        std::string err_string;

        // The code for the last known error (-10) is hardcoded in PETSc, 
        // in future releases it might change. It is defined in
        // src/ksp/ksp/interface/dlregisksp.c
        if (kspError >= -10 ) err_string = KSPConvergedReasons[kspError];
        else err_string = "Unknown KSP error code";             
                 
        err_string+= " in function '";
        err_string+= funct;
        err_string+= "' on line " ;
        err_string+= line;
        err_string+= " of file ";
        err_string+= file;
        
        EXCEPTION(err_string);
    }
}
