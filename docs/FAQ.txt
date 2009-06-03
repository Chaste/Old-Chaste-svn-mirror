
Q: I have problems compiling some libraries, the error looks to be caused by HDF5/MPI

A: Mpi redefines std:: c headers. Try this:
Add
  #undef SEEK_SET
  #undef SEEK_END
  #undef SEEK_CUR
before
  #include "mpi.h"
in $HOME/hdf5/include/H5public.h

