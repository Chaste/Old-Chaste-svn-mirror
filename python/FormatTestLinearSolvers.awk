BEGIN {USER=0}
IF $1=="RUN:" && USER==1 {print; print$0; COUNTER=0; USER=1}
IF $1=="RUN:" && USER==0 {print $0; COUNTER=0; USER=1}
COUNTER++ { }
IF COUNTER==13 && $1=="Petsc" {print "Error"}
IF COUNTER==13 && $1!="Petsc" {print $0}
IF $1=="user" {print $2; USER=0}
