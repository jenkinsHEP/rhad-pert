#
# usage:
# 	make                          (compile rhad)
#	make prog=<dum>               (compile and call executable x<dum>)
#	make clean                    (delete object files and executable)
#

ifndef FORTRAN
#  FORTRAN = f77
  FORTRAN = gfortran
endif

ifndef prog
  prog = example
endif

SHELL = /usr/bin/tcsh
OPTS =
SRC = ${prog}.f r012.f r34.f runal.f funcs.f vegas-rhad.f \
      parameters.f rhad.f
OBJ = ${SRC:.f=.o}


path=${HOME}/f
f77 = $(FORTRAN)

x${prog}: ${OBJ} common.f
	${f77} ${OPTS} -o $@ ${OBJ}

%.o: %.f
	${f77} ${OPTS} -c -o $@ $<

clean:
	if (-e x${prog}) rm x${prog};\
	rm *.o



