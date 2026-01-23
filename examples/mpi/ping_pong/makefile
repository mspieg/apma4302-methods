include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

all: ping_pong

ping_pong: ping_pong.o
	-${CLINKER} -o ping_pong ping_pong.o  ${PETSC_LIB}
	${RM} ping_pong.o


clean::
	@rm -f *~ ping_pong *tmp



