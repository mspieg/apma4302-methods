include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

all: reduce_example

reduce_example: reduce_example.o
	-${CLINKER} -o reduce_example reduce_example.o  ${PETSC_LIB}
	${RM} reduce_example.o

clean::
	@rm -f *~ reduce_example *tmp


