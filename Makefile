SHELL   = /bin/bash
CMP     = ifort # gfortran 
LIBS    = -mkl #-llapack -lblas #${LINK_LAPACK} #-mkl #-llapack -lblas
FRAMEWORK = # -framework Accelerate
FLAGS   = -O3 -r8 #-check all #-no-wrap-margin #-fbounds-check
DEBUG   = #-g -O0 #-ftrapuv #-check all -traceback #-O3
HAN     = -fp-model precise -heap-arrays
OBJrun  = dynamics_main.o arguments.o RK_help.o N_of_grids.o RK_module.o params.o coul_mat_module.o


run.x: dynamics_main.o
	${CMP} -o run.x ${OBJrun} ${LIBS} ${FRAMEWORK} ${OMP} ${DEBUG} ${FLAGS} ${HAN}

dynamics_main.o: dynamics_main.f90 arguments.o RK_help.o N_of_grids.o RK_module.o params.o coul_mat_module.o
	${CMP} -c ${DEBUG}   ${FLAGS} ${HAN} $<

params.o: params.f90 RK_help.o arguments.o N_of_grids.o
	${CMP} -c ${DEBUG}  ${FLAGS} ${HAN} $<

RK_module.o: RK_module.f90 RK_help.o arguments.o N_of_grids.o
	${CMP} -c ${DEBUG}  ${FLAGS} ${HAN} $<

RK_help.o: RK_help.f90 arguments.o N_of_grids.o
	${CMP} -c ${DEBUG}  ${FLAGS} ${HAN} $<

coul_mat_module.o: coul_mat_module.f90 arguments.o N_of_grids.o
	${CMP} -c ${DEBUG}  ${FLAGS} ${HAN} $<

arguments.o: arguments.f90 N_of_grids.o
	${CMP} -c ${DEBUG}  ${FLAGS} ${HAN} $<

N_of_grids.o: N_of_grids.f90 
	${CMP} -c ${DEBUG}  ${FLAGS} ${HAN} $<

clean:
	rm -f *.mod *.o *~ *.exe *.x 

core: run.x
	ulimit -c unlimited; ./run.x	
