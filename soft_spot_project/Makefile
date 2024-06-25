F90 = gfortran
F90FLAGS = -O3 -ftree-vectorize -ffast-math -fPIC -lblas -llapack -fopenmp


map    = 'fortran_lib/merge/kind_map'
kinds  = 'fortran_lib/merge/kinds.f95'
kinds2 = 'fortran_lib/merge/kinds2.f95'
system = 'fortran_lib/system/system.f95'
stools = 'fortran_lib/sim_tools/sim_tools.f95'
pot    = 'fortran_lib/potentials/potentials.f95'
cell   = 'fortran_lib/cell_list/cell.f95'
pair   = 'fortran_lib/pairwise/pairwise.f95'
bop    = 'fortran_lib/bond_order/bond_order.f95'



all: test

example.o:
	cat $(kinds) $(system) $(cell) $(stools) $(pot) $(pair) > src.f95
	$(F90) $(F90FLAGS) -c src.f95

mod: example.o
	f90wrap -k $(map) -m sim_lib src.f95
	f2py3 --fcompiler=$(F90)  --opt='$(F90FLAGS)' -lgomp -c -m _sim_lib f90wrap_*.f90 *.o
	f2py3 -c --fcompiler=gfortran  --opt='-O3 -ftree-vectorize -fopenmp' $(kinds2) $(stools) $(bop) -m lib_bond_order

	rm src.f95 *.o *.mod f90wrap*.f90
	sed -i.bu '2s/.*/from . import _sim_lib/' sim_lib.py

test: mod

clean:
	-rm f90wrap*.f90 *.o *.mod *.so

