f95 -O3 read_10.f metap.f iter_10.f ../WORK/Fortran/src/dranxor.f; 
f95 -O3 read_fit_log.f -o b.out;
ipython loop.py;
