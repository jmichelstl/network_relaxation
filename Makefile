CC=g++ --std=c++11 -O3 -fexceptions -fPIC -fopenmp -DGPU_BLAS
MKLCC=g++ --std=c++11 -O3 -fexceptions -m64
INC=-I/usr/include/mkl -I/usr/include/suitesparse
MKLINC=-I/usr/include/mkl
GPULINK= -lsuitesparseconfig -lSuiteSparse_GPURuntime -lGPUQREngine
MKLLINK= -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/mkl/

all:
	make fiber_relax
	make libflists.a
	make sp_relax
	make irelax
	make grelax

fiber_relax: libflists.a
	$(MKLCC) $(MKLINC)  src/relax_periodic_network.cpp -o fiber_relax \
	$(MKLLINK) -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_tbb_thread \
       	-lmkl_core -ltbb -lpthread -lm -ldl libflists.a -lboost_filesystem \
	-lboost_regex

sp_relax: libflists.a
	$(CC) $(INC) src/sparse_relax.cpp -o sp_relax libflists.a \
	-L/usr/local/lib -lcholmod -lspqr $(GPU_LINK) -lblas \
	-lboost_filesystem -lboost_regex

grelax: libflists.a
	nvcc -gencode=arch=compute_75,code=sm_75 -I/usr/include/mkl \
	src/g_relax_periodic_network.cu -std=c++11 libflists.a \
	-lboost_filesystem -lboost_regex -lcudart \
	-lcublas -lcusparse -L/usr/lib/cuda-10.2 $(MKLLINK) -lblas64 \
        -o grelax

libflists.a:
	g++ --std=c++11 -c src/file_lists.cpp -o libflists.a


irelax: libflists.a
	$(MKLCC) -I/usr/include/mkl src/mkl_sparse_relax.cpp -o irelax \
	libflists.a $(MKLLINK) /usr/lib/x86_64-linux-gnu/mkl/libblas.so \
	-lmkl_core  -lboost_filesystem -lboost_regex

lsrelax: libflists.a
	nvcc -gencode=arch=compute_75,code=sm_75 src/large_strain_relax.cu -std=c++11 -lcublas -lm libflists.a -lboost_filesystem -lboost_regex -o lsrelax

nmrelax: libflists.a
	nvcc -gencode=arch=compute_75,code=sm_75 src/network_mesh_relax.cu -std=c++11 -lcublas -lm libflists.a -lboost_filesystem -lboost_regex -o nmrelax

mkltest:
	$(MKLCC) -I/usr/include/mkl src/mkl_sp_qr_test.cpp -o mkltest \
	$(MKLLINK) -lmkl_core /usr/lib/x86_64-linux-gnu/mkl/libblas64.so

nacalc:
	g++ --std=c++11 -O3 -Werror -o nacalc src/na_calc.cpp \
	       	-lboost_filesystem -lboost_regex

straincalc:
	g++ --std=c++11 -O3 -Werror -o scalc src/strain_calc.cpp

nacorr:
	g++ --std=c++11 -O3 -Werror -o nacorr src/na_corr.cpp -lboost_regex -lboost_filesystem -lgsl -lblas -lfftw3

edensity:
	g++ --std=c++11 -O3 -Werror -o edensity src/edensity_analysis.cpp \
		-lboost_filesystem -lboost_regex

fclean:
	rm fiber_relax

sclean:
	rm sp_relax

lclean:
	rm libflists.a

iclean:
	rm irelax

lsclean:
	rm lsrelax

naclean:
	rm nacalc

stclean:
	rm straincalc

gclean:
	rm grelax

ncorrclean:
	rm nacorr

edclean:
	rm edensity

nmclean:
	rm nmrelax

clean:
	rm fiber_relax
	rm sp_relax
	rm libflists.a
	rm irelax
	rm nacalc
	rm straincalc
	rm grelax
	rm nacorr
	rm edensity
	rm nmrelax
