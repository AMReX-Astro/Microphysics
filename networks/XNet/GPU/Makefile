# Makefile for source/physics/sourceTerms/Burn/nuclearBurn/XNet/GPU

Burn += xnet_gpu_control.o cublasf.o cudaf.o

xnet_jacobian_dense.o : xnet_gpu_control.o xnet_types.o cublasf.o cudaf.o
cudaf.F90 : cudaDeviceProp.fh

## The cudaDeviceProp struct can change between CUDA versions.
## This can cause hard-to-detect stack corruption due to a difference in size
## with the Fortran interoperable derived type differs in size.
## Here, we grab the cudaDeviceProp sturct from the current CUDA header file and
## convert it to a Fortran interoperable derived type at compile time.
cudaDeviceProp.fh : ${CUDA_PATH}/include/driver_types.h
	@sed -n '/struct.*cudaDeviceProp/,/}/ p' ${CUDA_PATH}/include/driver_types.h | \
		sed -e '1 d' \
		    -e 's/^{$\/TYPE, BIND(C) :: cudaDeviceProp/' \
		    -e 's/^\s*\<char\>\s*/    CHARACTER(C_CHAR) :: /' \
		    -e 's/^\s*\<size_t\>\s*/    INTEGER(C_SIZE_T) :: /' \
		    -e 's/^\s*\<int\>\s*/    INTEGER(C_INT)    :: /' \
		    -e 's/\[\([0-9]\+\)\]/(\1)/' \
		    -e 's/;\(\s*\)\/\*\(.*\)\*\// \1!\2/' \
		    -e 's/;.*$\//' \
		    -e 's/^}$\/END TYPE cudaDeviceProp/' > $@
