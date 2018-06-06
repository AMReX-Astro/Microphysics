MODULE cublasf
!===============================================================================
! INTERFACE to CUBLAS routines
!===============================================================================
    USE, INTRINSIC :: ISO_C_BINDING

    ENUM, BIND(C) !:: cublasStatus_t
    ENUMERATOR :: CUBLAS_STATUS_SUCCESS=0
    ENUMERATOR :: CUBLAS_STATUS_NOT_INITIALIZED =1
    ENUMERATOR :: CUBLAS_STATUS_ALLOC_FAILED=3
    ENUMERATOR :: CUBLAS_STATUS_INVALID_VALUE=7
    ENUMERATOR :: CUBLAS_STATUS_ARCH_MISMATCH=8
    ENUMERATOR :: CUBLAS_STATUS_MAPPING_ERROR=11
    ENUMERATOR :: CUBLAS_STATUS_EXECUTION_FAILED=13
    ENUMERATOR :: CUBLAS_STATUS_INTERNAL_ERROR=14
    END ENUM !cublasStatus_t

    ENUM, BIND(C) !:: cublasFillMode_t
    ENUMERATOR :: CUBLAS_FILL_MODE_LOWER=0
    ENUMERATOR :: CUBLAS_FILL_MODE_UPPER=1
    END ENUM !cublasFillMode_t

    ENUM, BIND(C) !:: cublasDiag    TYPE_t
    ENUMERATOR :: CUBLAS_DIAG_NON_UNIT=0
    ENUMERATOR :: CUBLAS_DIAG_UNIT=1
    END ENUM !cublasDiag    TYPE_t

    ENUM, BIND(C) !:: cublasSideMode_t
    ENUMERATOR :: CUBLAS_SIDE_LEFT =0
    ENUMERATOR :: CUBLAS_SIDE_RIGHT=1
    END ENUM !cublasSideMode_t

    ENUM, BIND(C) !:: cublasOperation_t
    ENUMERATOR :: CUBLAS_OP_N=0
    ENUMERATOR :: CUBLAS_OP_T=1
    ENUMERATOR :: CUBLAS_OP_C=2
    END ENUM !cublasOperation_t

    INTERFACE

        FUNCTION cublasInit() &
        &   BIND(C, NAME="cublasInit")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasInit
        END FUNCTION cublasInit

        FUNCTION cublasShutdown() &
        &   BIND(C, NAME="cublasShutdown")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasShutdown
        END FUNCTION cublasShutdown

        FUNCTION cublasCreate_v2(handle) &
        &   BIND(C, NAME="cublasCreate_v2")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasCreate_v2
            TYPE(C_PTR) :: handle
        END FUNCTION cublasCreate_v2

        FUNCTION cublasDestroy_v2(handle) &
        &   BIND(C, NAME="cublasDestroy_v2")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDestroy_v2
            TYPE(C_PTR), VALUE :: handle
        END FUNCTION cublasDestroy_v2

        FUNCTION cublasGetStream_v2(handle, stream) &
        &   BIND(C, NAME="cublasGetStream_v2")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasGetStream_v2
            TYPE(C_PTR), VALUE :: handle
            TYPE(C_PTR) :: stream
        END FUNCTION cublasGetStream_v2

        FUNCTION cublasSetStream_v2(handle, stream) &
        &   BIND(C, NAME="cublasSetStream_v2")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasSetStream_v2
            TYPE(C_PTR), VALUE :: handle
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cublasSetStream_v2

        FUNCTION cublasGetVector(n, elemSize, dx_src, incx, hy_dst, incy) &
        &   BIND(C, NAME="cublasGetVector")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasGetVector
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: dx_src
            INTEGER(C_INT), VALUE :: incx
            TYPE(C_PTR), VALUE :: hy_dst
            INTEGER(C_INT), VALUE :: incy
        END FUNCTION cublasGetVector

        FUNCTION cublasGetVectorAsync(n, elemSize, dx_src, incx, hy_dst, incy, stream) &
        &   BIND(C, NAME="cublasGetVectorAsync")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasGetVectorAsync
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: dx_src
            INTEGER(C_INT), VALUE :: incx
            TYPE(C_PTR), VALUE :: hy_dst
            INTEGER(C_INT), VALUE :: incy
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cublasGetVectorAsync

        FUNCTION cublasSetVector(n, elemSize, hx_src, incx, dy_dst, incy) &
        &   BIND(C, NAME="cublasSetVector")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasSetVector
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: hx_src
            INTEGER(C_INT), VALUE :: incx
            TYPE(C_PTR), VALUE :: dy_dst
            INTEGER(C_INT), VALUE :: incy
        END FUNCTION cublasSetVector

        FUNCTION cublasSetVectorAsync(n, elemSize, hx_src, incx, dy_dst, incy, stream) &
        &   BIND(C, NAME="cublasSetVectorAsync")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasSetVectorAsync
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: hx_src
            INTEGER(C_INT), VALUE :: incx
            TYPE(C_PTR), VALUE :: dy_dst
            INTEGER(C_INT), VALUE :: incy
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cublasSetVectorAsync

        FUNCTION cublasSetMatrix(rows, cols, elemSize, hA_src, lda, dB_dst, lddb) &
        &   BIND(C, NAME="cublasSetMatrix")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasSetMatrix
            INTEGER(C_INT), VALUE :: rows
            INTEGER(C_INT), VALUE :: cols
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: hA_src
            INTEGER(C_INT), VALUE :: lda
            TYPE(C_PTR), VALUE :: dB_dst
            INTEGER(C_INT), VALUE :: lddb
        END FUNCTION cublasSetMatrix

        FUNCTION cublasSetMatrixAsync(rows, cols, elemSize, hA_src, lda, dB_dst, lddb, stream) &
        &   BIND(C, NAME="cublasSetMatrixAsync")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasSetMatrixAsync
            INTEGER(C_INT), VALUE :: rows
            INTEGER(C_INT), VALUE :: cols
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: hA_src
            INTEGER(C_INT), VALUE :: lda
            TYPE(C_PTR), VALUE :: dB_dst
            INTEGER(C_INT), VALUE :: lddb
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cublasSetMatrixAsync

        FUNCTION cublasSetBatchMatrixAsync(rows, cols, batch, elemSize, hA_src, lda, dB_dst, lddb, stream) &
        &   BIND(C, NAME="cublasSetBatchMatrixAsync")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasSetBatchMatrixAsync
            INTEGER(C_INT), VALUE :: rows
            INTEGER(C_INT), VALUE :: cols
            INTEGER(C_INT), VALUE :: batch
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: hA_src
            INTEGER(C_INT), VALUE :: lda
            TYPE(C_PTR), VALUE :: dB_dst
            INTEGER(C_INT), VALUE :: lddb
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cublasSetBatchMatrixAsync

        FUNCTION cublasGetMatrix(rows, cols, elemSize, dA_src, ldda, hB_dst, ldb) &
        &   BIND(C, NAME="cublasGetMatrix")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasGetMatrix
            INTEGER(C_INT), VALUE :: rows
            INTEGER(C_INT), VALUE :: cols
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: dA_src
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: hB_dst
            INTEGER(C_INT), VALUE :: ldb
        END FUNCTION cublasGetMatrix

        FUNCTION cublasGetMatrixAsync(rows, cols, elemSize, dA_src, ldda, hB_dst, ldb, stream) &
        &   BIND(C, NAME="cublasGetMatrixAsync")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasGetMatrixAsync
            INTEGER(C_INT), VALUE :: rows
            INTEGER(C_INT), VALUE :: cols
            INTEGER(C_SIZE_T), VALUE :: elemSize
            TYPE(C_PTR), VALUE :: dA_src
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: hB_dst
            INTEGER(C_INT), VALUE :: ldb
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cublasGetMatrixAsync

        FUNCTION cublasDgetrfBatched(handle, n, dA, ldda, dP, dInfo, nbatch) &
        &   BIND(C, NAME="cublasDgetrfBatched")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDgetrfBatched
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: n
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dP
            TYPE(C_PTR), VALUE :: dInfo
            INTEGER(C_INT), VALUE :: nbatch
        END FUNCTION cublasDgetrfBatched

        FUNCTION cublasDgetriBatched(handle, n, dA, ldda, dP, dC, lddc, dInfo, nbatch) &
        &   BIND(C, NAME="cublasDgetriBatched")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDgetriBatched
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: n
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dP
            TYPE(C_PTR), VALUE :: dC
            INTEGER(C_INT), VALUE :: lddc
            TYPE(C_PTR), VALUE :: dInfo
            INTEGER(C_INT), VALUE :: nbatch
        END FUNCTION cublasDgetriBatched

        FUNCTION cublasDtrsmBatched(handle, side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb, nbatch) &
        &   BIND(C, NAME="cublasDtrsmBatched")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDtrsmBatched
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: side
            INTEGER(C_INT), VALUE :: uplo
            INTEGER(C_INT), VALUE :: trans
            INTEGER(C_INT), VALUE :: diag
            INTEGER(C_INT), VALUE :: m
            INTEGER(C_INT), VALUE :: n
            REAL(C_DOUBLE) :: alpha
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dB
            INTEGER(C_INT), VALUE :: lddb
            INTEGER(C_INT), VALUE :: nbatch
        END FUNCTION cublasDtrsmBatched

        FUNCTION cublasDgemmBatched(handle, transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc, nbatch) &
        &   BIND(C, NAME="cublasDgemmBatched")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDgemmBatched
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: transa
            INTEGER(C_INT), VALUE :: transb
            INTEGER(C_INT), VALUE :: m
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_INT), VALUE :: k
            REAL(C_DOUBLE) :: alpha
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dB
            INTEGER(C_INT), VALUE :: lddb
            REAL(C_DOUBLE) :: beta
            TYPE(C_PTR), VALUE :: dC
            INTEGER(C_INT), VALUE :: lddc
            INTEGER(C_INT), VALUE :: nbatch
        END FUNCTION cublasDgemmBatched

        FUNCTION cublasDtrsv(uplo, trans, diag, n, dA, ldda, dx, incx) &
        &   BIND(C, NAME="cublasDtrsv")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDtrsv
            CHARACTER(C_CHAR), VALUE :: uplo
            CHARACTER(C_CHAR), VALUE :: trans
            CHARACTER(C_CHAR), VALUE :: diag
            INTEGER(C_INT), VALUE :: n
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dx
            INTEGER(C_INT), VALUE :: incx
        END FUNCTION cublasDtrsv

        FUNCTION cublasDtrsv_v2(handle, uplo, trans, diag, n, dA, ldda, dx, incx) &
        &   BIND(C, NAME="cublasDtrsv_v2")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDtrsv_v2
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: uplo
            INTEGER(C_INT), VALUE :: trans
            INTEGER(C_INT), VALUE :: diag
            INTEGER(C_INT), VALUE :: n
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dx
            INTEGER(C_INT), VALUE :: incx
        END FUNCTION cublasDtrsv_v2

        FUNCTION cublasDtrsm(uplo, side, trans, diag, m, n, alpha, dA, ldda, dB, lddb) &
        &   BIND(C, NAME="cublasDtrsm")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDtrsm
            CHARACTER(C_CHAR), VALUE :: uplo
            CHARACTER(C_CHAR), VALUE :: side
            CHARACTER(C_CHAR), VALUE :: trans
            CHARACTER(C_CHAR), VALUE :: diag
            INTEGER(C_INT), VALUE :: m
            INTEGER(C_INT), VALUE :: n
            REAL(C_DOUBLE), VALUE :: alpha
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dB
            INTEGER(C_INT), VALUE :: lddb
        END FUNCTION cublasDtrsm

        FUNCTION cublasDtrsm_v2(handle, uplo, side, trans, diag, m, n, alpha, dA, ldda, dB, lddb) &
        &   BIND(C, NAME="cublasDtrsm_v2")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDtrsm_v2
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: uplo
            INTEGER(C_INT), VALUE :: side
            INTEGER(C_INT), VALUE :: trans
            INTEGER(C_INT), VALUE :: diag
            INTEGER(C_INT), VALUE :: m
            INTEGER(C_INT), VALUE :: n
            REAL(C_DOUBLE), VALUE :: alpha
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dB
            INTEGER(C_INT), VALUE :: lddb
        END FUNCTION cublasDtrsm_v2

        FUNCTION cublasDgemm(transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc) &
        &   BIND(C, NAME="cublasDgemm")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDgemm
            CHARACTER(C_CHAR), VALUE :: transa
            CHARACTER(C_CHAR), VALUE :: transb
            INTEGER(C_INT), VALUE :: m
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_INT), VALUE :: k
            REAL(C_DOUBLE), VALUE :: alpha
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dB
            INTEGER(C_INT), VALUE :: lddb
            REAL(C_DOUBLE), VALUE :: beta
            TYPE(C_PTR), VALUE :: dC
            INTEGER(C_INT), VALUE :: lddc
        END FUNCTION cublasDgemm

        FUNCTION cublasDgemm_v2(handle, transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc) &
        &   BIND(C, NAME="cublasDgemm_v2")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDgemm_v2
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: transa
            INTEGER(C_INT), VALUE :: transb
            INTEGER(C_INT), VALUE :: m
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_INT), VALUE :: k
            REAL(C_DOUBLE), VALUE :: alpha
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dB
            INTEGER(C_INT), VALUE :: lddb
            REAL(C_DOUBLE), VALUE :: beta
            TYPE(C_PTR), VALUE :: dC
            INTEGER(C_INT), VALUE :: lddc
        END FUNCTION cublasDgemm_v2

        FUNCTION cublasDgetrsBatched(handle, trans, n, nrhs, dA, ldda, dP, dB, lddb, hInfo, nbatch) &
        &   BIND(C, NAME="cublasDgetrsBatched")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: cublasDgetrsBatched
            TYPE(C_PTR), VALUE :: handle
            INTEGER(C_INT), VALUE :: trans
            INTEGER(C_INT), VALUE :: n
            INTEGER(C_INT), VALUE :: nrhs
            TYPE(C_PTR), VALUE :: dA
            INTEGER(C_INT), VALUE :: ldda
            TYPE(C_PTR), VALUE :: dP
            TYPE(C_PTR), VALUE :: dB
            INTEGER(C_INT), VALUE :: lddb
            TYPE(C_PTR), VALUE :: hInfo
            INTEGER(C_INT), VALUE :: nbatch
        END FUNCTION cublasDgetrsBatched

    END INTERFACE

END MODULE cublasf
