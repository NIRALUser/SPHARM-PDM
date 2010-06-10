
INCLUDE_DIRECTORIES(${SPHARMPDM_SOURCE_DIR}/Libraries/SparseLibMVIml)
LINK_DIRECTORIES(${SPHARMP_BINARY_DIR}/bin)

# here prepared only for ISO C++ compatibles compilers
 
# IBM xlC  v. 1.1
# CCCFLAGS        = -DCOMPLEX=complex
# LDFLAGS         = -lm -lcomplex

# Sun C++ 4.0.1
# CCCFLAGS        = -DMV_VECTOR_BOUNDS_CHECK -g -DCOMPLEX_OSTREAM -DCOMPLEX=complex
# LDFLAGS         =   -lm -lcomplex

# g++ v. 2.6.3
# CCCFLAGS                =   -O5 '-DCOMPLEX=std::complex<double>'
# LDFLAGS                 =   -lm


IF(WIN32)
ADD_DEFINITIONS("-DCOMPLEX=std::complex<double>")
ELSE(WIN32)
ADD_DEFINITIONS('-DCOMPLEX=std::complex<double>')
ENDIF(WIN32)

LINK_LIBRARIES(SparseMatrixLib)

