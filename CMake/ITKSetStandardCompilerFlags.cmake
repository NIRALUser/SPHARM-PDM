# trict-null-sentinel Check the set of common warning flags supported by C and C++ compilers
# check_compiler_warning_flags(<c_flags_var> <cxx_flags_var>)
#  <c_flags_var> - variable to store valid C warning flags
#  <cxx_flags_var> - variable to store valid CXX warning flags
# This internally calls the check_c_compiler_flag and check_cxx_compiler_flag macros.


# To create a portable build system, it is best to not
# test for platforms, but to test for features.
#
# Instead of testing "if Windows then do this", test for
# "if the -Wno-invalid-offsetof flag works then use it".
#
# Typical use of this module is:
#
#  include(CheckCompilerWarningFlags)
#  check_compiler_warning_flags(C_WARNING_FLAGS CXX_WARNING_FLAGS)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_WARNING_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS}")


include(ITK_CheckCCompilerFlag)
include(ITK_CheckCXXCompilerFlag)

function(check_c_compiler_warning_flags c_flag_var)
  set(local_c_flags "")
  set(flag_list "${ARGN}")
  foreach(flag IN LISTS flag_list)
    ITK_CHECK_C_COMPILER_FLAG(${flag} C_HAS_WARNING${flag})
    if(${C_HAS_WARNING${flag}})
      set(local_c_flags "${local_c_flags} ${flag}")
    endif()
  endforeach()
  set(${c_flag_var} "${local_c_flags}" PARENT_SCOPE)
endfunction()


function(check_cxx_compiler_warning_flags cxx_flag_var)
  set(local_cxx_flags "")
  set(flag_list "${ARGN}")
  foreach(flag IN LISTS flag_list)
    ITK_CHECK_CXX_COMPILER_FLAG(${flag} CXX_HAS_WARNING${flag})
    if(${CXX_HAS_WARNING${flag}})
      set(local_cxx_flags "${local_cxx_flags} ${flag}")
    endif()
  endforeach()
  set(${cxx_flag_var} "${local_cxx_flags}" PARENT_SCOPE)
endfunction()


function(check_compiler_warning_flags c_warning_flags_var cxx_warning_flags_var)
  set(${c_warning_flags_var} "" PARENT_SCOPE)
  set(${cxx_warning_flags_var} "" PARENT_SCOPE)

  # Check this list on C compiler only
  set(c_flags
    -Wno-uninitialized
    -Wno-unused-parameter
  )

  ## On windows, the most verbose compiler options
  ## is reporting 1000's of wanings in windows
  ## header files, for now, limit the number of
  ## warnings to level 3
  if( WIN32 )
    set(VerboseWarningsFlag -W3 )
    ## A better solution would be to use -Wall,
    ## and then disable warnings one by one
    ## set(VerboseWarningsFlag -Wall -wd4820 -wd4682 )
  else()
    ## with Intel compiler, the -Wall compiler options
    ## is reporting 1000's of remarks of trivial items
    ## that will only slow day-to-day operations
    ## specify -w2 to restrict to only warnings and errors
    if (${CMAKE_C_COMPILER} MATCHES "icc.*$")
      set(USING_INTEL_ICC_COMPILER TRUE)
    endif()
    if (${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
      set(USING_INTEL_ICC_COMPILER TRUE)
    endif()
    if(USING_INTEL_ICC_COMPILER)
      # NOTE -w2 is close to gcc's -Wall warning level, -w5 is intels -Wall warning level, and it is too verbose.
      set(VerboseWarningsFlag -w2 -wd1268 -wd981 -wd383 -wd1418 -wd1419 -wd2259 -wd1572 -wd424 )
      #-wd424  #Needed for Intel compilers with remarki  #424: extra ";" ignored
      #-wd383  #Needed for Intel compilers with remark   #383: value copied to temporary, reference to temporary used
      #-wd981  #Needed for Intel compilers with remark   #981: operands are evaluated in unspecified order
      #-wd1418 #Needed for Intel compilers with remark  #1418: external function definition with no prior declaration
      #-wd1419 #Needed for Intel compilers with remark  #1419: external declaration in primary source file
      #-wd1572 #Needed for Intel compilers with remark  #1572: floating-point equality and inequality comparisons are unreliable
      #-wd2259 #Needed for Intel compilers with remark  #2259: non-pointer conversion from "itk::SizeValueType={unsigned long}" to "double" may lose significant bits
      #-wd1268 #Needed for Intel compliers with warning #1268: support for exported templates is disabled
    else()
      set(VerboseWarningsFlag -Wall )
    endif ()
  endif()

  # Check this list on both C and C++ compilers
  set(c_and_cxx_flags
    ${VerboseWarningsFlag}
# -Wno-long-double is only applicable with gcc on apple, and is generating
# warnings now on all platforms.
#    -Wno-long-double        #Needed on APPLE
    -Wcast-align
    -Wdisabled-optimization
    -Wextra
    -Wformat=2
    -Winvalid-pch
    -Wno-format-nonliteral
    -Wpointer-arith
    -Wshadow
    -Wunused
    -Wwrite-strings
    -funit-at-a-time
    -Wno-strict-overflow
  )

  # Check this list on C++ compiler only
  set(cxx_flags
    -Wno-deprecated
    -Wno-invalid-offsetof
    -Woverloaded-virtual
    -Wstrict-null-sentinel
  )
##-Wno-c++0x-static-nonintegral-init
    ## Clang compiler likes to warn about this feature that is technically only in
    ## c++0x, but works on many compilers, and if it fails, then alternate methods are used

  check_c_compiler_warning_flags(CMAKE_C_WARNING_FLAGS ${c_flags} ${c_and_cxx_flags})
  check_cxx_compiler_warning_flags(CMAKE_CXX_WARNING_FLAGS ${c_and_cxx_flags} ${cxx_flags})

  set(${c_warning_flags_var} "${CMAKE_C_WARNING_FLAGS}" PARENT_SCOPE)
  set(${cxx_warning_flags_var} "${CMAKE_CXX_WARNING_FLAGS}" PARENT_SCOPE)
endfunction()


macro(check_compiler_platform_flags)
  # On Visual Studio 8 MS deprecated C. This removes all 1.276E1265 security
  # warnings
  if(WIN32)
       if(NOT MINGW)
         if(NOT ITK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
           add_definitions(
             -D_CRT_FAR_MAPPINGS_NO_DEPRECATE
             -D_CRT_IS_WCTYPE_NO_DEPRECATE
             -D_CRT_MANAGED_FP_NO_DEPRECATE
             -D_CRT_NONSTDC_NO_DEPRECATE
             -D_CRT_SECURE_NO_DEPRECATE
             -D_CRT_SECURE_NO_DEPRECATE_GLOBALS
             -D_CRT_SETERRORMODE_BEEP_SLEEP_NO_DEPRECATE
             -D_CRT_TIME_FUNCTIONS_NO_DEPRECATE
             -D_CRT_VCCLRIT_NO_DEPRECATE
             -D_SCL_SECURE_NO_DEPRECATE
             )
         endif()
         # With MS compilers on Win64, we need the /bigobj switch, else generated
         # code results in objects with number of sections exceeding object file
         # format.
         # see http://msdn.microsoft.com/en-us/library/ms173499.aspx
         if(MSVC_VERSION GREATER 1310)
           set(ITK_REQUIRED_CXX_FLAGS "${ITK_REQUIRED_CXX_FLAGS} /bigobj")
         endif()
       endif()
  endif()
endmacro()#End the platform check function

#-----------------------------------------------------------------------------
#Check the set of warning flags the compiler supports
check_compiler_warning_flags(C_WARNING_FLAGS CXX_WARNING_FLAGS)

# Append ITK warnings to the CMake flags.
# We do not set them in ITK_REQUIRED FLAGS because all project which
# use ITK don't require these flags .
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_WARNING_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS}")

#-----------------------------------------------------------------------------
#Check the set of platform flags the compiler supports
check_compiler_platform_flags()
