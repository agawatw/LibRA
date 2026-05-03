cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
project(LibRATests VERSION 1.0
          DESCRIPTION "LibRA tests"
          LANGUAGES CXX)

link_directories("${CMAKE_INSTALL_PREFIX}/lib")

# test_utils.h is shared across tests — expose this directory as an include path
set(TESTS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
find_package(GTest REQUIRED)

# Mirrors APP_LINK_LIBRARIES_CORE — needed by tests whose direct deps
# (casacore_asp, casacore_restore, libracore, librautils) link casacore PRIVATE.
# WCSLIB_LIBRARIES is injected by the outer superbuild.
find_package(GSL REQUIRED)
find_package(CFITSIO REQUIRED)
set(CASACORE_LIBS
  casa_derivedmscal casa_msfits casa_ms casa_images casa_coordinates
  casa_lattices casa_fits casa_measures casa_meas casa_mirlib
  casa_tables casa_scimath casa_scimath_f casa_casa
  pthread parafeed ${WCSLIB_LIBRARIES} CFITSIO::CFITSIO GSL::gsl GSL::gslcblas m c)

# Helper: register one test executable and install it
macro(add_libra_test target source)
  add_executable(${target} ${source})
  target_include_directories(${target} PRIVATE ${TESTS_INCLUDE_DIR})
  target_link_libraries(${target} GTest::gtest_main stdc++fs ${ARGN})
  add_test(NAME ${target} COMMAND ${target} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS ${target} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/tests")
endmacro()

# test_asp and test_restore depend on translation/ targets (casacore_asp,
# casacore_restore) which are defined after Asp/ and Restore/ in the subdirectory
# order — so their build targets live here where all targets are available.
add_libra_test(test_asp    ${CMAKE_CURRENT_SOURCE_DIR}/../Asp/test_asp.cc
  casacore_asp casacpp_synthesis ${CASACORE_LIBS})

add_libra_test(test_restore  ${CMAKE_CURRENT_SOURCE_DIR}/../Restore/test_restore.cc
  restore_lib casacore_restore ${CASACORE_LIBS})

add_libra_test(test_libracoretypes  test_libracoretypes.cc
  libracore ${CASACORE_LIBS})

add_libra_test(test_librautils  test_librautils.cc
  librautils ${CASACORE_LIBS})
