cmake_minimum_required(VERSION 3.16)
project(display_gp)
#https://stackoverflow.com/questions/8774593/cmake-link-to-external-library
#add_library(bayesopt_tut_libs main.cpp)
find_package(OpenMP)

find_library(BAYESOPT_LIB bayesopt)
find_library(NLOPT_LIB nlopt)
find_library(MLPACK_LIB mlpack)
find_package(Boost)
find_package(CGAL)
message(bayeslibrarypath="${BAYESOPT_LIB}")
message(nloptlibrarypath="${NLOPT_LIB}")
message(cgullibrarypath="${CGAL_DIR}" )
enable_language(C)  #For OpenGL and bo_display demo


set ( ALGLIB_HEADER_FILES
libs/alglibinternal.h
libs/alglibmisc.h
libs/ap.h
libs/dataanalysis.h
libs/diffequations.h
libs/fasttransforms.h
libs/integration.h
libs/interpolation.h
libs/kernels_avx2.h
libs/kernels_fma.h
libs/kernels_sse2.h
libs/linalg.h
libs/optimization.h
libs/solvers.h
libs/specialfunctions.h
libs/statistics.h
libs/stdafx.h
)
set( ALGLIB_SOURCE_FILES
libs/alglibinternal.cpp
libs/alglibmisc.cpp
libs/ap.cpp
libs/dataanalysis.cpp
libs/diffequations.cpp
libs/fasttransforms.cpp
libs/integration.cpp
libs/interpolation.cpp
libs/kernels_avx2.cpp
libs/kernels_fma.cpp
libs/kernels_sse2.cpp
libs/linalg.cpp
libs/optimization.cpp
libs/solvers.cpp
libs/specialfunctions.cpp
libs/statistics.cpp
)


set(LIBS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libs")
message(libs_path=${LIBS_PATH})
message(ALGLIB_HEADER_FILES=${ALGLIB_HEADER_FILES})
add_library(alglib ${ALGLIB_SOURCE_FILES} ${ALGLIB_HEADER_FILES} )
target_include_directories(alglib PUBLIC ${LIBS_PATH})
#Display test
find_package(GLUT)
find_package(OpenGL)
if(OPENGL_FOUND AND GLUT_FOUND)

  SET(DISPLAY_SRC     
  /home/raphael/robolab/displaygp/utils/displaygp.cpp
    )
    message(DISPLAY_SRC="${DISPLAY_SRC}")

  if(GLUT_LIBRARY MATCHES freeglut_static.lib)
      add_definitions(-DFREEGLUT_STATIC)
endif()


  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/matplotpp
		      ${GLUT_INCLUDE_DIRS} 
		      ${GLUT_INCLUDE_DIR} 
		      ${OpenGL_INCLUDE_DIRS}
          ${Boost_INCLUDE_DIRS})			  

  link_directories(${GLUT_LIBRARY_DIRS}
		   ${OpenGL_LIBRARY_DIRS})

  add_definitions(${GLUT_DEFINITIONS}
		  ${OpenGL_DEFINITIONS})

  ADD_LIBRARY(matplotpp STATIC 
    matplotpp/matplotpp.cc 
    matplotpp/gl2ps.c
    )

  TARGET_LINK_LIBRARIES(matplotpp
    ${GLUT_LIBRARY} ${OPENGL_LIBRARY})

set(UTILS_SOURCES
    utils/displaygp.cpp
    utils/fileparser.cpp
    utils/param_loader.cpp
    utils/parser.cpp
    utils/ublas_extra.cpp
   )
set(UTILS_HEADERS
    utils/boundingbox.hpp
    utils/displaygp.hpp
    utils/fileparser.hpp
    utils/gridsampling.hpp
    utils/indexvector.hpp
    utils/lhs.hpp
    utils/log.hpp
    utils/param_loader.hpp
    utils/parser.hpp
    utils/testfunctions.hpp
    utils/ublas_cholesky.hpp
    utils/ublas_cholesky.hpp
    utils/ublas_elementwise.hpp
    utils/ublas_elementwise.hpp
    utils/ublas_extra.hpp
    utils/ublas_trace.hpp
    utils/prob_distribution.hpp
    utils/randgen.hpp
   
)   
set(SRCS
    src/display2dgp.cpp
    src/tumorModel.cpp
    src/fileParser.cpp
    src/meanShift.cpp
    src/contour.cpp
    )
set(HEADERS
    include/display2dgp.hpp
    include/meanShift.hpp
    include/fileParser.hpp
    include/contour.hpp
    include/tumorModel.hpp
    include/evaluation.hpp)


#create cgul lib
#find_library(GMP_LIBRARY NAMES gmp)
#find_library(GMPXX_LIBRARY NAMES gmpxx)
#add_library(cgul_lib INTERFACE)
#target_include_directories(cgul_lib INTERFACE ${CGAL_DIR})

#target_link_libraries(cgul_lib INTERFACE ${Boost_LIBRARIES} ${GMP_LIBRARY} ${GMPXX_LIBRARY})
#create library from additional sources and include headers
add_library(headers ${SRCS} ${HEADERS} )
target_link_libraries(headers ${BAYESOPT_LIB} ${NLOPT_LIB} matplotpp
                    ${GLUT_LIBRARY} ${OPENGL_LIBRARY} utils_lib alglib)

add_library(utils_lib ${UTILS_SOURCES} ${UTILS_HEADERS} )
target_include_directories(utils_lib PUBLIC /home/raphael/robolab/displaygp/utils)


add_executable(twod_display_gp src/twod_main.cpp )
add_executable(oned_display_gp src/oned_main.cpp )
add_executable(display_gp src/main.cpp )
#add_executable(compute_mean_shift src/meanShift.cpp )
add_executable(compute_contour src/contour_main.cpp)
add_executable(test_display_gp src/test.cpp )
add_executable(evaluate src/evaluation.cpp)
#target_link_libraries(bayesopt_tut bayesopt_tut_libs)
include_directories(include ${BAYESOPT_LIB_INCLUDE_DIRS}${NLOPT_LIB_INCLUDE_DIRS})

target_link_libraries(twod_display_gp ${BAYESOPT_LIB} ${NLOPT_LIB} matplotpp
${GLUT_LIBRARY} ${OPENGL_LIBRARY} utils_lib)


target_link_libraries(oned_display_gp ${BAYESOPT_LIB} ${NLOPT_LIB} matplotpp
${GLUT_LIBRARY} ${OPENGL_LIBRARY} utils_lib)


#target_link_libraries(compute_mean_shift ${BAYESOPT_LIB} ${NLOPT_LIB} ${MLPACK_LIB}   OpenMP::OpenMP_CXX
#${GLUT_LIBRARY} ${OPENGL_LIBRARY}  headers)

target_link_libraries(compute_contour alglib headers ${BAYESOPT_LIB} ${NLOPT_LIB} ${MLPACK_LIB} utils_lib  OpenMP::OpenMP_CXX)

target_link_libraries(display_gp ${BAYESOPT_LIB} ${NLOPT_LIB} matplotpp
${GLUT_LIBRARY} ${OPENGL_LIBRARY} utils_lib headers)


target_link_libraries(test_display_gp ${BAYESOPT_LIB} ${NLOPT_LIB} matplotpp
${GLUT_LIBRARY} ${OPENGL_LIBRARY} utils_lib)


target_link_libraries(evaluate alglib headers ${BAYESOPT_LIB} ${NLOPT_LIB} ${MLPACK_LIB} CGAL::CGAL ${Boost_LIBRARIES} utils_lib  )

endif()
