cmake_minimum_required(VERSION 3.16)
#set (CMAKE_CXX_FLAGS "-lstdc++fs -std=c++17")
set (CMAKE_CXX_FLAGS " -std=c++17")

project(display_gp)
#https://stackoverflow.com/questions/8774593/cmake-link-to-external-library
#add_library(bayesopt_tut_libs main.cpp)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
find_library(BAYESOPT_LIB bayesopt REQUIRED)
find_library(NLOPT_LIB nlopt REQUIRED)

find_library(MLPACK_LIB
  NAMES mlpack
  PATHS "$ENV{ProgramFiles}/mlpack/" /usr/lib64/ /usr/lib/ /usr/local/lib64/ /usr/local/
  REQUIRED
)
find_package(Boost REQUIRED)
find_package(CGAL REQUIRED)
message(OpenMPlibrarypath="${OpenMP_CXX_LIBRARIES}")
message(bayeslibrarypath="${BAYESOPT_LIB}")
message(nloptlibrarypath="${NLOPT_LIB}")
message(cgullibrarypath="${CGAL_DIR}" )
message(mlpacklibrary= "${MLPACK_LIB}")

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
find_package(GLUT REQUIRED)
message(GLUT_library "${GLUT_LIBRARY}")
find_package(OpenGL REQUIRED)
if(OPENGL_FOUND AND GLUT_FOUND)



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
    utils/helper.hpp
   
)   

set(SRCS
    src/display2dgp.cpp
    src/tumorModel.cpp
    src/meanShift.cpp
    src/contour.cpp
    src/helper.cpp
    src/parameters.cpp
    )
set(HEADERS
    include/display2dgp.hpp
    include/meanShift.hpp
    include/contour.hpp
    include/tumorModel.hpp
    include/evaluation.hpp
    include/helper.hpp
    include/parameters.hpp
    )



include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include )

add_library(utils_lib ${UTILS_SOURCES})
target_include_directories(utils_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utils)


add_library(headers ${SRCS} )

target_link_libraries(headers ${BAYESOPT_LIB} ${NLOPT_LIB} matplotpp ${MLPACK_LIB}
                    ${GLUT_LIBRARY} ${OPENGL_LIBRARY}  alglib utils_lib)

add_executable(display_gp src/main.cpp )
add_executable(compute_contour src/contour_main.cpp)
add_executable(evaluate src/evaluation.cpp)


target_link_libraries(compute_contour alglib  ${OpenMP_CXX_LIBRARIES} ${BAYESOPT_LIB} ${NLOPT_LIB} ${MLPACK_LIB} CGAL::CGAL ${Boost_LIBRARIES} headers utils_lib  )

target_link_libraries(display_gp ${OpenMP_CXX_LIBRARIES} ${BAYESOPT_LIB} ${NLOPT_LIB} matplotpp alglib ${MLPACK_LIB} 
${GLUT_LIBRARY} ${OPENGL_LIBRARY} headers utils_lib)


target_link_libraries(evaluate alglib  ${OpenMP_CXX_LIBRARIES} ${BAYESOPT_LIB} ${NLOPT_LIB} ${MLPACK_LIB} CGAL::CGAL ${Boost_LIBRARIES} headers utils_lib  )
target_compile_options(evaluate PRIVATE -fopenmp)
endif()
