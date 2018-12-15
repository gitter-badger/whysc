cmake_minimum_required(VERSION 2.8)

project(WHYSC)

set(CMAKE_BUILD_TYPE Debug)
set(WHYSC_VERSION_MAJOR 3)
set(WHYSC_VERSION_MINOR 0)


# Set path for dependencies
# TODO: make the code more flexible to find packages 
set(CGAL_DIR  "/home/why/local")
set(P4EST_DIR "/home/why/local")

# MPI
find_package(MPI REQUIRED)
#if(MPI_C_FOUND)
#    include_directories(${MPI_C_INCLUDE_PATH})
#    link_directories(${MPI_C_LINK_FLAG})
#    set(HAVE_MPI_C True)
#else()
#    set(HAVE_MPI_C False)
#endif()
if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    link_directories(${MPI_CXX_LINK_FLAG})
    set(HAVE_MPI_CXX True)
    message(STATUS "Configure package with mpi cxx!")
else()
    set(HAVE_MPI_CXX False)
endif()


# p4est
set(P4EST_FOUND True)
if(P4EST_FOUND)
    include_directories(${P4EST_DIR}/include)
    link_directories(${P4EST_DIR}/lib)
    set(SC_LIBS sc)
    set(P4EST_LIBS p4est sc ${MPI_CXX_LIBRARIES} gomp lapack cblas f77blas atlas z m)
    set(HAVE_P4EST True)
    message(STATUS "Configure package with p4est!")
else()
    set(HAVE_CGAL False)
    message(STATUS "Configure package without p4est!")
endif()

# CGAL
set(CGAL_DIR "${CGAL_DIR}/lib/CGAL")
message(STATUS ${CGAL_DIR})
find_package(CGAL)
if(CGAL_FOUND)
    include_directories(${CGAL_INCLUDE_DIRS})
    link_directories(${CGAL_LIBRARIES_DIR})
    set(CGAL_LIBRARIES "-lCGAL -lCGAL_Core -lgmp -lmpfr")
    set(HAVE_CGAL True)
    message(STATUS "Configure package with cgal!")
else()
    set(HAVE_CGAL False)
    message(STATUS "Configure package without cgal!")
endif()

add_subdirectory(test)
add_subdirectory(example)

