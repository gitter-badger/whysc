cmake_minimum_required(VERSION 3.1)

project(WHYSC)

set(CMAKE_BUILD_TYPE Debug)
set(WHYSC_VERSION_MAJOR 3)
set(WHYSC_VERSION_MINOR 0)
set(CMAKE_CXX_STANDARD 14)


include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)



# Set path for dependencies
# TODO: make the code more flexible to find packages 
set(CGAL_DIR  "/usr/local")
set(P4EST_DIR "/usr/local")
set(VTK_DIR "/usr/local/lib/cmake/vtk-8.1")

# MPI
find_package(MPI REQUIRED)
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
    set(P4EST_LIBS p4est sc ${MPI_CXX_LIBRARIES} gomp lapack z m)
    set(HAVE_P4EST True)
    message(STATUS "Configure package with p4est!")
else()
    set(HAVE_P4EST False)
    message(STATUS "Configure package without p4est!")
endif()

# CGAL
set(CGAL_DIR "${CGAL_DIR}/lib/CGAL")
message(STATUS ${CGAL_DIR})
find_package(CGAL)
if(CGAL_FOUND)
    include_directories(${CGAL_INCLUDE_DIRS})
    link_directories(${CGAL_LIBRARIES_DIR})
    set(CGAL_LIBS "-lCGAL -lCGAL_Core -lgmp -lmpfr")
    set(HAVE_CGAL True)
    message(STATUS "Configure package with cgal!")
else()
    set(HAVE_CGAL False)
    message(STATUS "Configure package without cgal!")
endif()

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
if(OpenGL_FOUND)
    include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
    set(HAVE_OpenGL True)
    message(STATUS "Configure package with OpenGL!")
    #${OPENGL_LIBRARIES} ${GLUT_LIBRARY}
else()
    set(HAVE_OpenGL False)
    message(STATUS "Configure package without OpenGL!")
endif()

# VTK
find_package(VTK REQUIRED)
if(VTK_FOUND)
    include(${VTK_USE_FILE})
    include_directories(${VTK_INCLUDE_DIRS})
    link_directories(${VTK_LIBRARIES_DIR})
    message(STATUS "Configure package with VTK!")
    set(HAVE_VTK True)
else()
    set(HAVE_VTK False)
    message(STATUS "Configure package without VTK!")
endif()
# VTK_LIBRARIES



set(WHYSC_LIBRARIES whysc)

add_subdirectory(test)
add_subdirectory(example)
#add_subdirectory(src)

