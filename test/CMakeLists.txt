
if(HAVE_P4EST AND HAVE_CGAL)
    add_executable(hello_p4est hello_p4est.cpp)
    target_link_libraries(hello_p4est ${P4EST_LIBS})

    add_executable(test_interface_mesh_alg_2 test_interface_mesh_alg_2.cpp)
    target_link_libraries(test_interface_mesh_alg_2 ${P4EST_LIBS} ${CGAL_LIBS} ${VTK_LIBRARIES})

    add_executable(test_interface_mesh_generator_3d test_interface_mesh_generator_3d.cpp)
    target_link_libraries(test_interface_mesh_generator_3d ${P4EST_LIBS} ${VTK_LIBRARIES})
endif()


if(HAVE_VTK)
    add_executable(test_vtk test_vtk.cpp)
    target_link_libraries(test_vtk ${VTK_LIBRARIES})     

    add_executable(test_polydata test_polydata.cpp)
    target_link_libraries(test_polydata ${VTK_LIBRARIES})     

    add_executable(test_pwriter test_pwriter.cpp)
    target_link_libraries(test_pwriter ${VTK_LIBRARIES})     

    add_executable(test_vtkpolyhedron test_vtkpolyhedron.cpp)
    target_link_libraries(test_vtkpolyhedron ${VTK_LIBRARIES})     

    add_executable(test_vtkhexahedron test_vtkhexahedron.cpp)
    target_link_libraries(test_vtkhexahedron ${VTK_LIBRARIES})     
endif()

if(HAVE_OpenGL)
    add_executable(test_opengl test_opengl.cpp)
    target_link_libraries(test_opengl ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})     
endif()

add_executable(test_array test_array.cpp)
add_executable(test_edges test_edges.cpp)
add_executable(test_geo test_geo.cpp)

add_executable(test_mesh3d test_mesh3d.cpp)
#add_executable(test_algebra test_algebra.cpp)
#target_link_libraries(test_mesh3d ${WHYSC_LIBRARIES})     
