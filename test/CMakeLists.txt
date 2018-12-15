
if(HAVE_P4EST AND HAVE_CGAL)
    add_executable(hello_p4est hello_p4est.cpp)
    target_link_libraries(hello_p4est ${P4EST_LIBS})

    add_executable(test_interface_mesh_alg_2 test_interface_mesh_alg_2.cpp)
    target_link_libraries(test_interface_mesh_alg_2 ${P4EST_LIBS}
        ${CGAL_LIBS})
endif()
