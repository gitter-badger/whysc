if(HAVE_P4EST AND HAVE_CGAL)
    add_executable(adaptive_interface_fitted_mesh_generation 
        adaptive_interface_fitted_mesh_generation.cpp)
    target_link_libraries(adaptive_interface_fitted_mesh_generation 
        ${P4EST_LIBS})
endif()
