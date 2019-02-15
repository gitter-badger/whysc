
#file(GLOB_RECURSE SOURCE ./*.cpp)
set(LIBNAME whysc)

add_library(${LIBNAME} SHARED ${SOURCE})
add_library(${LIBNAME}_static STATIC ${SOURCE})
set_target_properties(${LIBNAME} PROPERTIES VERSION 0.0 SOVERTSION 0)
set_target_properties(${LIBNAME}_static PROPERTIES OUTPUT_NAME ${LIBNAME})
install(TARGETS  ${LIBNAME} ${LIBNAME}_static
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
