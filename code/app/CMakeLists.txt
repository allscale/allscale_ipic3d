if(DEFINED ALLSCALECC)
	message("Using ${ALLSCALECC} as compiler")
	set(CMAKE_C_COMPILER ${ALLSCALECC})
	set(CMAKE_CXX_COMPILER ${ALLSCALECC})
endif()

add_module_library(app)

add_dependencies(app allscale)
target_include_directories(app PUBLIC ${ALLSCALE_API_INCLUDE_PATH})

glob_executables(app_exes src)
foreach(exe ${app_exes})
	add_module_executable(app ${exe})
endforeach(exe)

glob_tests(app_tests test)
foreach(test ${app_tests})
	add_module_unittest(app ${test})
endforeach(test)

add_definitions( -DPATH_TO_INPUTS="${PROJECT_SOURCE_DIR}/../inputs" )

configure_file(../../scripts/verify_output.sh.in ${PROJECT_BINARY_DIR}/verify_output.sh @ONLY IMMEDIATE)
