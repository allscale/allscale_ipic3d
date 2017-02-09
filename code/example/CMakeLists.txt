if(USE_ALLSCALECC)
	set(CMAKE_C_COMPILER ${ALLSCALECC})
	set(CMAKE_CXX_COMPILER ${ALLSCALECC})
endif()

add_module_library(example)

add_dependencies(example allscale)
target_include_directories(example PUBLIC ${ALLSCALE_API_INCLUDE_PATH})

glob_executables(example_exes src)
foreach(exe ${example_exes})
	add_module_executable(example ${exe})
endforeach(exe)

glob_tests(example_tests test)
foreach(test ${example_tests})
	add_module_unittest(example ${test})
endforeach(test)
