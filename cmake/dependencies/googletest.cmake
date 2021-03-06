if(BUILD_TESTS AND NOT TARGET googletest)
	include(ExternalProject)

	if(MSVC)
		# disable checked iterators
		# this should be equal to the settings defined in msvc_settings.cmake
		set(_additional_c_flags /D_ITERATOR_DEBUG_LEVEL=0)
		set(_additional_cxx_flags /D_ITERATOR_DEBUG_LEVEL=0)
		# disable warning regarding deprecated tr1 namespace (introduced with MSVC 2017 15.5.1)
		set(_additional_c_flags "${_additional_c_flags} /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING=1")
		set(_additional_cxx_flags "${_additional_cxx_flags} /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING=1")

		set(_additional_flags -DCMAKE_C_FLAGS=${_additional_c_flags} -DCMAKE_CXX_FLAGS=${_additional_cxx_flags})
	endif()

	# gtest should be build with the same compiler as the project using it
	ExternalProject_Add(
		googletest
		URL http://insieme-compiler.org/ext_libs/gtest-1.8.0.tar.gz
		URL_HASH SHA256=58a6f4277ca2bc8565222b3bbd58a177609e9c488e8a72649359ba51450db7d8
		INSTALL_COMMAND ""
		CMAKE_ARGS ${CMAKE_EXTERNALPROJECT_FORWARDS} ${_additional_flags}
		DOWNLOAD_NO_PROGRESS 1
		EXCLUDE_FROM_ALL 1
	)
	ExternalProject_Get_Property(googletest source_dir binary_dir)

	# workaround https://itk.org/Bug/view.php?id=15052
	file(MAKE_DIRECTORY ${source_dir}/googletest/include)

	if(BUILD_SHARED_LIBS)
		set(_prefix ${CMAKE_SHARED_LIBRARY_PREFIX})
		set(_suffix ${CMAKE_SHARED_LIBRARY_SUFFIX})
		set(_linking SHARED)
	else()
		set(_prefix ${CMAKE_STATIC_LIBRARY_PREFIX})
		set(_suffix ${CMAKE_STATIC_LIBRARY_SUFFIX})
		set(_linking STATIC)
	endif()

	if(MSVC)
		set(GTEST_LIBRARY_PATH ${binary_dir}/googlemock/gtest/$(Configuration))
	else()
		set(GTEST_LIBRARY_PATH ${binary_dir}/googlemock/gtest/${LIBRARY_OUTPUT_DIRECTORY})
	endif()

	# import libgtest
	add_library(gtest ${_linking} IMPORTED)
	set(gtest ${_prefix}gtest${_suffix})
	set_target_properties(gtest PROPERTIES
		IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}/${gtest}

		# attach include path
		INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/googletest/include

		# attach pthread dependency
		INTERFACE_LINK_LIBRARIES Threads::Threads
	)
	add_dependencies(gtest googletest)

	# import libgtest_main
	add_library(gtest_main ${_linking} IMPORTED)
	set(gtest_main ${_prefix}gtest_main${_suffix})
	set_target_properties(gtest_main PROPERTIES IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}/${gtest_main})
	add_dependencies(gtest_main googletest)
endif()
