if(NOT TARGET allscale)
	include(ExternalProject)

	if(USE_ALLSCALECC)
		if(NOT EXISTS ${THIRD_PARTY_DIR})
			message(STATUS
				"=================================================================\n"
				"No third_party directory found, will set it up for you in 5 seconds:\n"
				"====================================================================\n"
			)
			execute_process(COMMAND ${CMAKE_COMMAND} -E sleep 5)
			execute_process(
				COMMAND bash ${PROJECT_SOURCE_DIR}/../scripts/dependencies/installer
				WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
			)
			execute_process(
				COMMAND bash ${PROJECT_SOURCE_DIR}/../scripts/dependencies/third_party_linker
				WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
			)
		endif()

		ExternalProject_Add(
			allscale
			GIT_REPOSITORY https://github.com/allscale/allscale_compiler
			GIT_TAG b3b7c10afee156826cf9fa5fd87ad395b5edb3c7
			CMAKE_COMMAND
				${CMAKE_COMMAND} -E env
				"INSIEME_LIBS_HOME=${THIRD_PARTY_DIR}"
				${CMAKE_COMMAND}
			CMAKE_ARGS
				${CMAKE_EXTERNALPROJECT_FORWARDS}
				-DINSIEME_C_BACKEND_COMPILER=${CMAKE_C_COMPILER}
				-DINSIEME_CXX_BACKEND_COMPILER=${CMAKE_CXX_COMPILER}
				-DINVOKED_AS_EXTERNAL_PROJECT=ON
				-DTHIRD_PARTY_DIR=${THIRD_PARTY_DIR}
			BUILD_COMMAND $(MAKE) compiler_allscalecc
			INSTALL_COMMAND ""
			EXCLUDE_FROM_ALL 1
			DOWNLOAD_NO_PROGRESS 1
		)
		ExternalProject_Get_Property(allscale source_dir binary_dir)

		set(ALLSCALECC ${binary_dir}/code/compiler/allscalecc)
		set(ALLSCALE_API_INCLUDE_PATH ${source_dir}/api/code/api/include ${source_dir}/api/code/utils/include)
	elseif(NOT OVERRIDE_ALLSCALE_API)
		ExternalProject_Add(
			allscale
			GIT_REPOSITORY https://github.com/allscale/allscale_api
			GIT_TAG c81015dc5339a9af834c5d760fe3ca18636bbf5c
			CONFIGURE_COMMAND ""
			BUILD_COMMAND ""
			INSTALL_COMMAND ""
			EXCLUDE_FROM_ALL 1
			DOWNLOAD_NO_PROGRESS 1
		)
		ExternalProject_Get_Property(allscale source_dir binary_dir)

		set(ALLSCALE_API_INCLUDE_PATH ${source_dir}/code/api/include ${source_dir}/code/utils/include)
	else()
		# dummy target
		add_custom_target(allscale)
	endif()

	if(OVERRIDE_ALLSCALECC)
		set(ALLSCALECC ${OVERRIDE_ALLSCALECC})
	endif()

	if(OVERRIDE_ALLSCALE_API)
		set(ALLSCALE_API_INCLUDE_PATH ${OVERRIDE_ALLSCALE_API}/code/api/include ${OVERRIDE_ALLSCALE_API}/code/utils/include)
	endif()
endif()
