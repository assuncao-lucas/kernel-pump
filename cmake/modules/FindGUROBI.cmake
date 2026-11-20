# FindGUROBI.cmake -- this module finds the Gurobi Optimizer C library
#
# Input:  GUROBI_HOME (CMake variable or environment variable, e.g. C:/gurobi1400/win64 or /opt/gurobi1400/linux64)
# Output: GUROBI_FOUND, GUROBI_INCLUDE_DIR, GUROBI_LIBRARY, GUROBI_RUNTIME_LIBRARY (Windows DLL),
#         imported target Gurobi::Gurobi (shared library; the C++ wrapper is not used)
#
# Finding and linking Gurobi does not check out a license; only starting an environment does.

if (NOT GUROBI_HOME AND DEFINED ENV{GUROBI_HOME})
	set(GUROBI_HOME "$ENV{GUROBI_HOME}")
endif()
file(TO_CMAKE_PATH "${GUROBI_HOME}" GUROBI_HOME)

find_path(GUROBI_INCLUDE_DIR "gurobi_c.h" PATHS "${GUROBI_HOME}/include" NO_DEFAULT_PATH)

# The library name carries the version (gurobi140.lib, libgurobi140.so): glob for it, newest first,
# excluding the C++ wrapper (gurobi_c++*) and the .NET/light variants.
if (WIN32)
	file(GLOB GUROBI_LIBRARY_CANDIDATES "${GUROBI_HOME}/lib/gurobi[0-9]*.lib")
	file(GLOB GUROBI_RUNTIME_CANDIDATES "${GUROBI_HOME}/bin/gurobi[0-9]*.dll")
	list(FILTER GUROBI_RUNTIME_CANDIDATES EXCLUDE REGEX "_light|NET")
else()
	file(GLOB GUROBI_LIBRARY_CANDIDATES "${GUROBI_HOME}/lib/libgurobi[0-9]*${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()
list(FILTER GUROBI_LIBRARY_CANDIDATES EXCLUDE REGEX "gurobi_c\\+\\+|netstandard|_light")
if (GUROBI_LIBRARY_CANDIDATES)
	list(SORT GUROBI_LIBRARY_CANDIDATES ORDER DESCENDING)
	list(GET GUROBI_LIBRARY_CANDIDATES 0 _gurobi_lib)
	set(GUROBI_LIBRARY "${_gurobi_lib}" CACHE FILEPATH "Gurobi C library (import library on Windows)")
endif()
if (GUROBI_RUNTIME_CANDIDATES)
	list(SORT GUROBI_RUNTIME_CANDIDATES ORDER DESCENDING)
	list(GET GUROBI_RUNTIME_CANDIDATES 0 _gurobi_dll)
	set(GUROBI_RUNTIME_LIBRARY "${_gurobi_dll}" CACHE FILEPATH "Gurobi runtime DLL")
endif()

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_RUNTIME_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY GUROBI_INCLUDE_DIR)

if (GUROBI_FOUND AND NOT TARGET Gurobi::Gurobi)
	add_library(Gurobi::Gurobi SHARED IMPORTED)
	set_target_properties(Gurobi::Gurobi PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIR}"
	)
	if (WIN32)
		set_target_properties(Gurobi::Gurobi PROPERTIES IMPORTED_IMPLIB "${GUROBI_LIBRARY}")
		if (GUROBI_RUNTIME_LIBRARY)
			set_target_properties(Gurobi::Gurobi PROPERTIES IMPORTED_LOCATION "${GUROBI_RUNTIME_LIBRARY}")
		endif()
	else()
		set_target_properties(Gurobi::Gurobi PROPERTIES IMPORTED_LOCATION "${GUROBI_LIBRARY}")
	endif()
endif()
