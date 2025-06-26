#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "SimTKcommon_recorder" for configuration "RelWithDebInfo"
set_property(TARGET SimTKcommon_recorder APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(SimTKcommon_recorder PROPERTIES
  IMPORTED_IMPLIB_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/SimTKcommon_recorder.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELWITHDEBINFO "lapack;blas;pthread"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/SimTKcommon_recorder.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS SimTKcommon_recorder )
list(APPEND _IMPORT_CHECK_FILES_FOR_SimTKcommon_recorder "${_IMPORT_PREFIX}/lib/SimTKcommon_recorder.lib" "${_IMPORT_PREFIX}/bin/SimTKcommon_recorder.dll" )

# Import target "SimTKmath_recorder" for configuration "RelWithDebInfo"
set_property(TARGET SimTKmath_recorder APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(SimTKmath_recorder PROPERTIES
  IMPORTED_IMPLIB_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/SimTKmath_recorder.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELWITHDEBINFO "SimTKcommon_recorder;lapack;blas;pthread"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/SimTKmath_recorder.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS SimTKmath_recorder )
list(APPEND _IMPORT_CHECK_FILES_FOR_SimTKmath_recorder "${_IMPORT_PREFIX}/lib/SimTKmath_recorder.lib" "${_IMPORT_PREFIX}/bin/SimTKmath_recorder.dll" )

# Import target "SimTKsimbody_recorder" for configuration "RelWithDebInfo"
set_property(TARGET SimTKsimbody_recorder APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(SimTKsimbody_recorder PROPERTIES
  IMPORTED_IMPLIB_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/SimTKsimbody_recorder.lib"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELWITHDEBINFO "SimTKmath_recorder;SimTKcommon_recorder;lapack;blas;pthread"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/SimTKsimbody_recorder.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS SimTKsimbody_recorder )
list(APPEND _IMPORT_CHECK_FILES_FOR_SimTKsimbody_recorder "${_IMPORT_PREFIX}/lib/SimTKsimbody_recorder.lib" "${_IMPORT_PREFIX}/bin/SimTKsimbody_recorder.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
