# Apple Silicon GNU Fortran does not support -mcmodel=large.
# Override generic GNU defaults for mac builds.
string(REPLACE " -mcmodel=large" " -mcmodel=small" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
string(REPLACE " -mcmodel=large" " -mcmodel=small" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

# Point E3SM CMake netcdf discovery at Homebrew defaults on macOS.
if (NOT DEFINED ENV{NETCDF_C_PATH})
	set(ENV{NETCDF_C_PATH} "/opt/homebrew/opt/netcdf")
endif()
if (NOT DEFINED ENV{NETCDF_FORTRAN_PATH})
	set(ENV{NETCDF_FORTRAN_PATH} "/opt/homebrew/opt/netcdf-fortran")
endif()
