# macro for building PREQx primitive equations executable using Kokkos for on-node parallelism

MACRO(PREQX_KOKKOS_SETUP)

  SET (KOKKOS_TARGET TRUE)

  INCLUDE_DIRECTORIES ( SYSTEM ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
  LINK_DIRECTORIES    ( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

  SET(PREQX_SHARE_DIR       ${HOMME_SOURCE_DIR}/src/preqx/share)
  SET(TARGET_DIR            ${HOMME_SOURCE_DIR}/src/preqx_kokkos)
  SET(SRC_DIR               ${HOMME_SOURCE_DIR}/src)
  SET(SRC_SHARE_DIR         ${HOMME_SOURCE_DIR}/src/share)
  SET(TRILINOS_ZOLTAN_DIR   ${HOMME_SOURCE_DIR}/src/zoltan)
  SET(TEST_SRC_DIR          ${HOMME_SOURCE_DIR}/src/test_src)

  SET(UTILS_TIMING_SRC_DIR  ${HOMME_SOURCE_DIR}/utils/cime/CIME/non_py/src/timing)
  SET(UTILS_TIMING_DIR      ${HOMME_BINARY_DIR}/utils/cime/CIME/non_py/src/timing)

  # Make INCLUDE_DIRS global so the tests can access it
  SET (EXEC_LIB_INCLUDE_DIRS
    ${TARGET_DIR}/cxx
    ${UTILS_TIMING_DIR}
    ${UTILS_TIMING_SRC_DIR}
    ${PREQX_SHARE_DIR}
    ${SRC_SHARE_DIR}/cxx
    ${HOMME_BINARY_DIR}/src/share/cxx
  )

  SET (EXEC_INCLUDE_DIRS
    ${EXEC_LIB_INCLUDE_DIRS}
    ${PIO_INCLUDE_DIRS}
  )

  SET (TARGET_F90
    ${TARGET_DIR}/bndry_mod.F90
    ${TARGET_DIR}/derivative_mod.F90
    ${TARGET_DIR}/edge_mod.F90
    ${TARGET_DIR}/element_state.F90
    ${TARGET_DIR}/model_init_mod.F90
    ${TARGET_DIR}/prim_advection_mod.F90
    ${TARGET_DIR}/prim_driver_mod.F90
    ${TARGET_DIR}/preqx_f2c_mod.F90
    ${TARGET_DIR}/viscosity_mod.F90
    ${PREQX_SHARE_DIR}/element_ops.F90
    ${PREQX_SHARE_DIR}/prim_advance_mod.F90
    ${PREQX_SHARE_DIR}/prim_state_mod.F90
    ${PREQX_SHARE_DIR}/viscosity_preqx_base.F90
    ${PREQX_SHARE_DIR}/vertremap_mod.F90
  )

  SET (SRC_SHARE_F90
    ${SRC_SHARE_DIR}/bndry_mod_base.F90
    ${SRC_SHARE_DIR}/cg_mod.F90
    ${SRC_SHARE_DIR}/control_mod.F90
    ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
    ${SRC_SHARE_DIR}/cube_mod.F90
    ${SRC_SHARE_DIR}/derivative_mod_base.F90
    ${SRC_SHARE_DIR}/dimensions_mod.F90
    ${SRC_SHARE_DIR}/dof_mod.F90
    ${SRC_SHARE_DIR}/domain_mod.F90
    ${SRC_SHARE_DIR}/edge_mod_base.F90
    ${SRC_SHARE_DIR}/edgetype_mod.F90
    ${SRC_SHARE_DIR}/element_mod.F90
    ${SRC_SHARE_DIR}/gllfvremap_mod.F90
    ${SRC_SHARE_DIR}/gllfvremap_util_mod.F90
    ${SRC_SHARE_DIR}/global_norms_mod.F90
    ${SRC_SHARE_DIR}/gridgraph_mod.F90
    ${SRC_SHARE_DIR}/hybrid_mod.F90
    ${SRC_SHARE_DIR}/hybvcoord_mod.F90
    ${SRC_SHARE_DIR}/interpolate_mod.F90
    ${SRC_SHARE_DIR}/kinds.F90
    ${SRC_SHARE_DIR}/ll_mod.F90
    ${SRC_SHARE_DIR}/mass_matrix_mod.F90
    ${SRC_SHARE_DIR}/mesh_mod.F90
    ${SRC_SHARE_DIR}/sort_mod.F90
    ${SRC_SHARE_DIR}/metagraph_mod.F90
    ${SRC_SHARE_DIR}/metis_mod.F90
    ${SRC_SHARE_DIR}/namelist_mod.F90
    ${SRC_SHARE_DIR}/parallel_mod.F90
    ${SRC_SHARE_DIR}/params_mod.F90
    ${SRC_SHARE_DIR}/physical_constants.F90
    ${SRC_SHARE_DIR}/physics_mod.F90
    ${SRC_SHARE_DIR}/prim_advection_base.F90
    ${SRC_SHARE_DIR}/prim_derived_type_mod.F90
    ${SRC_SHARE_DIR}/prim_driver_base.F90
    ${SRC_SHARE_DIR}/prim_implicit_mod.F90
    ${SRC_SHARE_DIR}/prim_si_mod.F90
    ${SRC_SHARE_DIR}/quadrature_mod.F90
    ${SRC_SHARE_DIR}/reduction_mod.F90
    ${SRC_SHARE_DIR}/schedtype_mod.F90
    ${SRC_SHARE_DIR}/schedule_mod.F90
    ${SRC_SHARE_DIR}/solver_mod.F90
    ${SRC_SHARE_DIR}/spacecurve_mod.F90
    ${SRC_SHARE_DIR}/thread_mod.F90
    ${SRC_SHARE_DIR}/time_mod.F90
    ${SRC_SHARE_DIR}/unit_tests_mod.F90
    ${SRC_SHARE_DIR}/vertremap_base.F90
    ${SRC_SHARE_DIR}/viscosity_base.F90
    ${SRC_SHARE_DIR}/zoltan_mod.F90
    ${SRC_SHARE_DIR}/sl_advection.F90
    ${SRC_SHARE_DIR}/compose_mod.F90
    ${SRC_SHARE_DIR}/compose_test_mod.F90
    ${SRC_SHARE_DIR}/scalable_grid_init_mod.F90
    ${SRC_SHARE_DIR}/cxx/prim_cxx_driver_base.F90
    ${SRC_SHARE_DIR}/cxx/utilities/bfb_mod.F90
    ${SRC_SHARE_DIR}/planar_mod.F90
    ${SRC_SHARE_DIR}/geometry_mod.F90
    ${SRC_SHARE_DIR}/planar_mesh_mod.F90
  )

  SET(TEST_SRC_F90
    ${TEST_SRC_DIR}/asp_tests.F90
    ${TEST_SRC_DIR}/baroclinic_inst_mod.F90
    ${TEST_SRC_DIR}/dcmip12_wrapper.F90
    ${TEST_SRC_DIR}/dcmip16_wrapper.F90
    ${TEST_SRC_DIR}/dcmip2012_test1_2_3.F90
    ${TEST_SRC_DIR}/dcmip2012_test1_conv_mod.F90
    ${TEST_SRC_DIR}/dcmip2012_test4.F90
    ${TEST_SRC_DIR}/dcmip2012_test5.F90
    ${TEST_SRC_DIR}/dcmip2016-baroclinic.F90
    ${TEST_SRC_DIR}/dcmip2016-kessler.F90
    ${TEST_SRC_DIR}/dcmip2016-physics-z.F90
    ${TEST_SRC_DIR}/dcmip2016-simple-physics-v5.F90
    ${TEST_SRC_DIR}/dcmip2016-simple-physics-v6.F90
    ${TEST_SRC_DIR}/dcmip2016-supercell.F90
    ${TEST_SRC_DIR}/dcmip2016-terminator.F90
    ${TEST_SRC_DIR}/dcmip2016-tropical-cyclone.F90
    ${TEST_SRC_DIR}/held_suarez_mod.F90
    ${TEST_SRC_DIR}/dry_planar.F90
    ${TEST_SRC_DIR}/moist_planar.F90
    ${TEST_SRC_DIR}/mtests.F90
  )

  SET(PREQX_DEPS_F90
    ${TARGET_F90}
    ${SRC_SHARE_F90}
    ${TEST_SRC_F90}
    ${SRC_DIR}/checksum_mod.F90
    ${SRC_DIR}/common_io_mod.F90
    ${SRC_DIR}/common_movie_mod.F90
    ${SRC_DIR}/interpolate_driver_mod.F90
    ${SRC_DIR}/interp_movie_mod.F90
    ${SRC_DIR}/netcdf_io_mod.F90
    ${SRC_DIR}/pio_io_mod.F90
    ${SRC_DIR}/prim_movie_mod.F90
    ${SRC_DIR}/prim_restart_mod.F90
    ${SRC_DIR}/repro_sum_mod.F90
    ${SRC_DIR}/restart_io_mod.F90
    ${SRC_DIR}/surfaces_mod.F90
    ${SRC_DIR}/test_mod.F90
  )

  SET(PREQX_DEPS_CXX
    ${TARGET_DIR}/cxx/cxx_f90_interface_preqx.cpp
    ${TARGET_DIR}/cxx/CamForcing.cpp
    ${TARGET_DIR}/cxx/Diagnostics.cpp
    ${TARGET_DIR}/cxx/ElementsForcing.cpp
    ${TARGET_DIR}/cxx/ElementsState.cpp
    ${TARGET_DIR}/cxx/HyperviscosityFunctorImpl.cpp
    ${TARGET_DIR}/cxx/prim_advance_exp.cpp
    ${SRC_SHARE_DIR}/cxx/CaarFunctor.cpp
    ${SRC_SHARE_DIR}/cxx/Context.cpp
    ${SRC_SHARE_DIR}/cxx/Elements.cpp
    ${SRC_SHARE_DIR}/cxx/ElementsDerivedState.cpp
    ${SRC_SHARE_DIR}/cxx/ElementsGeometry.cpp
    ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
    ${SRC_SHARE_DIR}/cxx/EulerStepFunctor.cpp
    ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
    ${SRC_SHARE_DIR}/cxx/FunctorsBuffersManager.cpp
    ${SRC_SHARE_DIR}/cxx/GllFvRemap.cpp
    ${SRC_SHARE_DIR}/cxx/GllFvRemapImpl.cpp
    ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
    ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
    ${SRC_SHARE_DIR}/cxx/HyperviscosityFunctor.cpp
    ${SRC_SHARE_DIR}/cxx/ReferenceElement.cpp
    ${SRC_SHARE_DIR}/cxx/Tracers.cpp
    ${SRC_SHARE_DIR}/cxx/VerticalRemapManager.cpp
    ${SRC_SHARE_DIR}/cxx/mpi/BoundaryExchange.cpp
    ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
    ${SRC_SHARE_DIR}/cxx/mpi/Connectivity.cpp
    ${SRC_SHARE_DIR}/cxx/mpi/MpiBuffersManager.cpp
    ${SRC_SHARE_DIR}/cxx/mpi/mpi_cxx_f90_interface.cpp
    ${SRC_SHARE_DIR}/cxx/prim_advec_tracers_remap.cpp
    ${SRC_SHARE_DIR}/cxx/prim_driver.cpp
    ${SRC_SHARE_DIR}/cxx/prim_step.cpp
    ${SRC_SHARE_DIR}/cxx/vertical_remap.cpp
    ${SRC_SHARE_DIR}/cxx/utilities/BfbUtils.cpp
    ${SRC_SHARE_DIR}/cxx/utilities/InternalDiagnostics.cpp
    ${SRC_SHARE_DIR}/cxx/utilities/Hash.cpp
  )

  IF (HOMME_USE_TRILINOS)
    SET (PREQX_SRCS_ZOLTAN   ${TRILINOS_ZOLTAN_DIR}/zoltan_interface.c  ${TRILINOS_ZOLTAN_DIR}/zoltan_cppinterface.cpp)
    SET_SOURCE_FILES_PROPERTIES( ${TRILINOS_ZOLTAN_DIR}/zoltan_cppinterface.cpp PROPERTIES LANGUAGE CXX )
  ENDIF()

  # If the user specified a file for custom compiler options use those
  IF (DEFINED PREQX_CUSTOM_FLAGS_FILE)
    setCustomCompilerFlags(PREQX_CUSTOM_FLAGS_FILE PREQX_SRCS_F90)
  ENDIF ()

  # We pulled prim_main off the deps cause it would create a main function
  # clash when building the unit test (a main in prim_main.F90 and one in tester.cpp)
  SET(PREQX_DEPS
    ${TARGET_F90}
    ${PREQX_DEPS_F90}
    ${PREQX_DEPS_CXX}
  )

  SET(EXEC_SOURCES ${PREQX_DEPS} ${SRC_DIR}/prim_main.F90)

  # Set up defaults if they are not given
  IF (NOT PREQX_NP)
    SET (PREQX_NP 4)
  ENDIF ()
  IF (NOT PREQX_PLEV)
    SET (PREQX_PLEV 20)
  ENDIF ()
  IF (NOT PREQX_USE_PIO)
    SET (PREQX_USE_PIO FALSE)
  ENDIF ()
  IF (NOT PREQX_USE_ENERGY)
    SET (PREQX_USE_ENERGY FALSE)
  ENDIF ()
  #this should be removed later
  IF (NOT PREQX_NC)
    SET (PREQX_NC 4)
  ENDIF ()
  IF (NOT QSIZE_D)
    SET (QSIZE_D 4)
  ENDIF ()
  # for some of matt's old scripts which use preqx_qsize_d
  IF (PREQX_QSIZE_D)      
    SET (QSIZE_D ${PREQX_QSIZE_D})
  ENDIF ()

  SET(USE_OPENACC FALSE)

ENDMACRO(PREQX_KOKKOS_SETUP)
