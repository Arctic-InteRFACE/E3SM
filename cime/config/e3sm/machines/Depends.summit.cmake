list(APPEND NOOPT_FILES
  cam/src/dynamics/eul/dyn_comp.F90
  cam/src/dynamics/fv/dyn_comp.F90
  cam/src/dynamics/se/dyn_comp.F90
  cam/src/dynamics/sld/dyn_comp.F90
  cam/src/physics/cam/microp_aero.F90)

set(FILES_NEED_OPENACC_FLAGS
  cam/src/physics/rrtmgp/external/rte/mo_fluxes.F90
  cam/src/physics/rrtmgp/external/rte/mo_optical_props.F90
  cam/src/physics/rrtmgp/external/rte/mo_rte_kind.F90
  cam/src/physics/rrtmgp/external/rte/mo_rte_lw.F90
  cam/src/physics/rrtmgp/external/rte/mo_rte_sw.F90
  cam/src/physics/rrtmgp/external/rte/mo_rte_util_array.F90
  cam/src/physics/rrtmgp/external/rte/mo_source_functions.F90
  cam/src/physics/rrtmgp/external/rte/kernels/mo_fluxes_broadband_kernels.F90
  cam/src/physics/rrtmgp/external/rte/kernels/mo_optical_props_kernels.F90
  cam/src/physics/rrtmgp/external/rte/kernels/mo_rte_solver_kernels.F90
  cam/src/physics/rrtmgp/external/rte/kernels-openacc/mo_optical_props_kernels.F90
  cam/src/physics/rrtmgp/external/rte/kernels-openacc/mo_rte_solver_kernels.F90
  cam/src/physics/rrtmgp/external/rrtmgp/mo_gas_concentrations.F90
  cam/src/physics/rrtmgp/external/rrtmgp/mo_gas_optics.F90
  cam/src/physics/rrtmgp/external/rrtmgp/mo_gas_optics_rrtmgp.F90
  cam/src/physics/rrtmgp/external/rrtmgp/mo_rrtmgp_constants.F90
  cam/src/physics/rrtmgp/external/rrtmgp/mo_rrtmgp_util_reorder.F90
  cam/src/physics/rrtmgp/external/rrtmgp/mo_rrtmgp_util_string.F90
  cam/src/physics/rrtmgp/external/rrtmgp/kernels/mo_gas_optics_kernels.F90
  cam/src/physics/rrtmgp/external/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
  cam/src/physics/rrtmgp/external/rrtmgp/kernels-openacc/mo_gas_optics_kernels.F90
  cam/src/physics/rrtmgp/external/extensions/mo_compute_bc.F90
  cam/src/physics/rrtmgp/external/extensions/mo_fluxes_byband.F90
  cam/src/physics/rrtmgp/external/extensions/mo_fluxes_byband_kernels.F90
  cam/src/physics/rrtmgp/external/extensions/mo_fluxes_bygpoint.F90
  cam/src/physics/rrtmgp/external/extensions/mo_heating_rates.F90
  cam/src/physics/rrtmgp/external/extensions/mo_rrtmgp_clr_all_sky.F90
  cam/src/physics/rrtmgp/radiation.F90
  cam/src/physics/rrtmgp/radconstants.F90
  cam/src/physics/rrtmgp/cam_optics.F90
  cam/src/physics/rrtmgp/mcica_subcol_gen.F90
  cam/src/physics/rrtmgp/cloud_rad_props.F90
  cam/src/physics/rrtmgp/ebert_curry.F90
  cam/src/physics/rrtmgp/slingo.F90)

foreach(ITEM IN LISTS FILES_NEED_OPENACC_FLAGS)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS "-Minline -ta=nvidia,cc70,fastmath,loadcache:L1,unroll,fma,managed,ptxinfo -Mcuda -Minfo=accel")
endforeach()