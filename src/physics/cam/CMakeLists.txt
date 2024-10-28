# Append sources from this directory to the cam_sources list.
list(APPEND cam_sources micro_mg_data.F90 micro_mg_utils.F90
  vdiff_lu_solver.F90 corrector.F90 conv_state_swap.F90)

sourcelist_to_parent(cam_sources)
