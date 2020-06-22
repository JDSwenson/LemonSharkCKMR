# LemonSharkCKMR
Use the following naming convention when naming scripts:
[prior]_[script_type(data, null, simulation description)]_loop_[CURRENT or date]

e.g.: HS_lemon_data_loop_CURRENT; PO_prop_mat_vs_samp_size_loop_12.18.19; HS_time_series_null_loop_CURRENT

Workflow:
1) Write and test scripts in respective subfolders (e.g. time_series, constant_abundance, etc)
2) Use Viz_script in misc_scripts to plot, and Model_analysis script to investigate results more closely
3) Copy and paste plots into Results_and_figures doc on Box (Manuscripts/CKMR-Validation/Results) and write caption
4) Copy simulation script and viz script into CKMR_script_repository on Box (Manuscripts/CKMR-Validation/Scripts)
