# Oliktok_Point_KAZR_spectra
Declutter spectra, shift-and-average spectra, and calculate multiple peak moments
35-GHz Vertically Pointing Radar Processing at Oliktok Point, Alaska 

Author: 	Christopher R. Williams
Updated: 	01-August-2018

This repository processes vertically pointing radar Doppler velocity spectra and estimates high-order moments of multiple peaks. The processing is described in the manuscript “Clutter Mitigation, Multiple Peaks, and High-Order Spectral Moments in 35-GHz Vertically Pointing Radar Velocity Spectra”, submitted to Atmospheric Measurement Techniques (amt-2018-66).As an overview, this processing performs:
  1.	Remove ground clutter in velocity spectra
  2.	Estimate multiple peaks and high order moments
  3.	Average velocity spectra into 15-s intervals
  4.	Estimate multiple peaks and high order moments of 15-s averaged spectra 
  
# Processing Steps
The processing consists of three different steps.
- 1. Pre-processing routine
  - a. Read partial hours of raw spectra in netCDF format and save as hourly Matlab files
- 2. Main processing routine 
  - a. Read hourly raw spectra files
  - b. Declutter spectra
  - c. Estimate multiple-peak moments on decluttered spectra at original time resolution
  - d. Construct 15-s shift-and-average spectra
  - e. Estimate multiple-peak moments on 15-s averaged spectra
  - f. Concatenate hourly files into daily Matlab files
- 3. Post-process routines
  - a. Read daily matlab processed files and save in netCDF format
  - b. Generate daily and hourly tif images using the netCDF data files

# Directory Structure

The processing routines expect to find intermediate and final files in specific directories. Thus, the following directories need to be defined before running the routines: 

- /images_15sec_ave_moments 	- images generated after processing day of data
- /images_from_nc_files/daily 	- daily images generated from netCDF files
- /images_from_nc_files/hourly	- hourly images generated from netCDF files
- /m_files 			- all matlab m-files (routines) are in this directory
- /mat_15sec_ave_moments  	- hourly files, moments estimated from 15s averaged spectra
- /mat_15sec_ave_spc 		- hourly files, 15s averaged spectra
- /mat_clutter_stats 		- hourly files, clutter statistics
- /mat_daily_15sec_ave_moments 	- daily files, concatenation of hourly moment files
- /mat_dcl_mom 		- hourly files, decluttered moments for each profile
- /mat_hourly_spc_files 	- hourly raw spectra files in Matlab format (input files for main routine)
- /mat_orig_mom 		- hourly files, moments before decluttering spectra
- /nc_daily_15sec_ave_moments 	- daily netCDF files of 15sec averaged moments
- /raw_netCDF 		- raw spectra in netCDF format (downloaded from ARM archive)
- /temp	- temporary files are stored and then deleted after processing each day

# Processing Files

All processing files are written in Matlab and are in the directory:
- /m_files 			- All Matlab m-files (routines) are in this directory

# Pre-Processing Routine

Routine: **main_make_oli_hourly_mat_spc_kazr_ge_copol_2018_0606.m**
This pre-processing routine reads partial hours of raw spectra written in netCDF format and saves the spectra in hourly matlab files. These intermediate hourly matlab files are read by the main processing routine so that it can work with complete hours of data. 

Functions needed for this routine:
- func_convert_netCDF_time_to_vector.m
- func_copy_spc_from_archive_to_local_disk_kazr_ge_2018_0427.m
- func_delete_local_copy_of_spc_kazr_ge_2018_0427.m
- func_get_hourly_spc_into_mat_kazr_ge_copol_2018_0427.m
- func_read_kazr_netCDF_file_no_spectra.m
- func_read_kazr_spc_time_in_netCDF_file.m

Directories needed for this routine:
- /raw_netCDF 	- directory containing multiple days of raw spectra in netCDF format (downloaded from ARM archive)
- /temp	- directory to store daily raw spectra. Spectra deleted from this directory after processing day of spectra.

Data files from ARM Archive. These routines process KAZR spectra from the general, co-pol mode with the following filename structure:
- olikazrspeccmaskgecopolM1.a0.YYYYMMDD.HHMMSS.nc

Place multiple days of raw spectra files in the directory:
- /raw_netCDF 	- directory containing multiple days of raw spectra in netCDF format (downloaded from ARM archive)

# Main Processing Routine

Routine: **main_make_oli_hourly_mat_spc_kazr_ge_copol_2018_0606.m**
This routine processes the Oliktok Point KAZR spectra as described in the Williams et al. (2018) AMT manuscript.

Functions needed for this routine:
- func_calc_multi_peak_ge_moments.m
- func_calc_Vmean_var.m
- func_declutter_spc.m 	
- func_find_mean_HS_noise.m
- func_find_moments.m
- func_find_multi_mom_3spc_Vmean_prior.m
- func_find_noise_adjusted_zdb_and_snr.m
- func_find_single_peak.m
- func_find_single_peak_Vmean_prior_valley.m
- func_ge_ave_spc_to_15sec_and_calc_mom_and_save_spc.m
- func_incoherent_ave_spc_valid_obs.m
- func_make_daily_mat_files.m
- func_plot_15sec_kazr_ZVSkew_mom_ge_copol.m
- func_plot_15sec_kazr_VZW_mom_ge_copol.m
- func_shift_then_ave_spc.m

# Post-Processing Routines

Routine: **main_save_daily_ge_15sec_all_mom_as_netCDF_2018_0718.m**
This post-processing routine reads the daily matlab file and generates a daily netCDF file.

Functions needed for this routine:
- func_remove_isolated_pixels_3x3.m
- func_replace_Nstd_nos_with_NaN.m
- func_save_ge_all_mom_in_netCDF.m

Routine: **main_plot_ge_15sec_mom_nc_all_ZVW_2018_0718.m**
This post-processing routine reads the daily netCDF file and generates daily and hourly images. 

Functions needed for this routine:
- func_fill_time_gaps_with_NaN_profiles.m
- func_read_nc_file.m

# Final Products
The main output products from this processing are:
- /nc_daily_15sec_ave_moments 	- daily netCDF files of 15sec averaged moments
- /images_from_nc_files/daily 		- daily images generated from netCDF files
- /images_from_nc_files/hourly		- hourly images generated from netCDF files

# Data Files
Due to the size of the raw data files, the user must download raw data files from the DOE ARM Archive. The Matlab code is written to process the KAZR spectra from Oliktok Point using the general, co-pol mode. The filenames for this file type is: “olikazrspeccmaskgecopolM1.a0.yyyymmdd.hhmmss.nc”.
Before running this Matlab code, the use must download the raw spectra data for 20-June-2016 and place those files in the directory:
- /raw_netCDF 	- directory containing multiple days of raw spectra in netCDF format (downloaded from ARM archive)
