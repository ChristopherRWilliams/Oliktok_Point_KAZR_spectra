function [done] = func_copy_spc_from_archive_to_local_disk_kazr_ge_2018_0427(...
   start_year, month_loop, dayofmonth_loop, ...
   raw_netCDF_spc_source_directory, raw_netCDF_spc_local_directory,...
   filename_root, suffix)

% updated = '24-April-2018'

% This routine copys kazr ge spectra from oliktok pt archive to local disk.
% **********************************************************************


%% Define the days to process

year_str       = num2str(start_year);
month_num_str  = num2str_2digits(month_loop);
dayofmonth_str = num2str_2digits(dayofmonth_loop);

disp(['...starting to process day: ',year_str,'-',month_num_str,'-',dayofmonth_str,'...']);

%% Copy the raw spectra to a local directory

%netCDF_input_directory  = [netCDF_drive,'\Oliktok_Point\ARM_Archive\KAZR\2017_01\olikazrspeccmaskgecopolM1\'];
%raw_netCDF_spc_source_directory  = 'J:\oli\kazr\olikazrspeccmaskgecopolM1\';
%raw_netCDF_spc_local_directory   = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\temp_job2_nc_files\';

%filename_root           = 'olikazrspeccmaskgecopolM1.a0.';

wildcard_name           = [raw_netCDF_spc_source_directory,filename_root,year_str,month_num_str,dayofmonth_str,suffix];
files_to_process        = dir(wildcard_name);

% the filenames to process are accessed using:
% files_to_process(i).name

num_files_to_process    = length(files_to_process);

if(num_files_to_process > 0)
   
   % copy the files from the source to a local directory
   
   for file_loop = 1:num_files_to_process
      %for file_loop = 33:33
      
      % get the next filename to process
      filename_no_dir      = files_to_process(file_loop).name;
      source_spc_filename  = [raw_netCDF_spc_source_directory,filename_no_dir];
      local_spc_filename   = [raw_netCDF_spc_local_directory,filename_no_dir];
      
      disp(['copying ',filename_no_dir,' to local directory...']);
      
      copyfile(source_spc_filename,local_spc_filename);
      
   end % end for file_loop
   
end % end if(num_files_to_process > 0)

%% Copy the raw spectra from the previous day to a local directory

%netCDF_input_directory  = [netCDF_drive,'\Oliktok_Point\ARM_Archive\KAZR\2017_01\olikazrspeccmaskgecopolM1\'];
%raw_netCDF_spc_source_directory  = 'J:\oli\kazr\olikazrspeccmaskgecopolM1\';
%raw_netCDF_spc_local_directory   = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\temp_job2_nc_files\';

%filename_root           = 'olikazrspeccmaskgecopolM1.a0.';

prior_day            = dayofmonth_loop - 1;
prior_dayofmonth_str = num2str_2digits(prior_day);
prior_month          = month_loop;
prior_month_num_str  = num2str_2digits(prior_month);

if(prior_day < 1)
   prior_month = month_loop - 1;
   prior_day   = 31;
   prior_month_num_str  = num2str_2digits(prior_month);
   prior_dayofmonth_str = num2str_2digits(prior_day);
   
   wildcard_name           = [raw_netCDF_spc_source_directory,filename_root,year_str,prior_month_num_str,prior_dayofmonth_str,'.23',suffix];
   files_to_process        = dir(wildcard_name);
   num_files_to_process    = length(files_to_process);
   
   if(num_files_to_process == 0)
      prior_day   = 30;
      prior_dayofmonth_str = num2str_2digits(prior_day);
   end % end if(num_files_to_process == 0)
   
end % end if(prior_day < 1)

% get the prior day's filenames
wildcard_name           = [raw_netCDF_spc_source_directory,filename_root,year_str,prior_month_num_str,prior_dayofmonth_str,'.23',suffix];
files_to_process        = dir(wildcard_name);
num_files_to_process    = length(files_to_process);

if(num_files_to_process > 0)
   
   % copy the files from the source to a local directory
   
   for file_loop = 1:num_files_to_process
      %for file_loop = 33:33
      
      % get the next filename to process
      filename_no_dir      = files_to_process(file_loop).name;
      source_spc_filename  = [raw_netCDF_spc_source_directory,filename_no_dir];
      local_spc_filename   = [raw_netCDF_spc_local_directory,filename_no_dir];
      
      disp(['copying ',filename_no_dir,' to local directory...']);
      
      copyfile(source_spc_filename,local_spc_filename);
      
   end % end for file_loop
   
end % end if(num_files_to_process > 0)

%% Set output flag

done = 1;

