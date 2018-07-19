function [done] = func_delete_local_copy_of_spc_kazr_ge_2018_0427(...
    raw_netCDF_spc_local_directory, suffix)

% updated = '24-April-2018'

% This routine copys kazr ge spectra from oliktok pt archive to local disk.
% **********************************************************************

%% Delete raw spectra that were copied to local directory

%netCDF_input_directory  = [netCDF_drive,'\Oliktok_Point\ARM_Archive\KAZR\2017_01\olikazrspeccmaskgecopolM1\'];
%raw_netCDF_spc_source_directory  = 'J:\oli\kazr\olikazrspeccmaskgecopolM1\';
%raw_netCDF_spc_local_directory   = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\temp_job1_nc_files\';

wildcard_name           = [raw_netCDF_spc_local_directory,suffix];
files_to_process        = dir(wildcard_name);

% the filenames to process are accessed using:
% files_to_process(i).name

num_files_to_process    = length(files_to_process);

if(num_files_to_process > 0)
   
   % delete the files from the local directory
   disp('...deleting files from local directory...')
   
   for file_loop = 1:num_files_to_process
      %for file_loop = 33:33
      
      % get the next filename to process
      filename_no_dir      = files_to_process(file_loop).name;
      %source_spc_filename  = [raw_netCDF_spc_source_directory,filename_no_dir];
      local_spc_filename   = [raw_netCDF_spc_local_directory,filename_no_dir];
      
      %disp(['deleting ',filename_no_dir,' to local directory...']);
      
      delete(local_spc_filename);
      
   end % end for file_loop
   
end % end if(num_files_to_process > 0)

%% Set output flag

done = 1;

