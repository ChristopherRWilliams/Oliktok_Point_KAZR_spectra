% Program = 'main_make_oli_hourly_mat_spc_kazr_ge_copol_2018_0606'
% updated = '06-June-2018'

% This routine process KAZR spectra and estimates time-averaged moments.

% **********************************************************************

%% Define the directories letters

%netCDF_drive            = 'K:'; % home computer
%mat_drive               = 'C:'; % home computer

%netCDF_drive            = 'I:'; % office computer
%mat_drive               = 'E:'; % office computer

%make_video_flag         = 0;

%version_suffex          = '_20170921';

%% Define directories:

% Define archive directory with raw spectra (from ARM Archive)
%raw_netCDF_spc_source_directory  = '..\raw_netCDF\';
raw_netCDF_spc_source_directory  = 'D:\oli\kazr\olikazrspeccmaskgecopolM1\';
%raw_netCDF_spc_source_directory  = 'E:\Projects\ARM_Summer_Camp\KAZR\raw_netCDF\';

% define local directory to temporarily store a day's worth of spectra
% files in the directory are *deleted* after processing each day
raw_netCDF_spc_local_directory   = '..\temp\';
%raw_netCDF_spc_local_directory   = 'E:\Projects\ARM_Summer_Camp\KAZR\temp\';

% define where the hourly (*.mat) spectra will be stored 
% (input and output directories are the same)
input_directory        = '..\mat_hourly_spc_files\';
output_directory       = '..\mat_hourly_spc_files\';

% define the filename_root for the hourly files         
% (input and output filename_root names are the same)
input_filename_root    = 'oli_kaz_ge_hourly_raw_spc_';
output_filename_root   = 'oli_kaz_ge_hourly_raw_spc_';

%% Define the days to process

start_year  = 2017;
start_month =   11;
end_month   =   11;
start_day   =    1;
end_day     =    1;
start_hour  =    0;
end_hour    =   23;

%% Loop through each day in list

   
for month_loop = start_month:end_month

   for dayofmonth_loop = start_day:end_day

      year_str       = num2str(start_year);
      month_num_str  = num2str_2digits(month_loop);
      dayofmonth_str = num2str_2digits(dayofmonth_loop);
      
      %% copy raw spectra from archive to a local directory
      
      filename_root           = 'olikazrspeccmaskgecopolM1.a0.';
      suffix                  = '*.nc';

      [~] = func_copy_spc_from_archive_to_local_disk_kazr_ge_2018_0427(...
         start_year, month_loop, dayofmonth_loop, ...
         raw_netCDF_spc_source_directory, raw_netCDF_spc_local_directory,...
         filename_root, suffix);
      
      for hour_loop  = start_hour:end_hour
         
         %% Define the day to be processed
         
         %year_str       = num2str(start_year);
         %month_num_str  = num2str_2digits(month_loop);
         %dayofmonth_str = num2str_2digits(dayofmonth_loop);
         hour_str       = num2str_2digits(hour_loop);
         
         date_str       = [year_str,'-',month_num_str,'-',dayofmonth_str,', ',hour_str,' UTC'];
         disp(['...loading spectra for ',date_str,'...']);
         
         %% load the hourly spc file
         
         % check to see if the hourly spc file has already been made.
         %input_directory        = 'C:\oli\kazr\mat_hourly_spc_kazr_ge_raw\';
         %input_filename_root    = 'oli_kaz_ge_hourly_raw_spc_';
         input_filename         = [input_directory,input_filename_root,year_str,month_num_str,dayofmonth_str,'_',hour_str,'.mat'];
         
         good_hour_spc_mat_file = exist(input_filename,'file');
         
         if(good_hour_spc_mat_file == 2)
            
            %             % load the mat file
            %             load(input_filename);
            %             % the variables are:
            %             %   kaz_ge_Nspc                            1x1                           8  double
            %             %   kaz_ge_Vd                            512x1                        4096  double
            %             %   kaz_ge_Z_near_field_correction         1x607                      4856  double
            %             %   kaz_ge_orig_spc_lin                 1771x607x512            4403187712  double
            %             %   kaz_ge_orig_time                    1771x7                       99176  double
            %             %   kaz_ge_radar_const_dB                  1x1                           8  double
            %             %   kaz_ge_range                         607x1                        4856  double
            
         else
                        
            %% Construct hourly spectra files
            
            [kaz_ge_orig_spc_lin, kaz_ge_orig_time,kaz_ge_Vd,...
               kaz_ge_range, kaz_ge_radar_const_dB, kaz_ge_Z_near_field_correction,...
               kaz_ge_Nspc] = func_get_hourly_spc_into_mat_kazr_ge_copol_2018_0427(...
               start_year, month_loop, dayofmonth_loop, hour_loop, ...
               raw_netCDF_spc_local_directory, suffix);
            
            % The variables are:
            %   kaz_ge_Nspc                            1x1                           8  double
            %   kaz_ge_Vd                            512x1                        4096  double
            %   kaz_ge_Z_near_field_correction         1x607                      4856  double
            %   kaz_ge_orig_spc_lin                 1771x607x512            4403187712  double
            %   kaz_ge_orig_time                    1771x7                       99176  double
            %   kaz_ge_radar_const_dB                  1x1                           8  double
            %   kaz_ge_range                         607x1                        4856  double
            
            %% verify that there are observations for this hour
            
            [mm, nn, pp]     = size(kaz_ge_orig_spc_lin);
            
            if( (mm > 1) && (nn > 1) && (pp > 1) )
               process_hour_flag = 1;
            else
               process_hour_flag = 0;
            end % end if(pp > 1)
            
            if(process_hour_flag)
               
               %% Estimate the moments before doing spectrum declutter
               
               % These moments will be used to compare with the spc
               % decluttered moments.
               
               disp(' ')
               disp('Saving the original spectra in matlab format...')
               
               %output_directory        = 'I:\oli\kazr\mat_hourly_spc_kazr_ge_raw\';
               %output_filename_root    = 'oli_kaz_ge_hourly_raw_spc_';
               
               output_filename         = [output_directory,output_filename_root,year_str,month_num_str,dayofmonth_str,'_',hour_str,'.mat'];
               
               save(output_filename,'kaz_ge_Nspc','kaz_ge_Vd',...
                  'kaz_ge_Z_near_field_correction','kaz_ge_orig_spc_lin',...
                  'kaz_ge_orig_time','kaz_ge_radar_const_dB',...
                  'kaz_ge_range','-v7.3')
               
            end % end if(process_hour_flag)
            
            
         end % end if(good_hour_spc_mat_file == 2)
         
      end % end for hour_loop      

      %% deletey raw spectra from local directory
      
      [~] = func_delete_local_copy_of_spc_kazr_ge_2018_0427(...
         raw_netCDF_spc_local_directory, suffix);
      
   end % end for day_loop
   
   
end % end for month_loop

