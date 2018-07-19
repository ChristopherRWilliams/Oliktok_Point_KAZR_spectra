% Program = 'main_sgp_15s_dwell_spc_kazr_ge_copol_2017_1021'
% updated = '21-October-2017'

% This routine estimates CNR, residual SNR, and clutter velocity stats.

% **********************************************************************

%% Define the directories letters

%netCDF_drive            = 'K:'; % home computer
%mat_drive               = 'C:'; % home computer

%netCDF_drive            = 'I:'; % office computer
%mat_drive               = 'E:'; % office computer

%make_video_flag         = 0;

site_prefix             = 'sgp';
version_suffex          = '_v15';

%% Define the days to process

start_year  = 2011;
start_month =    5;
end_month   =    5;
start_day   =   20;
end_day     =   20;
start_hour  =   10;
end_hour    =   10;

%% Loop through each day in list


for month_loop = start_month:end_month
   
   for dayofmonth_loop = start_day:1:end_day
      
      year_str       = num2str(start_year);
      month_num_str  = num2str_2digits(month_loop);
      dayofmonth_str = num2str_2digits(dayofmonth_loop);
      
      day_flag = 0;
      
      for hour_loop  = start_hour:end_hour
         
         %% Define the day to be processed
         
         %year_str       = num2str(start_year);
         %month_num_str  = num2str_2digits(month_loop);
         %dayofmonth_str = num2str_2digits(dayofmonth_loop);
         hour_str       = num2str_2digits(hour_loop);
         
         date_str       = [year_str,'-',month_num_str,'-',dayofmonth_str,', ',hour_str,' UTC'];
         %disp(['...loading spectra for ',date_str,'...']);
         
         %% Has this hour already been processed?
         
         % define the input hourly files
         output_file_directory      = '..\mat_15sec_ave_moments\';
         %output_file_daily          = [output_file_directory,'oli_kazr_ge_15sec_mom_daily_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         
         % define the output daily file for only the single peak moments
         output_file_daily_single   = [output_file_directory,site_prefix,'_kazr_ge_15sec_mom_',year_str,month_num_str,dayofmonth_str,'_',hour_str,version_suffex,'.mat'];
         
         good_output_file_daily_single = exist(output_file_daily_single,'file');
                  
         if(good_output_file_daily_single ~= 2)
      
            day_flag = 1;
            
            %% load the hourly spc file
            
            % check to see if the hourly spc file has already been made.
            input_directory        = '..\mat_hourly_spc_files\';
            input_filename_root    = [site_prefix,'_kaz_ge_hourly_raw_spc_'];
            input_filename         = [input_directory,input_filename_root,year_str,month_num_str,dayofmonth_str,'_',hour_str,'.mat'];
            
            good_hour_spc_mat_file = exist(input_filename,'file');
            
            if(good_hour_spc_mat_file == 2)
               
               disp(['...loading matlab spectra, ',date_str,'...']);
               
               % load the mat file
               load(input_filename);
               % the variables are:
               %   kaz_ge_Nspc                            1x1                           8  double
               %   kaz_ge_Vd                            512x1                        4096  double
               %   kaz_ge_Z_near_field_correction         1x607                      4856  double
               %   kaz_ge_orig_spc_lin                 1771x607x512            4403187712  double
               %   kaz_ge_orig_time                    1771x7                       99176  double
               %   kaz_ge_radar_const_dB                  1x1                           8  double
               %   kaz_ge_range                         607x1                        4856  double            
               
               %% Estimate the moments before doing spectrum declutter
               
               % These moments will be used to compare with the spc
               % decluttered moments.
               
               % Check if original moments have already been calculated.
               % If not, then calculate orig moments.
               
               %output_directory  = [mat_drive,'\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\mat_orig_mom\'];
               output_directory  = '..\mat_orig_mom\';
               output_filename   = [output_directory,site_prefix,'_kazr_ge_mom_',year_str,month_num_str,dayofmonth_str,'_',hour_str,version_suffex,'.mat'];
               
               good_orig_mom_filename = exist(output_filename,'file');
               
               if(good_orig_mom_filename ~= 2)
                  
                  disp(' ')
                  disp(['Calculating moments of original spectra...',date_str,'...'])
                  
                  % Get the original spectra in linear units
                  %spc_orig_time_lin = kaz_ge_orig_spc_lin;
                  
                  spc_input         = kaz_ge_orig_spc_lin;
                  [mm,nn,~]         = size(spc_input);
                  cnt_spc_input     = ones(mm,nn);
                  kaz_ge_time_stamp = kaz_ge_orig_time;
                  %kaz_ge_range      = kaz_ge_range;
                  
                  %max_cnt_spc    = max(max(cnt_spc_30sec));
                  cnt_threshold     = 1;
                  
                  % input variables
                  Vd_orig        = kaz_ge_Vd;
                  Nspc           = kaz_ge_Nspc;
                  ht_range       = kaz_ge_range;
                  ht_range_20dB  = 20 .* log10(ht_range);
                  
                  % Define how many valid points are needed before
                  % estimating moments
                  valid_pts_thres = 5;
                  
                  %radar_const_dB = kaz_ge_radar_const_dB;
                  %Z_near_field_correction = kaz_ge_Z_near_field_correction;
                  
                  
                  % Return selected moments which will will be used later:
                  % spc_dcl_snr & spc_dcl_Vpeak - used in clutter moment filter
                  % spc_dcl_Vmean - used for 15 sec dwell processing
                  
                  [spc_orig_snr, spc_orig_Vmean, spc_orig_Vpeak] = ...
                     func_calc_multi_peak_ge_moments(spc_input, cnt_spc_input, ...
                     Vd_orig, Nspc, valid_pts_thres, ht_range_20dB, cnt_threshold, kaz_ge_radar_const_dB, ...
                     kaz_ge_Z_near_field_correction, output_filename, kaz_ge_time_stamp, kaz_ge_range);
                  
               end % end if(good_orig_mom_filename ~= 2)
                              
               %% Process 15sec dwell spectra
               % ============================
               
               % define the number of valid points in spectra before
               % estimating moments
               %valid_pts_thres   = 5;
               
               % save just the 15sec moments
               %[~] = func_ave_spc_to_15sec_and_calc_mom(spc_orig_time_lin, ...
               %   kaz_ge_orig_time, spc_dcl_Vmean, hour_loop,...
               %   mat_drive, year_str, month_num_str, dayofmonth_str, hour_str, Vd_orig, ...
               %   Nspc, ht_range_20dB, kaz_ge_radar_const_dB, ...
               %   kaz_ge_Z_near_field_correction, kaz_ge_range, version_suffex);
               %disp('paused...')
               %pause

               % save both the 15sec moments and the 15sec spectra
               % kaz_ge_orig_spc_lin
               [~] = func_ge_ave_spc_to_15sec_and_calc_mom_and_save_spc(kaz_ge_orig_spc_lin, ...
                  kaz_ge_orig_time, spc_orig_Vmean, hour_loop,...
                  year_str, month_num_str, dayofmonth_str, hour_str, Vd_orig, ...
                  Nspc, valid_pts_thres, ht_range_20dB, kaz_ge_radar_const_dB, ...
                  kaz_ge_Z_near_field_correction, kaz_ge_range, site_prefix, version_suffex);
               
            end % end if(good_hour_spc_mat_file == 2)
            
         end % end if(good_output_file_daily_single == 2)
         
      end % end for hour_loop
      
      % For this day, generate daily mat files
      
      %% Construct daily files  -  Straight Average Dwells
      %  =================================================
      
      if(day_flag)
         
         %% Combine hourly files and save as a daily file: 15sec dwell
         
         disp(' ')
         disp('...combining 15sec hourly files into one daily file...');
         
         % Process the straight average
         % =======================
         
         % define the input hourly files
         mat_input_directory        = '..\mat_15sec_ave_moments\';
         filename_root              = [site_prefix,'_kazr_ge_15sec_mom_',year_str,month_num_str,dayofmonth_str];
         
         % define the output daily file for all variables
         output_file_directory      = '..\mat_daily_15sec_ave_moments\';
         output_file_daily          = [output_file_directory,site_prefix,'_kazr_ge_15sec_mom_daily_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         
         % define the output daily file for only the single peak moments
         output_file_daily_single   = [output_file_directory,site_prefix,'_kazr_ge_15sec_mom_daily_single_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         
         % call the function that combines these files
         [~] = func_make_daily_mat_files(mat_input_directory, filename_root, ...
            output_file_daily, output_file_daily_single);
         
         % Process the shift & average
         % ======================
         
         % define the input hourly files
         mat_input_directory        = '..\mat_15sec_ave_moments\';
         filename_root              = [site_prefix,'_kazr_ge_15sec_shift_mom_',year_str,month_num_str,dayofmonth_str];
         
         % define the output daily file for all variables
         output_file_directory      = '..\mat_daily_15sec_ave_moments\';
         output_file_daily          = [output_file_directory,site_prefix,'_kazr_ge_15sec_shift_mom_daily_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         
         % define the output daily file for only the single peak moments
         output_file_daily_single   = [output_file_directory,site_prefix,'_kazr_ge_15sec_shift_mom_daily_single_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         
         % call the function that combines these files
         [~] = func_make_daily_mat_files(mat_input_directory, filename_root, ...
            output_file_daily, output_file_daily_single);
         
         %% Combine hourly files and save as a daily file: original dwell
         
         %    disp(' ')
         %    disp('...combining original dwell hourly files into one daily file...');
         %
         %    % Process the straight average
         %    % =======================
         %
         %    % define the input hourly files
         %    mat_input_directory        = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\mat_orig_moments\';
         %    filename_root              = ['oli_kazr_ge_orig_mom_',year_str,month_num_str,dayofmonth_str];
         %
         %    % define the output daily file for all variables
         %    output_file_directory      = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\mat_daily_orig_moments\';
         %    output_file_daily          = [output_file_directory,'oli_kazr_ge_orig_mom_daily_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         %
         %    % define the output daily file for only the single peak moments
         %    output_file_daily_single   = [output_file_directory,'oli_kazr_ge_orig_mom_daily_single_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         %
         %    % call the function that combines these files
         %    [~] = func_make_daily_mat_files(mat_input_directory, filename_root, ...
         %       output_file_daily, output_file_daily_single);
         
         %% Plot the 15sec moments
         
         disp(' ')
         disp('...generating daily plots of Z, Vmean, Vsig, and Vskew...')
         
         max_ht      = 8;
         hourly_flag = 0;
         
         % define the output daily file for all variables
         input_file_directory       = '..\mat_daily_15sec_ave_moments\';
         input_filename             = [input_file_directory,site_prefix,'_kazr_ge_15sec_shift_mom_daily_single_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         
         output_directory           = '..\images_15sec_ave_moments\';
         
         % plot the Z, V, and std(V)
         output_filename_root       = [site_prefix,version_suffex,'_kazr_ge_15sec_ZVW_shift_'];
         [~] = func_plot_15sec_kazr_ZVW_mom_ge_copol(input_filename, output_directory, output_filename_root, site_prefix, max_ht, hourly_flag);
         
         % plot the Z, V, and skewness
         output_filename_root       = [site_prefix,version_suffex,'_kazr_ge_15sec_ZVSkew_shift_'];
         [~] = func_plot_15sec_kazr_ZVSkew_mom_ge_copol(input_filename, output_directory, output_filename_root, site_prefix, max_ht, hourly_flag);
         
         % define the output daily file for all variables
         input_file_directory       = '..\mat_daily_15sec_ave_moments\';
         input_filename             = [input_file_directory,site_prefix,'_kazr_ge_15sec_mom_daily_single_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
         
         % plot the Z, V, and skewness
         output_filename_root       = [site_prefix,version_suffex,'_kazr_ge_15sec_ZVSkew_'];
         [~] = func_plot_15sec_kazr_ZVSkew_mom_ge_copol(input_filename, output_directory, output_filename_root, site_prefix, max_ht, hourly_flag);
         
      end % end if(day_flag)
      
   end % end for day_loop
   
   
end % end for month_loop

