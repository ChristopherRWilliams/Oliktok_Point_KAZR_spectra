% Program = 'main_oli_dcl_15s_spc_kazr_ge_copol_v15'
% updated = '18-July-2018'

% This routine performs decluttering, spectra shift-then-average, and
% multiple peak moment estimation.

% **********************************************************************

%% Define the directories letters

%make_video_flag         = 0;

site_prefix             = 'oli';
version_suffex          = '_v15';

%% Define the days to process

start_year  = 2016;
start_month =    6;
end_month   =    6;
start_day   =   19;
end_day     =   19;
start_hour  =   12;
end_hour    =   13;

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
         
         % process this hour if it has not been processed, yet
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
               
               %% The following code is the decluttering code
               % ============================================
               
               %% Do decluttering of spectra
               
               disp(' ')
               disp(['...Decluttering spectra...',date_str,'...'])
               
               [m,n,p]                                = size(kaz_ge_orig_spc_lin);
               
               n_clt                                  = 100;
               % Define a decluttered version of the spectra
               %kaz_ge_dcl_spc_dB                   = ones(m,n,p) .* NaN;
               kaz_ge_dcl_spc_lin                     = ones(m,n,p) .* NaN;
               
               kaz_ge_mean_noise                      = ones(m,n_clt) .* NaN;
               kaz_ge_HS_thres_noise                  = ones(m,n_clt) .* NaN;
               kaz_ge_std_noise                       = ones(m,n_clt) .* NaN;
               
               kaz_ge_clt_1bin_peak_dB                = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_drop_dB                = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_cnr_dB                 = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_snr_dB        = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_Vpeak_index   = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_Vpeak         = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_drop_dB       = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_num_valid_pts = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_lt_index      = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_rt_index      = ones(m,n_clt) .* NaN;
               kaz_ge_clt_1bin_residual_sacr_flag     = ones(m,n_clt) .* NaN;
               
               % Define some constants
               % input variables
               Vd             = kaz_ge_Vd;
               Nspc           = kaz_ge_Nspc;
               %ht_range       = kaz_ge_range;
               
               % Define how many valid points are needed before
               % estimating moments
               valid_pts_thres = 5; % this many or more
               
               % Define the drop to the first neighbor threshold
               drop_dB_thres  = 3;
               
               %pause
               
               for r = 1:m
                  %for r = 1:1
                  
                  if(rem(r,250) == 0)
                     disp(['...processing: ',num2str(r),' out of ',num2str(m), ' profiles...']);
                  end % end if(rem(r,250) == 0)
                  
                  % for each profile, do the decluttering
                  prof_spc_orig_lin    = squeeze(kaz_ge_orig_spc_lin(r,:,:));
                  
                  % declutter profile
                  % testing variables
                  %Vd    = Vd_orig;
                  
                  % This code calculates stats for 1bin drop & calculates
                  % decluttered spectra
                  
                  [prof_noise_stats, ...
                     prof_clt_peak_dB, prof_clt_drop_dB, prof_clt_cnr_dB, ...
                     prof_clt_residual_snr_dB, prof_clt_residual_Vpeak,...
                     prof_clt_residual_Vpeak_index, ...
                     prof_clt_residual_drop_dB, ...
                     prof_clt_residual_num_valid_pts,...
                     prof_clt_residual_lt_index, prof_clt_residual_rt_index,...
                     prof_clt_residual_sacr_flag, prof_spc_dcl_lin] = ...
                     func_declutter_spc(prof_spc_orig_lin, Vd, Nspc, n_clt, ...
                     drop_dB_thres, valid_pts_thres);
                  
                  
                  % save the values for this spectrum
                  % =================================
                  % save the noise stats
                  kaz_ge_mean_noise(r,:)      = 10.*log10(prof_noise_stats(:,1));
                  kaz_ge_HS_thres_noise(r,:)  = 10.*log10(prof_noise_stats(:,2));
                  kaz_ge_std_noise(r,:)       = prof_noise_stats(:,3);
                  
                  kaz_ge_clt_1bin_peak_dB(r,:)                 = prof_clt_peak_dB;
                  kaz_ge_clt_1bin_drop_dB(r,:)                 = prof_clt_drop_dB;
                  kaz_ge_clt_1bin_cnr_dB(r,:)                  = prof_clt_cnr_dB;
                  kaz_ge_clt_1bin_residual_snr_dB(r,:)         = prof_clt_residual_snr_dB;
                  kaz_ge_clt_1bin_residual_Vpeak_index(r,:)    = prof_clt_residual_Vpeak_index;
                  kaz_ge_clt_1bin_residual_Vpeak(r,:)          = prof_clt_residual_Vpeak;
                  kaz_ge_clt_1bin_residual_drop_dB(r,:)        = prof_clt_residual_drop_dB;
                  kaz_ge_clt_1bin_residual_num_valid_pts(r,:)  = prof_clt_residual_num_valid_pts;
                  kaz_ge_clt_1bin_residual_lt_index(r,:)       = prof_clt_residual_lt_index;
                  kaz_ge_clt_1bin_residual_rt_index(r,:)       = prof_clt_residual_rt_index;
                  kaz_ge_clt_1bin_residual_sacr_flag(r,:)      = prof_clt_residual_sacr_flag;
                  kaz_ge_dcl_spc_lin(r,:,:)                    = prof_spc_dcl_lin;
                  
               end % end for r loop
               
               %% Save the clutter statistics
               
               %output_directory  = [mat_drive,'\Projects\Oliktok_Point\KAZR\multi_peak_clutter_ave_spc_2017_0921\mat_clutter_stats\'];
               output_directory  = '..\mat_clutter_stats\';
               output_filename   = [output_directory,site_prefix,'_kazr_ge_1bin_stats_',year_str,month_num_str,dayofmonth_str,'_',hour_str,'.mat'];
               
               kaz_ge_clt_1bin_range   = kaz_ge_range(1:n_clt);
               
               save(output_filename,'kaz_ge_orig_time','kaz_ge_Vd',...
                  'kaz_ge_clt_1bin_range','kaz_ge_mean_noise',...
                  'kaz_ge_HS_thres_noise','kaz_ge_std_noise',...
                  'kaz_ge_clt_1bin_peak_dB','kaz_ge_clt_1bin_drop_dB',...
                  'kaz_ge_clt_1bin_cnr_dB','kaz_ge_clt_1bin_residual_snr_dB',...
                  'kaz_ge_clt_1bin_residual_Vpeak_index','kaz_ge_clt_1bin_residual_Vpeak',...
                  'kaz_ge_clt_1bin_residual_drop_dB','kaz_ge_clt_1bin_residual_num_valid_pts',...
                  'kaz_ge_clt_1bin_residual_lt_index','kaz_ge_clt_1bin_residual_rt_index',...
                  'kaz_ge_clt_1bin_residual_sacr_flag')
               
               %% Save the decluttered spectra
               
               %                output_directory  = '..\mat_dcl_spc\';
               %                output_filename   = [output_directory,site_prefix,'_kazr_ge_dcl_spc_',year_str,month_num_str,dayofmonth_str,'_',hour_str,'.mat'];
               %
               %                % Replace noise values with NaN
               %                %input_spc_lin  = kaz_ge_dcl_spc_lin;
               %                %[m,n,~]        = size(kaz_ge_dcl_spc_lin);
               %                %input_Nspc     = ones(m,n) .* kaz_ge_Nspc;
               %                %input_Nstd     = 2;
               %
               %                %[kaz_ge_mean_noise, kaz_ge_std_noise, ...
               %                %   kaz_ge_HS_thres_noise, kaz_ge_dcl_spc_NaN_lin] = ...
               %                %   func_replace_Nstd_nos_with_NaN(input_spc_lin, input_Nspc, input_Nstd);
               %
               %                ht_range_20dB  = 20 .* log10(kaz_ge_range);
               %
               %                save(output_filename,'kaz_ge_orig_time','kaz_ge_mean_noise',...
               %                   'kaz_ge_HS_thres_noise','kaz_ge_std_noise',...
               %                   'kaz_ge_dcl_spc_lin','kaz_ge_Vd','kaz_ge_range',...
               %                   'ht_range_20dB','kaz_ge_Nspc','kaz_ge_radar_const_dB',...
               %                   'kaz_ge_Z_near_field_correction','-v7.3')
                              
               %% Estimate the moments using the raw (original) dwell times
               
               disp(' ')
               disp(['Calculating moments of decluttered spectra...',date_str,'...'])
               
               % Get the original spectra in linear units
               %spc_orig_time_lin = 10.^(kaz_ge_dcl_spc_dB./10);
               spc_orig_time_lin = kaz_ge_dcl_spc_lin;
               
               spc_input         = spc_orig_time_lin;
               [mm,nn,~]         = size(spc_input);
               cnt_spc_input     = ones(mm,nn);
               kaz_ge_time_stamp = kaz_ge_orig_time;
               %kaz_ge_range      = kaz_ge_range;
               
               Vd_orig        = kaz_ge_Vd;
               Nspc           = kaz_ge_Nspc;
               ht_range       = kaz_ge_range;
               ht_range_20dB  = 20 .* log10(ht_range);
               
               % Define how many valid points are needed before
               % estimating moments
               %valid_pts_thres = 5;
               
               %max_cnt_spc    = max(max(cnt_spc_30sec));
               cnt_threshold     = 1;
               
               %output_directory  = [mat_drive,'\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\mat_orig_mom\'];
               output_directory  = '..\mat_dcl_mom\';
               %output_filename   = [output_directory,'oli_kazr_ge_dcl_vpt5_mom_',year_str,month_num_str,dayofmonth_str,'_',hour_str,version_suffex,'.mat'];
               output_filename   = [output_directory,site_prefix,'_kazr_ge_dcl_mom_',year_str,month_num_str,dayofmonth_str,'_',hour_str,version_suffex,'.mat'];
               
               good_dcl_mom_filename = exist(output_filename,'file');
               
               % if this hour has been processed before, get the moments,
               % else, process this hour
               if(good_dcl_mom_filename == 2)
                  
                  disp(' ');
                  disp('...loading pre-calculated dcl moments: snr, Vmean, & Vpeak...');
                  
                  load(output_filename,'kaz_ge_single_snr','kaz_ge_single_Vmean','kaz_ge_single_Vpeak');
                  
                  % rename the variables
                  spc_dcl_snr    = kaz_ge_single_snr;
                  spc_dcl_Vmean  = kaz_ge_single_Vmean;
                  spc_dcl_Vpeak  = kaz_ge_single_Vpeak;
                  
               else
                  
                  disp(' ')
                  disp('...Calculating dcl moments: snr, Vmean, & Vpeak...');
                  
                  % Return selected moments which will will be used later:
                  % spc_dcl_snr & spc_dcl_Vpeak - used in clutter moment filter
                  % spc_dcl_Vmean - used for 15 sec dwell processing
                  
                  [spc_dcl_snr, spc_dcl_Vmean, spc_dcl_Vpeak] = ...
                     func_calc_multi_peak_ge_moments(spc_input, cnt_spc_input, ...
                     Vd_orig, Nspc, valid_pts_thres, ht_range_20dB, cnt_threshold, kaz_ge_radar_const_dB, ...
                     kaz_ge_Z_near_field_correction, output_filename, kaz_ge_time_stamp, kaz_ge_range);
                  
               end % end if(good_dcl_mom_filename == 2)
               
               %% Do 3x3 filter of the velocities
               
               %                % generate a valid_obs matrix: 1=valid, NaN=not good.
               %                [m,n]    = size(spc_dcl_Vmean);
               %                spc_dcl_valid_obs    = ones(m,n);
               %                spc_dcl_valid_obs(isnan(spc_dcl_Vmean)) = spc_dcl_valid_obs(isnan(spc_dcl_Vmean)) .* NaN;
               %
               %                % Remove isolated pixels (need 3 of the 8 neighbors in 3x3)
               %                [spc_dcl_valid_obs] = func_remove_isolated_pixels_3x3(spc_dcl_valid_obs);
               %                %[spc_dcl_valid_obs_6_5x5] = func_remove_isolated_pixels_5x5(spc_dcl_valid_obs);
               %
               %                % At this point, spc_dcl_valid_obs_3x3 identifies valid profiles
               %
               %                % The shifting logic requires that the spc_dcl_Vmean be valid
               %                % to include in the shifted-averaged spectrum. So, NaN-out
               %                % spc_dcl_Vmean where ever spc_dcl_valid_obs = NaN;
               %                spc_dcl_Vmean(isnan(spc_dcl_valid_obs)) = spc_dcl_Vmean(isnan(spc_dcl_valid_obs)) .* NaN;
               
               %% Save the valid obs

               % append this valid obs variable to the raw kazr moments
               % kaz_ge_single_valid_obs    = spc_dcl_valid_obs;               
               % save(output_filename,'kaz_ge_single_valid_obs','-append');
                
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
               [~] = func_ge_ave_spc_to_15sec_and_calc_mom_and_save_spc(spc_orig_time_lin, ...
                  kaz_ge_orig_time, spc_dcl_Vmean, hour_loop, ...
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
         
         max_ht      = 4;
         hourly_flag = 1;
         
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

