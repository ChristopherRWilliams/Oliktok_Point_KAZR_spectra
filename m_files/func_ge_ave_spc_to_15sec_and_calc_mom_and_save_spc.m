function [done] = func_ge_ave_spc_to_15sec_and_calc_mom_and_save_spc(spc_orig_time_lin, ...
   kaz_ge_orig_time, spc_dcl_Vmean, hour_loop, ...
   year_str, month_num_str, dayofmonth_str, hour_str, Vd_orig, ...
   Nspc, valid_pts_thres, ht_range_20dB, kaz_ge_radar_const_dB, ...
   kaz_ge_Z_near_field_correction, kaz_ge_range, site_prefix, version_suffex)


%% Process 15sec dwell spectra
% ============================

%% Define the 60 second incoherent averaging time stamp

delta_sec      = 15;
num_profiles   = (60*60) / delta_sec; % (number of seconds in hour)/(dwell)

time_stamp_15sec    = zeros(num_profiles,7);

time_stamp_15sec(:,1)  = ones(num_profiles,1) .* kaz_ge_orig_time(1,1); % year
time_stamp_15sec(:,2)  = ones(num_profiles,1) .* kaz_ge_orig_time(1,2); % month
time_stamp_15sec(:,3)  = ones(num_profiles,1) .* kaz_ge_orig_time(1,3); % dayofmonth
time_stamp_15sec(:,4)  = ones(num_profiles,1) .* hour_loop; % hour
time_stamp_15sec(1,5)  = 0; % minute
time_stamp_15sec(1,6)  = 0; % second

for r = 2:num_profiles
   
   time_stamp_15sec(r,5)      = time_stamp_15sec(r-1,5);
   time_stamp_15sec(r,6)      = time_stamp_15sec(r-1,6) + delta_sec;
   if(time_stamp_15sec(r,6) >= 60)
      time_stamp_15sec(r,6)   = 0;
      time_stamp_15sec(r,5)   = time_stamp_15sec(r,5) + 1;
      
   end % end if(new_sec >= 60)
   
end % end for r loop

% define the time in min.fraction
time_orig_minsec  = kaz_ge_orig_time(:,5) + ...
   kaz_ge_orig_time(:,6)./(60) + (kaz_ge_orig_time(:,7)./1000)./(60);

time_15sec_minsec   = time_stamp_15sec(:,5) + ...
   time_stamp_15sec(:,6)./(60) + (time_stamp_15sec(:,7)./1000)./(60);

%% Estimate the velocity variance during the 15 sec dwell

m              = length(time_15sec_minsec);
[~,n]          = size(spc_dcl_Vmean);

kaz_ge_single_Vmean_var       = ones(m,n) .* NaN;
kaz_ge_single_cnt_Vmean_var   = ones(m,n) .* NaN;

time_desired_minsec  = time_15sec_minsec;
   
for c = 1:n

   % Get the Vmean values
   Vmean_input    = spc_dcl_Vmean(:,c);

   [Vmean_var_output, cnt_Vmean_var_output] = ...
      func_calc_Vmean_var(Vmean_input, ...
      time_orig_minsec, time_desired_minsec);
   
   kaz_ge_single_Vmean_var(:,c)          = Vmean_var_output;
   kaz_ge_single_cnt_Vmean_var(:,c)      = cnt_Vmean_var_output;
   
end % end for c loop

%% Average spectra over 15 sec dwell

disp(' ')
disp(' Averaging spectra with 15sec dwells...');

% define the 15sec spectra
[~,n,p]           = size(spc_orig_time_lin);
m                 = length(time_15sec_minsec);
spc_15sec_lin     = ones(m,n,p) .* NaN;
cnt_spc_15sec     = ones(m,n) .* NaN;

for c = 1:n
   %for c = 10:10
   
   %if(rem(c,100) == 0)
   %   disp(['15sec, range gate #: ',num2str(c)]);
   %end % end if(rem(c,10) == 0)
   
   % Get the inputs
   spc_input            = squeeze(spc_orig_time_lin(:,c,:));
   valid_obs_input      = spc_dcl_Vmean(:,c);
   
   %Vmean_input          = squeeze(spc_orig_Vmean(:,c));
   %time_desired_minsec  = time_15sec_minsec;
   
   % Calculate the shifted and averaged spectra
   %[spc_desired_output, cnt_spc_output, Vmean_shift_output] = ...
   %   func_shift_then_ave_spc(spc_input, Vmean_input, Vd_orig, ...
   %   time_orig_minsec, time_desired_minsec);
   
   [spc_desired_output, cnt_spc_output]   = ...
      func_incoherent_ave_spc_valid_obs(spc_input, valid_obs_input, time_orig_minsec, time_15sec_minsec);
   
   spc_15sec_lin(:,c,:)    = spc_desired_output;
   cnt_spc_15sec(:,c)      = cnt_spc_output;
   
end % end for c loop

%% Estimate the moments using 15 sec dwell

disp(' ')
disp('Calculating moments of 15sec dwell spectra...')

spc_input         = spc_15sec_lin;
cnt_spc_input     = cnt_spc_15sec;
kaz_ge_time_stamp = time_stamp_15sec;
%kaz_ge_range      = kaz_ge_range;

% Determine how many samples are needed per dwell interval
max_cnt_spc    = max(max(cnt_spc_15sec));
cnt_threshold  = round(max_cnt_spc/2);

%output_directory  = [mat_drive,'\Projects\Oliktok_Point\KAZR\multi_peak_clutter_ave_spc_2017_0921\mat_15sec_ave_moments\'];
output_directory  = '..\mat_15sec_ave_moments\';
output_filename   = [output_directory,site_prefix,'_kazr_ge_15sec_mom_',year_str,month_num_str,dayofmonth_str,'_',hour_str,version_suffex,'.mat'];

[~, spc_15sec_Vmean, ~]   = func_calc_multi_peak_ge_moments(spc_input, cnt_spc_input, ...
   Vd_orig, Nspc, valid_pts_thres, ht_range_20dB, cnt_threshold, kaz_ge_radar_const_dB, ...
   kaz_ge_Z_near_field_correction, output_filename, kaz_ge_time_stamp, kaz_ge_range);

% Append the Vmean_var variable
save(output_filename,...
   'kaz_ge_single_Vmean_var','kaz_ge_single_cnt_Vmean_var',...
   '-append');

%% Shift each spectrum and then average spectra over 15 sec dwell

disp(' ')
disp(' shifting and then averaging spc over 15sec dwell...');

% define the 15sec spectra
[~,n,p]           = size(spc_orig_time_lin);
m                 = length(time_15sec_minsec);
spc_15sec_lin     = ones(m,n,p) .* NaN;
cnt_spc_15sec     = ones(m,n) .* NaN;

for c = 1:n
   %for c = 10:10
   
   %if(rem(c,100) == 0)
   %   disp(['15sec, range gate #: ',num2str(c)]);
   %end % end if(rem(c,10) == 0)
   
   % Get the inputs
   spc_input            = squeeze(spc_orig_time_lin(:,c,:));
   %valid_obs_input      = spc_dcl_valid_obs(:,c);

   Vmean_input          = squeeze(spc_dcl_Vmean(:,c));
   Vmean_desired_input  = squeeze(spc_15sec_Vmean(:,c));
   time_desired_minsec  = time_15sec_minsec;
   
   % Calculate the shifted and averaged spectra
   % The shifting requires that Vmean_input ~= NaN. So, do not need the
   % valid_obs_flag
   
   [spc_desired_output, cnt_spc_output] = ...
      func_shift_then_ave_spc(spc_input, Vmean_input,...
      Vmean_desired_input, Vd_orig, ...
      time_orig_minsec, time_desired_minsec);
   
   spc_15sec_lin(:,c,:)    = spc_desired_output;
   cnt_spc_15sec(:,c)      = cnt_spc_output;
   
end % end for c loop


%% Replace noise values with NaN

input_spc_lin  = spc_15sec_lin;
input_Nspc     = cnt_spc_15sec .* Nspc;
input_Nstd     = 2;

[kaz_ge_15sec_nos_mean, kaz_ge_15sec_nos_std, kaz_ge_15sec_nos_thres, ...
   kaz_ge_15sec_spc_lin] = ...
   func_replace_Nstd_nos_with_NaN(input_spc_lin, input_Nspc, input_Nstd); %#ok<ASGLU>

kaz_ge_15sec_spc_cnt    = input_Nspc; %#ok<NASGU>

%% Save the 15sec averaged spectra

output_directory  = '..\mat_15sec_ave_spc\';
output_filename   = [output_directory,site_prefix,'_kazr_ge_15sec_spc_',year_str,month_num_str,dayofmonth_str,'_',hour_str,version_suffex,'.mat'];

% save the spectra:
save(output_filename,'kaz_ge_15sec_nos_mean','kaz_ge_15sec_nos_std', ...
   'kaz_ge_15sec_nos_thres','kaz_ge_15sec_spc_lin','kaz_ge_15sec_spc_cnt',...
   'time_stamp_15sec','Vd_orig','kaz_ge_range');

%% Estimate the moments using 15 sec dwell

disp(' ')
disp('Calculating moments of shifted and averaged 15sec spectra...')

spc_input      = spc_15sec_lin;
cnt_spc_input  = cnt_spc_15sec;
kaz_ge_time_stamp = time_stamp_15sec;
%kaz_ge_range      = kaz_ge_range;

% Determine how many samples are needed per dwell interval
max_cnt_spc    = max(max(cnt_spc_15sec));
cnt_threshold  = round(max_cnt_spc/2);

%output_directory  = [mat_drive,'\Projects\Oliktok_Point\KAZR\multi_peak_clutter_ave_spc_2017_0921\mat_15sec_ave_moments\'];
output_directory  = '..\mat_15sec_ave_moments\';
output_filename   = [output_directory,site_prefix,'_kazr_ge_15sec_shift_mom_',year_str,month_num_str,dayofmonth_str,'_',hour_str,version_suffex,'.mat'];

[~, ~, ~]   = func_calc_multi_peak_ge_moments(spc_input, cnt_spc_input, ...
   Vd_orig, Nspc, valid_pts_thres, ht_range_20dB, cnt_threshold, kaz_ge_radar_const_dB, ...
   kaz_ge_Z_near_field_correction, output_filename, kaz_ge_time_stamp, kaz_ge_range);

% Append the Vmean_var variable
save(output_filename,...
   'kaz_ge_single_Vmean_var','kaz_ge_single_cnt_Vmean_var',...
   '-append');

%% Define output flag

done = 1;
