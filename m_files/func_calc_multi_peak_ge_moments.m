function [kaz_ge_single_snr, kaz_ge_single_Vmean, kaz_ge_single_Vpeak] = ...
   func_calc_multi_peak_ge_moments(spc_input, cnt_spc_input, ...
   Vd_orig, Nspc, valid_pts_thres, ht_range_20dB, cnt_threshold, kaz_ge_radar_const_dB, ...
   kaz_ge_Z_near_field_correction, output_filename, kaz_ge_time_stamp, kaz_ge_range)

% This function calculates the moments for this input spectrum
% matrix

% updated: 29-July-2017
% =======================================================================

%% Construct the expaned velocity bins

dVd         = Vd_orig(6) - Vd_orig(5);

% The Nyquist is in velocity_bins(1);
% velocity_bins(1)        =  -5.9634
% velocity_bins(2)        =  -5.9168
% velocity_bins(npts/2)   =  -0.0466
% velocity_bins(npts/2+1) =   0
% velocity_bins(npts)     =  5.9168
% dV                      =  0.0466

% Expand the velocity_bins
Vd          = [Vd_orig ; Vd_orig ; Vd_orig];

new_npts = length(Vd);
% The DC point is in location (new_npts/2 + 1)
% Relabel the velocity bins:
for c = (new_npts/2):-1:1
   Vd(c) = Vd(c+1) - dVd;
end % end for c loop

for c = (new_npts/2+1+1):new_npts
   Vd(c) = Vd(c-1) + dVd;
end % end for c loop

%% Process each minute of observations

% Define some parameters
npts           = length(Vd_orig);
index_orig(1)  = npts+1;
index_orig(2)  = index_orig(1) + npts - 1;

%% Define the output variables

% Pre-define the moment matrices
[m, nhts, ~]         = size(spc_input);
%nhts                    = length(ht_range);

sub_peak_max_num          = 4;
separate_peak_max_num     = 4;

% define the dominant, single peak
kaz_ge_single_snr            = ones(m,nhts) .* NaN;
kaz_ge_single_snr_adj        = ones(m,nhts) .* NaN;
kaz_ge_single_zdb_arm        = ones(m,nhts) .* NaN;
kaz_ge_single_zdb_arm_near   = ones(m,nhts) .* NaN;
kaz_ge_single_zdb_crw        = ones(m,nhts) .* NaN;
kaz_ge_single_nos            = ones(m,nhts) .* NaN;
kaz_ge_single_nos_adj        = ones(m,nhts) .* NaN;
kaz_ge_single_nos_profile    = ones(m,1) .* NaN;
kaz_ge_single_Vmean          = ones(m,nhts) .* NaN;
kaz_ge_single_Vsig           = ones(m,nhts) .* NaN;
kaz_ge_single_Vskew          = ones(m,nhts) .* NaN;
kaz_ge_single_Vkurt          = ones(m,nhts) .* NaN;
kaz_ge_single_Vpeak          = ones(m,nhts) .* NaN;
kaz_ge_single_Ppeak          = ones(m,nhts) .* NaN;
kaz_ge_single_Vslope_lt      = ones(m,nhts) .* NaN;
kaz_ge_single_Vslope_rt      = ones(m,nhts) .* NaN;
kaz_ge_single_index_lt       = ones(m,nhts) .* NaN;
kaz_ge_single_index_rt       = ones(m,nhts) .* NaN;
kaz_ge_single_Vleft          = ones(m,nhts) .* NaN;
kaz_ge_single_Vright         = ones(m,nhts) .* NaN;
kaz_ge_single_nos_mean       = ones(m,nhts) .* NaN;
kaz_ge_single_nos_max        = ones(m,nhts) .* NaN;
kaz_ge_single_nos_std        = ones(m,nhts) .* NaN;
kaz_ge_single_std_1bin       = ones(m,nhts) .* NaN;
kaz_ge_single_std_2bin       = ones(m,nhts) .* NaN;
kaz_ge_single_max_1bin_delta = ones(m,nhts) .* NaN;
kaz_ge_single_max_2bin_delta = ones(m,nhts) .* NaN;

% define the sub peaks
kaz_ge_sub_snr            = ones(m,nhts,sub_peak_max_num) .* NaN;
%kaz_ge_sub_snr_adj        = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_zdb_arm        = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_zdb_crw        = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_nos            = ones(m,nhts,sub_peak_max_num) .* NaN;
%kaz_ge_sub_nos_adj        = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vmean          = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vsig           = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vskew          = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vkurt          = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vpeak          = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Ppeak          = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vslope_lt      = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vslope_rt      = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_index_lt       = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_index_rt       = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vleft          = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_Vright         = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_std_1bin       = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_std_2bin       = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_max_1bin_delta = ones(m,nhts,sub_peak_max_num) .* NaN;
kaz_ge_sub_max_2bin_delta = ones(m,nhts,sub_peak_max_num) .* NaN;

% define the separate peaks
kaz_ge_separate_snr            = ones(m,nhts,separate_peak_max_num) .* NaN;
%kaz_ge_separate_snr_adj        = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_zdb_arm        = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_zdb_crw        = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_nos            = ones(m,nhts,separate_peak_max_num) .* NaN;
%kaz_ge_separate_nos_adj        = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vmean          = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vsig           = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vskew          = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vkurt          = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vpeak          = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Ppeak          = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vslope_lt      = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vslope_rt      = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_index_lt       = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_index_rt       = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vleft          = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_Vright         = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_std_1bin       = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_std_2bin       = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_max_1bin_delta = ones(m,nhts,separate_peak_max_num) .* NaN;
kaz_ge_separate_max_2bin_delta = ones(m,nhts,separate_peak_max_num) .* NaN;

%disp('paused...')
%pause

%% Process each time interval

for prof_num = 1:m
   %for prof_num = 1:1
   
   %% display a message
   if(rem(prof_num,60) == 0)
      disp(['...processing profile #',num2str(prof_num),' out of ',num2str(m),' profiles...']);
   end % end if(rem(prof_num,60) == 0)
   
   %% Calculate the Vmean of averaged spectrum
   
   % Find moments of averaged spectrum
   % ==================================
   % Start from the top, and find the moments of each spectrum
   % Use the previous above estimate for initial guess.
   
   % set the Vmean_prior to NaN
   Vmean_prior    = NaN;
   
   % Process each range gate, starting from the top
   for c = nhts:-1:1
      
      spc_lin        = make_column_vector(squeeze(spc_input(prof_num,c,:)));
      
      % Expand the spectra to correct for aliasing
      spc_ref        = [spc_lin ; spc_lin; spc_lin];
      
      % get the range in 20dB for reflectivity calculations
      ht_ref_dB      = ht_range_20dB(c);
      
      % Calculate the moments for this spectrum
      if(cnt_spc_input(prof_num,c) >= cnt_threshold)
         
         total_Nspc  = Nspc * cnt_spc_input(prof_num,c);
         
         [dominant_moments, sub_peak_moments, ...
            separate_peak_moments, noise_stats] = ...
            func_find_multi_mom_3spc_Vmean_prior(spc_ref, Vd, total_Nspc, valid_pts_thres, ...
            index_orig, Vmean_prior, ht_ref_dB, kaz_ge_radar_const_dB, ...
            sub_peak_max_num, separate_peak_max_num);

         %if(c == 77)
         %   disp('paused.')
         %   pause
         %end % end if(c == 77)
         
         % Format is:
         %                   dominant_moments(1)  = snr_out;
         %                   dominant_moments(2)  = nos_out;
         %                   dominant_moments(3)  = zdb_out;
         %                   dominant_moments(4)  = Vmean_out;
         %                   dominant_moments(5)  = Vsig_out;
         %                   dominant_moments(6)  = Vskew_out;
         %                   dominant_moments(7)  = Vkurt_out;
         %                   dominant_moments(8)  = Vpeak_out;
         %                   dominant_moments(9)  = Ppeak_out;
         %                   dominant_moments(10) = Vslope_lt_out;
         %                   dominant_moments(11) = Vslope_rt_out;
         %                   dominant_moments(12) = dominant_index_limit(k,1);
         %                   dominant_moments(13) = dominant_index_limit(k,2);
         %                   dominant_moments(14) = std_1bin_out;
         %                   dominant_moments(15) = std_2bin_out;
         %                   dominant_moments(16) = max_1bin_delta_out;
         %                   dominant_moments(17) = max_2bin_delta_out;
         %
         
         % Save the noise stats
         kaz_ge_single_nos_mean(prof_num,c)       = noise_stats(1);
         kaz_ge_single_nos_max(prof_num,c)        = noise_stats(2);
         kaz_ge_single_nos_std(prof_num,c)        = noise_stats(3);
         
         % Save the dominant peak moments (if any)
         % =======================================
         if(~isnan(dominant_moments(1)))
            kaz_ge_single_snr(prof_num,c)            = dominant_moments(1);
            kaz_ge_single_nos(prof_num,c)            = dominant_moments(2);
            kaz_ge_single_zdb_arm(prof_num,c)        = dominant_moments(3);
            kaz_ge_single_Vmean(prof_num,c)          = dominant_moments(4);
            kaz_ge_single_Vsig(prof_num,c)           = dominant_moments(5);
            kaz_ge_single_Vskew(prof_num,c)          = dominant_moments(6);
            kaz_ge_single_Vkurt(prof_num,c)          = dominant_moments(7);
            kaz_ge_single_Vpeak(prof_num,c)          = dominant_moments(8);
            kaz_ge_single_Ppeak(prof_num,c)          = dominant_moments(9);
            kaz_ge_single_Vslope_lt(prof_num,c)      = dominant_moments(10);
            kaz_ge_single_Vslope_rt(prof_num,c)      = dominant_moments(11);
            kaz_ge_single_index_lt(prof_num,c)       = dominant_moments(12);
            kaz_ge_single_index_rt(prof_num,c)       = dominant_moments(13);
            kaz_ge_single_std_1bin(prof_num,c)       = dominant_moments(14);
            kaz_ge_single_std_2bin(prof_num,c)       = dominant_moments(15);
            kaz_ge_single_max_1bin_delta(prof_num,c) = dominant_moments(16);
            kaz_ge_single_max_2bin_delta(prof_num,c) = dominant_moments(17);
            
            kaz_ge_single_Vleft(prof_num,c)          = Vd(dominant_moments(12));
            kaz_ge_single_Vright(prof_num,c)         = Vd(dominant_moments(13));
            
            % set the Vmean_prior to a valid number
            if(kaz_ge_single_snr(prof_num,c) > -10)
               Vmean_prior    = dominant_moments(4);
            end % end if(kaz_ge_snr(prof_num,c) > -10)
            
         end % end if(~isnan(dominant_moments(1)))
                  
         % Save the sub peak moments (if any)
         % ==================================
         [mm,~]    = size(sub_peak_moments);
         
         for k = 1:mm
            
            if(~isnan(sub_peak_moments(k,1)))
               kaz_ge_sub_snr(prof_num,c,k)          = sub_peak_moments(k,1);
               kaz_ge_sub_nos(prof_num,c,k)          = sub_peak_moments(k,2);
               kaz_ge_sub_zdb_arm(prof_num,c,k)      = sub_peak_moments(k,3);
               kaz_ge_sub_Vmean(prof_num,c,k)        = sub_peak_moments(k,4);
               kaz_ge_sub_Vsig(prof_num,c,k)         = sub_peak_moments(k,5);
               kaz_ge_sub_Vskew(prof_num,c,k)        = sub_peak_moments(k,6);
               kaz_ge_sub_Vkurt(prof_num,c,k)        = sub_peak_moments(k,7);
               kaz_ge_sub_Vpeak(prof_num,c,k)        = sub_peak_moments(k,8);
               kaz_ge_sub_Ppeak(prof_num,c,k)        = sub_peak_moments(k,9);
               kaz_ge_sub_Vslope_lt(prof_num,c,k)    = sub_peak_moments(k,10);
               kaz_ge_sub_Vslope_rt(prof_num,c,k)    = sub_peak_moments(k,11);
               kaz_ge_sub_index_lt(prof_num,c,k)     = sub_peak_moments(k,12);
               kaz_ge_sub_index_rt(prof_num,c,k)     = sub_peak_moments(k,13);
               kaz_ge_sub_std_1bin(prof_num,c)       = sub_peak_moments(14);
               kaz_ge_sub_std_2bin(prof_num,c)       = sub_peak_moments(15);
               kaz_ge_sub_max_1bin_delta(prof_num,c) = sub_peak_moments(16);
               kaz_ge_sub_max_2bin_delta(prof_num,c) = sub_peak_moments(17);
               kaz_ge_sub_Vleft(prof_num,c,k)        = Vd(sub_peak_moments(k,12));
               kaz_ge_sub_Vright(prof_num,c,k)       = Vd(sub_peak_moments(k,13));
               
            end % end if(~isnan(sub_peak_moments(k,1)))
         end % end for k loop
         
         % Save the separate peak moments (if any)
         % ==================================
         [mm,~]    = size(separate_peak_moments);
         
         for k = 1:mm
            
            if(~isnan(separate_peak_moments(k,1)))
               kaz_ge_separate_snr(prof_num,c,k)           = separate_peak_moments(k,1);
               kaz_ge_separate_nos(prof_num,c,k)           = separate_peak_moments(k,2);
               kaz_ge_separate_zdb_arm(prof_num,c,k)       = separate_peak_moments(k,3);
               kaz_ge_separate_Vmean(prof_num,c,k)         = separate_peak_moments(k,4);
               kaz_ge_separate_Vsig(prof_num,c,k)          = separate_peak_moments(k,5);
               kaz_ge_separate_Vskew(prof_num,c,k)         = separate_peak_moments(k,6);
               kaz_ge_separate_Vkurt(prof_num,c,k)         = separate_peak_moments(k,7);
               kaz_ge_separate_Vpeak(prof_num,c,k)         = separate_peak_moments(k,8);
               kaz_ge_separate_Ppeak(prof_num,c,k)         = separate_peak_moments(k,9);
               kaz_ge_separate_Vslope_lt(prof_num,c,k)     = separate_peak_moments(k,10);
               kaz_ge_separate_Vslope_rt(prof_num,c,k)     = separate_peak_moments(k,11);
               kaz_ge_separate_index_lt(prof_num,c,k)      = separate_peak_moments(k,12);
               kaz_ge_separate_index_rt(prof_num,c,k)      = separate_peak_moments(k,13);
               kaz_ge_separate_std_1bin(prof_num,c)        = separate_peak_moments(14);
               kaz_ge_separate_std_2bin(prof_num,c)        = separate_peak_moments(15);
               kaz_ge_separate_max_1bin_delta(prof_num,c)  = separate_peak_moments(16);
               kaz_ge_separate_max_2bin_delta(prof_num,c)  = separate_peak_moments(17);
               kaz_ge_separate_Vleft(prof_num,c,k)         = Vd(separate_peak_moments(k,12));
               kaz_ge_separate_Vright(prof_num,c,k)        = Vd(separate_peak_moments(k,13));
               
            end % end if(~isnan(separate_peak_moments(k,1)))
         end % end for k loop
         
      end % end if(cnt_spc_input(prof_num,c) > 1)
      
   end % end for c loop
   
end % end for prof_num loop

%% Now that you have the profile, do the noise adjustment

for prof_num = 1:m
   
   % Adjust the snr by the noise & calculate reflectivity
   % Dominant peak
   % Average profile
   % ======================================================
   snr_profile    = kaz_ge_single_snr(prof_num,:);
   nos_profile    = kaz_ge_single_nos(prof_num,:);
   ht_profile     = kaz_ge_range;
   
   [zdb_adj, snr_adj, nos_adj, new_nos_power] = func_find_noise_adjusted_zdb_and_snr(snr_profile, nos_profile, ht_profile, kaz_ge_radar_const_dB);
   
   % do the near-field correction
   kaz_zdb  = zdb_adj + abs(kaz_ge_Z_near_field_correction);
   
   kaz_ge_single_zdb_arm_near(prof_num,:) = kaz_ge_single_zdb_arm(prof_num,:) + abs(kaz_ge_Z_near_field_correction);
   
   % Save these values
   kaz_ge_single_snr_adj(prof_num,:)      = snr_adj;
   kaz_ge_single_nos_adj(prof_num,:)      = nos_adj;
   kaz_ge_single_nos_profile(prof_num,1)  = new_nos_power;
   kaz_ge_single_zdb_crw(prof_num,:)      = kaz_zdb;
   
   % Adjust the snr by the noise & calculate reflectivity
   % sub-peak
   % ======================================================
   for k = 1:sub_peak_max_num
      snr_profile    = kaz_ge_sub_snr(prof_num,:,k);
      
      [zdb_adj, ~, ~, ~] = ...
         func_find_noise_adjusted_zdb_and_snr(snr_profile, nos_profile, ht_profile, kaz_ge_radar_const_dB);
      
      % do the near-field correction
      kaz_ge_sub_zdb_crw(prof_num,:,k) = zdb_adj + abs(kaz_ge_Z_near_field_correction);
      
   end % end for k loop
   
   % Adjust the snr by the noise & calculate reflectivity
   % separate peak
   % ======================================================
   for k = 1:separate_peak_max_num
      snr_profile    = kaz_ge_separate_snr(prof_num,:,k);
      
      [zdb_adj, ~, ~, ~] = ...
         func_find_noise_adjusted_zdb_and_snr(snr_profile, nos_profile, ht_profile, kaz_ge_radar_const_dB);
      
      % do the near-field correction
      kaz_ge_separate_zdb_crw(prof_num,:,k) = zdb_adj + abs(kaz_ge_Z_near_field_correction);
      
   end % end for k loop
      
end % end for r loop

%% Save the decluttered moments

save(output_filename,'kaz_ge_time_stamp',...
   'kaz_ge_range',...
   'kaz_ge_single_snr','kaz_ge_single_snr_adj',...
   'kaz_ge_single_zdb_crw','kaz_ge_single_zdb_arm',...
   'kaz_ge_single_zdb_arm_near',...
   'kaz_ge_single_nos','kaz_ge_single_nos_adj',...
   'kaz_ge_single_nos_profile',...
   'kaz_ge_single_Vmean','kaz_ge_single_Vsig',...
   'kaz_ge_single_Vskew','kaz_ge_single_Vkurt',...
   'kaz_ge_single_Vpeak','kaz_ge_single_Ppeak',...
   'kaz_ge_single_Vslope_lt','kaz_ge_single_Vslope_rt',...
   'kaz_ge_single_index_lt','kaz_ge_single_index_rt',...
   'kaz_ge_single_Vleft','kaz_ge_single_Vright',...
   'kaz_ge_single_nos_mean','kaz_ge_single_nos_max',...
   'kaz_ge_single_nos_std',...
   'kaz_ge_single_std_1bin','kaz_ge_single_std_2bin',...
   'kaz_ge_single_max_1bin_delta','kaz_ge_single_max_2bin_delta',...
   'kaz_ge_sub_snr','kaz_ge_sub_nos',...
   'kaz_ge_sub_zdb_crw','kaz_ge_sub_zdb_arm',...
   'kaz_ge_sub_Vmean','kaz_ge_sub_Vsig',...
   'kaz_ge_sub_Vskew','kaz_ge_sub_Vkurt',...
   'kaz_ge_sub_Vpeak','kaz_ge_sub_Ppeak',...
   'kaz_ge_sub_Vslope_lt','kaz_ge_sub_Vslope_rt',...
   'kaz_ge_sub_index_lt','kaz_ge_sub_index_rt',...
   'kaz_ge_sub_Vleft','kaz_ge_sub_Vright',...
   'kaz_ge_sub_std_1bin','kaz_ge_sub_std_2bin',...
   'kaz_ge_sub_max_1bin_delta','kaz_ge_sub_max_2bin_delta',...
   'kaz_ge_separate_snr','kaz_ge_separate_nos',...
   'kaz_ge_separate_zdb_crw','kaz_ge_separate_zdb_arm',...
   'kaz_ge_separate_Vmean','kaz_ge_separate_Vsig',...
   'kaz_ge_separate_Vskew','kaz_ge_separate_Vkurt',...
   'kaz_ge_separate_Vpeak','kaz_ge_separate_Ppeak',...
   'kaz_ge_separate_Vslope_lt','kaz_ge_separate_Vslope_rt',...
   'kaz_ge_separate_index_lt','kaz_ge_separate_index_rt',...
   'kaz_ge_separate_Vleft','kaz_ge_separate_Vright',...
   'kaz_ge_separate_std_1bin','kaz_ge_separate_std_2bin',...
   'kaz_ge_separate_max_1bin_delta','kaz_ge_separate_max_2bin_delta',...
   'kaz_ge_radar_const_dB','kaz_ge_Z_near_field_correction');

%% Create a dumby variable for sending back to calling program

%done = 1;

