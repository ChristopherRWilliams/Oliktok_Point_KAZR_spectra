function [dominant_moments, sub_peak_moments, ...
   separate_peak_moments, noise_stats] = ...
   func_find_multi_mom_3spc_Vmean_prior(spc_ref, Vd, Nspc, valid_pts_thres, ...
   index_orig, Vmean_prior, ht_ref_dB, radar_const_dB, ...
   sub_peak_max_num, separate_peak_max_num)

% Inputs:
% spc_ref     - reference spectrum
% Vd          - Doppler velocity bins
% Nspc        - number of incoherent spectral averages
% index_orig  - first and last index of original spc
% Vmean_prior - Vmean as first guess for location of this Vmean
% ht_ref_dB   - ht expressed as 20log10(ht)
% radar_const_dB - calibration constant
% sub_peak_max_num   - max number of sub-peaks to find
% separate_peak_max_num    - max number of separte peaks to find

% Outputs:
% dominant_moments   -
% sub_peak_moments   - 
% separate_peak_moments    - 
% noise_stats              - mean, max and std of noise

% Steps:
%  0. Set the output variables to NaNs
%  1. Work with spc_ref
%  2. Find mean and max noise
%  3. Weight spectrum to isolate downward peak
%  4. Find single peak
%  5. Calculate the TDA correction
%  6. Do the TDA correction
%  7. Find moments

% updated: 21-September-2017
% ========================================================================

%% 0. Set the output variables to NaNs

noise_stats    = [NaN NaN NaN];

% single peak indices
dominant_moments     = ones(1,17) .* NaN;
dominant_index_limit = ones(1,2) .* NaN;

% set the number of multiple subpeaks (of dominant peak)
%sub_peak_max_num       = 4;
%sub_peak_num           = NaN;
sub_peak_index_limit    = ones(sub_peak_max_num,2) .* NaN;
sub_peak_moments        = ones(sub_peak_max_num,17) .* NaN;

%separate_peak_max_num     = 4;
%separate_peak_num         = NaN;
separate_peak_index_limit  = ones(separate_peak_max_num,2) .* NaN;
separate_peak_moments      = ones(separate_peak_max_num,17) .* NaN;

%spc_declutter_dB = spc_ref .* NaN;

% define the drop thresholds
%drop_1st_neighbor_thres    = 3;
%drop_2nd_neighbor_thres    = 9;

% Define the number of consecutive points needed for valid spectrum
%valid_pts_thres            = 5;
%valid_pts_thres            = 3;

%%  1. rename the spectrum variable to keep the original safe

spc = spc_ref;

%%  2. Find noise statistics: mean, max and std

f = ~isnan(spc);
if(sum(f) > 0)
      
   % Get the original centeral spectrum
   spc_orig    = spc(index_orig(1):index_orig(2));
   
   % Perform the HS noise search
   initial_seed_length = length(spc_orig)/8;
   [nos_mean, nos_max] = func_find_mean_HS_noise(spc_orig, Nspc, initial_seed_length);

   % find the noise standard deviation
   f  = spc_orig < nos_max;
   if(sum(f) > 3)
      nos_std  = std(spc_orig(f));
   else
      nos_std  = NaN;
   end % end if(sum(f) > 3)
   
   % Set the noise threshold
   nos_threshold = nos_max;

   % save the noise stats
   noise_stats(1)    = nos_mean;
   noise_stats(2)    = nos_max;
   noise_stats(3)    = nos_std;
      
else
   
   % If not enough valid points, then exit
   return
end % end if(sum(f) > 0)

%%    3. Determine if enough points are above the noise thres.

% Determine if any points above the max noise
f = spc > nos_threshold;

% Exit this routine if there aren't enough points above the noise
if(sum(f) < valid_pts_thres)
   return
end % end if(sum(f) < 2)

%%    4. Use the Vmean_Prior to isolate the Vmean

if(~isnan(Vmean_prior))
   [~, peak_index]   = min(abs(Vd - Vmean_prior));
else
   [~, peak_index]   = min(abs(Vd - 0));
end % end if(~isnan(Vmean_prior))

% keep only the middle portion
left_index  = peak_index - (length(Vd)/3)/2;
right_index = peak_index + (length(Vd)/3)/2;
end_index   = length(spc);

%spc_center = spc;

if(left_index >= 1)
   spc(1:left_index) = ones(length(1:left_index),1) .* nos_mean;
end % end if(left_index >= 1)

if(right_index <= length(spc))
   spc(right_index:end_index) = ones(length(right_index:end_index),1) .* nos_mean;
end % end if(left_index >= 1)

%%  5. Find single peak


% Need to set points below the mean_noise floor to NaN
f              = spc < nos_threshold;
spc(f)         = spc(f) .* NaN;

% save the spc for future moment estimation
spc_multi_peak = spc;   % in linear units

% find a single point by setting valley threshold to a very large number
valley_thres   = 100;
[~, index_limit] = func_find_single_peak_Vmean_prior_valley(spc_multi_peak, valley_thres);
   
%%    6. Are there enough points in this single peak to keep it?

% how many points are above the noise threshold?
dominant_num_valid_pts  = index_limit(2) - index_limit(1) + 1;

% Set flags to zero
dominant_peak_flag         = 0;
find_multiple_peak_flag    = 0;
%valid_pts_thres            = 5;

if(dominant_num_valid_pts >= valid_pts_thres)
   
   dominant_index_limit(1)   = index_limit(1);
   dominant_index_limit(2)   = index_limit(2);
   
   % set flag to look for multiple peaks
   find_multiple_peak_flag = 1;
   
   % set the flag to calculate moments for this peak
   dominant_peak_flag = 1;
   
end % end if(num_valid_pts >= valid_pts_thres)

%%    7. Find multiple peaks

% Set the valley (in dB)
valley_thres   = 6;

% Get the spectra into a new variable.
% After finding a peak, those points are set to NaN in this variable
spc_multi_peak    = spc;   % in linear units
%spc_wt_multi_peak = spc_wt;   % in linear units

% set the current_peak_num
sub_peak_num        = 0;
separate_peak_num   = 0;

while(find_multiple_peak_flag)
   
   % Find a peak
   % ===========
   [~, index_limit] = func_find_single_peak_Vmean_prior_valley(spc_multi_peak, valley_thres);
      
   % make sure there are enough points in spectrum
   num_valid_pts  = index_limit(2) - index_limit(1) + 1;
   
   % Process if there are enough points in spectrum
   if(num_valid_pts >= valid_pts_thres)
      
      % Check to see if this peak is the dominant peak.
      % ===============================================
      % If yes, then do not save it as a sub peak.
      same_lt_index = index_limit(1) == dominant_index_limit(1);
      same_rt_index = index_limit(2) == dominant_index_limit(2);
      
      if(same_lt_index && same_rt_index)
         
         % do not save this peak.
         % set these values to NaN in the spectrum
         spc_multi_peak(index_limit(1):index_limit(2))      =    spc_multi_peak(index_limit(1):index_limit(2)) .* NaN;
         %spc_wt_multi_peak(index_limit(1):index_limit(2))   = spc_wt_multi_peak(index_limit(1):index_limit(2)) .* NaN;
         
      else
         % This is either a sub-peak or a separate peak
         % ============================================
         
         % Sub-Peak: Is this peak within the dominant peak domain?
         % =======================================================
         good_lt_index = index_limit(1) >= dominant_index_limit(1);
         good_rt_index = index_limit(2) <= dominant_index_limit(2);
         
         if(good_lt_index && good_rt_index)
            
            % this peak is a sub peak of the dominant peak
            sub_peak_num  = sub_peak_num + 1;
            
            % save the indices
            sub_peak_index_limit(sub_peak_num,1)  = index_limit(1);
            sub_peak_index_limit(sub_peak_num,2)  = index_limit(2);
            
            % set these values to NaN in the spectrum
            spc_multi_peak(index_limit(1):index_limit(2))      =    spc_multi_peak(index_limit(1):index_limit(2)) .* NaN;
            %spc_wt_multi_peak(index_limit(1):index_limit(2))   = spc_wt_multi_peak(index_limit(1):index_limit(2)) .* NaN;
            
         else
            % Separte Peak: if we get here, then peak is not within dominant peak
            % ============================================================
            
            % this peak is a separate peak outside of the dominant peak
            separate_peak_num  = separate_peak_num + 1;
            
            % save the indices
            separate_peak_index_limit(separate_peak_num,1)  = index_limit(1);
            separate_peak_index_limit(separate_peak_num,2)  = index_limit(2);
            
            % set these values to NaN in the spectrum
            spc_multi_peak(index_limit(1):index_limit(2))      =    spc_multi_peak(index_limit(1):index_limit(2)) .* NaN;
            %spc_wt_multi_peak(index_limit(1):index_limit(2))   = spc_wt_multi_peak(index_limit(1):index_limit(2)) .* NaN;
            
         end % end if(good_lt_index && good_rt_index)
         
      end % end if(same_lt_index && same_rt_index)
        
   else
         
      % If we get here, there are not enough points in peak to keep it.
      % ===============================================================
         
      % stop looking for peaks: set flag to zero
      find_multiple_peak_flag = 0;
         
   end % end if(num_valid_pts >= valid_pts_thres)
      
   % To prevent infinite loop, set limit on number of multiple peaks
   if((sub_peak_num >= sub_peak_max_num) || (separate_peak_num >= separate_peak_max_num))
      
      % stop looking for peaks: set flag to zero
      find_multiple_peak_flag = 0;
      
   end % end if(current_peak >= max_num_multi_peaks)
   
end % end while(find_multiple_peak_flag)

%%    8. Find the moments for each peak

% At this point, there are three different types of peaks:
% 1. Dominant peak (single peak method)
% 2. Sub-peaks of dominant peak
% 3. Separate peaks that have noise between them and the dominant peak

% process the dominant peak if dominant_peak_process_flag is set
% process sub-peaks if sub_peaks_num > 0
% process separate peaks if separate_peaks_num > 0

% Process if there is a valid dominant peak
if(dominant_peak_flag)
   
   % Dominant peak processing:
   % ========================
   
   % Construct short arrays
   spk_short = spc(dominant_index_limit(1):dominant_index_limit(2));
   Vd_short  = Vd(dominant_index_limit(1):dominant_index_limit(2));
   
   Npts        = length(Vd)/3;
   
   % get the moments
   [snr_out, nos_out, zdb_out, Vmean_out, Vsig_out, ...
      Vskew_out, Vkurt_out, Vpeak_out, Ppeak_out, ...
      Vslope_lt_out, Vslope_rt_out, ...
      std_1bin_out, std_2bin_out, ...
      max_1bin_delta_out, max_2bin_delta_out] = ...
      func_find_moments(spk_short, Vd_short, nos_mean, Npts, ...
      ht_ref_dB, radar_const_dB);
   
   dominant_moments(1)  = snr_out;
   dominant_moments(2)  = nos_out;
   dominant_moments(3)  = zdb_out;
   dominant_moments(4)  = Vmean_out;
   dominant_moments(5)  = Vsig_out;
   dominant_moments(6)  = Vskew_out;
   dominant_moments(7)  = Vkurt_out;
   dominant_moments(8)  = Vpeak_out;
   dominant_moments(9)  = Ppeak_out;
   dominant_moments(10) = Vslope_lt_out;
   dominant_moments(11) = Vslope_rt_out;
   dominant_moments(12) = dominant_index_limit(1);
   dominant_moments(13) = dominant_index_limit(2);
   dominant_moments(14) = std_1bin_out;
   dominant_moments(15) = std_2bin_out;
   dominant_moments(16) = max_1bin_delta_out;
   dominant_moments(17) = max_2bin_delta_out;
   
   % Sub-peak processing
   % ===================
   if(sub_peak_num > 0)
      
      sub_peak_moments  = ones(sub_peak_num,10) .* NaN;
      
      for k = 1:sub_peak_num
         
         % Construct short arrays
         spk_short = spc(sub_peak_index_limit(k,1):sub_peak_index_limit(k,2));
         Vd_short  =  Vd(sub_peak_index_limit(k,1):sub_peak_index_limit(k,2));
         
         Npts        = length(Vd)/3;
         
         % get the moments
         [snr_out, nos_out, zdb_out, Vmean_out, Vsig_out, ...
            Vskew_out, Vkurt_out, Vpeak_out, Ppeak_out, ...
            Vslope_lt_out, Vslope_rt_out, ...
            std_1bin_out, std_2bin_out, ...
            max_1bin_delta_out, max_2bin_delta_out] = ...
            func_find_moments(spk_short, Vd_short, nos_mean, Npts, ...
            ht_ref_dB, radar_const_dB);
         
         sub_peak_moments(k,1)   = snr_out;
         sub_peak_moments(k,2)   = nos_out;
         sub_peak_moments(k,3)   = zdb_out;
         sub_peak_moments(k,4)   = Vmean_out;
         sub_peak_moments(k,5)   = Vsig_out;
         sub_peak_moments(k,6)   = Vskew_out;
         sub_peak_moments(k,7)   = Vkurt_out;
         sub_peak_moments(k,8)   = Vpeak_out;
         sub_peak_moments(k,9)   = Ppeak_out;
         sub_peak_moments(k,10)  = Vslope_lt_out;
         sub_peak_moments(k,11)  = Vslope_rt_out;
         sub_peak_moments(k,12)  = sub_peak_index_limit(k,1);
         sub_peak_moments(k,13)  = sub_peak_index_limit(k,2);
         sub_peak_moments(k,14)  = std_1bin_out;
         sub_peak_moments(k,15)  = std_2bin_out;
         sub_peak_moments(k,16)  = max_1bin_delta_out;
         sub_peak_moments(k,17)  = max_2bin_delta_out;
         
      end % end for k loop
   end % end if(sub_peak_num > 0)
   
   % Separate peak processing
   % ========================
   if(separate_peak_num > 0)
      
      separate_peak_moments  = ones(separate_peak_num,10) .* NaN;
      
      for k = 1:separate_peak_num
         
         % Construct short arrays
         spk_short = spc(separate_peak_index_limit(k,1):separate_peak_index_limit(k,2));
         Vd_short  =  Vd(separate_peak_index_limit(k,1):separate_peak_index_limit(k,2));
         
         Npts        = length(Vd)/3;
         
         % get the moments
         [snr_out, nos_out, zdb_out, Vmean_out, Vsig_out, ...
            Vskew_out, Vkurt_out, Vpeak_out, Ppeak_out, ...
            Vslope_lt_out, Vslope_rt_out, ...
            std_1bin_out, std_2bin_out, ...
            max_1bin_delta_out, max_2bin_delta_out] = ...
            func_find_moments(spk_short, Vd_short, nos_mean, Npts, ...
            ht_ref_dB, radar_const_dB);
         
         separate_peak_moments(k,1)   = snr_out;
         separate_peak_moments(k,2)   = nos_out;
         separate_peak_moments(k,3)   = zdb_out;
         separate_peak_moments(k,4)   = Vmean_out;
         separate_peak_moments(k,5)   = Vsig_out;
         separate_peak_moments(k,6)   = Vskew_out;
         separate_peak_moments(k,7)   = Vkurt_out;
         separate_peak_moments(k,8)   = Vpeak_out;
         separate_peak_moments(k,9)   = Ppeak_out;
         separate_peak_moments(k,10)  = Vslope_lt_out;
         separate_peak_moments(k,11)  = Vslope_rt_out;
         separate_peak_moments(k,12)  = separate_peak_index_limit(k,1);
         separate_peak_moments(k,13)  = separate_peak_index_limit(k,2);
         separate_peak_moments(k,14)  = std_1bin_out;
         separate_peak_moments(k,15)  = std_2bin_out;
         separate_peak_moments(k,16)  = max_1bin_delta_out;
         separate_peak_moments(k,17)  = max_2bin_delta_out;
         
      end % end for k loop
   end % end if(separate_peak_num > 0)
     
end % end if(dominant_peak_flag)
  