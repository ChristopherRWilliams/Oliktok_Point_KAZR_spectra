function [prof_noise_stats, ...
   prof_clt_peak_dB, prof_clt_drop_dB, prof_clt_cnr_dB, ...
   prof_clt_residual_snr_dB, prof_clt_residual_Vpeak,...
   prof_clt_residual_Vpeak_index, ...
   prof_clt_residual_drop_dB,...
   prof_clt_residual_num_valid_pts,...
   prof_clt_residual_lt_index, prof_clt_residual_rt_index,...
   prof_clt_residual_sacr_flag, prof_spc_dcl_lin] = ...
   func_declutter_spc(prof_spc_orig_lin, Vd, Nspc, n_clt, ...
   drop_dB_thres, valid_pts_thres)

% This function performs the de-cluttering of spectra following the logic
% of Figure 10 in Williams et al. (2018) AMT

% updated: 27-June-2018
% ========================================================================

%% Define the output variables

[~,p]                            = size(prof_spc_orig_lin);

prof_spc_dcl_lin                 = prof_spc_orig_lin;

prof_noise_stats                 = ones(n_clt,3) .* NaN;
prof_clt_peak_dB                 = ones(n_clt,1) .* NaN;
prof_clt_drop_dB                 = ones(n_clt,1) .* NaN;
prof_clt_cnr_dB                  = ones(n_clt,1) .* NaN;
prof_clt_residual_snr_dB         = ones(n_clt,1) .* NaN;
prof_clt_residual_Vpeak          = ones(n_clt,1) .* NaN;
prof_clt_residual_Vpeak_index    = ones(n_clt,1) .* NaN;
prof_clt_residual_drop_dB        = ones(n_clt,1) .* NaN;
prof_clt_residual_num_valid_pts  = ones(n_clt,1) .* NaN;
prof_clt_residual_lt_index       = ones(n_clt,1) .* NaN;
prof_clt_residual_rt_index       = ones(n_clt,1) .* NaN;
prof_clt_residual_sacr_flag      = ones(n_clt,1) .* NaN;

%% Get the indices used in the interpolation of the clutter

% find the indices for velocities within 2 velocity bins of zero
[~, DC_index]  = min(abs(Vd));

%   (DC-2)  (DC-1)  (DC)  (DC+1)  (DC+2)
% In the interpolation, valuesa at indices (DC-2) and (DC+2) do not change.

% define the interpolation indices
lt_clt_index   = DC_index - 2;
rt_clt_index   = DC_index + 2;

%% Process each range gate

% Do the decluttering for range gates 1 to n_clt

for c = 1:n_clt
   
   %% Get the spectrum from this range gate
   
   spc   = prof_spc_orig_lin(c,:);
   
   %% calculate the noise stats
   
   % set flag for future processing
   nos_threshold  = NaN;
   
   % Process spectrum if there are any valid values (not NaN)
   f = ~isnan(spc);
   if(sum(f) > 0)
      
      % Perform the HS noise search
      initial_seed_length = length(spc)/8;
      [nos_mean, nos_max] = func_find_mean_HS_noise(spc, Nspc, initial_seed_length);
      
      % find the noise standard deviation
      f  = spc < nos_max;
      if(sum(f) > 3)
         nos_std  = std(spc(f));
      end % end if(sum(f) > 3)
      
      % Set the noise threshold
      nos_threshold = nos_max;
      
      % save the noise stats
      prof_noise_stats(c,1)    = nos_mean;
      prof_noise_stats(c,2)    = nos_max;
      prof_noise_stats(c,3)    = nos_std;
      
   else
      % If not enough valid points, then exit
      % return
   end % end if(sum(f) > 0)
   
   %% If spectrum is valid, process this spectrum
   
   if(~isnan(nos_threshold))
      
      %%  Put spectrum into log units
      
      spc_dB            = 10.*log10(spc);
      spc_dcl_lin       = spc;
      nos_threshold_dB  = 10.*log10(nos_threshold);
      
      %%  Verify that there is power at DC
      
      if((spc(DC_index) > nos_threshold))
         
         % If you get here, then there could be a clutter peak.
         
         %%  Calculate the magnitude of the DC peak (relative to noise thres).
         %   (Box #2 in Fig. 10)
         
         prof_clt_peak_dB(c,1) = spc_dB(DC_index) - nos_threshold_dB;
         
         %%  Estimate the drop to 1st nearest neighbor
         
         % ensure spc(DC) is larger than both neighbors
         if((spc_dB(DC_index) > spc_dB(DC_index - 1)) && ...
               (spc_dB(DC_index) > spc_dB(DC_index + 1)))
            
            drop_lt   = spc_dB(DC_index) - spc_dB(DC_index - 1);
            drop_rt   = spc_dB(DC_index) - spc_dB(DC_index + 1);
            
            drop_1bin               = max([drop_lt drop_rt]);
            prof_clt_drop_dB(c,1)   = drop_1bin;
            
         else
            prof_clt_drop_dB(c,1)   = -1;
         end % end if spc(DC) > neighbors
         
         %% Do the process if the drop is large enough
         %  (Box #3 in Fig. 10)
         
         if(drop_1bin > drop_dB_thres)
            
            %%  Calculate the CNR for three bins around DC
            
            % The clutter that will be removed is relative to the remaining signal.
            spc_bottom  = spc(lt_clt_index:rt_clt_index);
            nn          = length(spc_bottom);
            delta_spc   = (spc_bottom(nn) - spc_bottom(1)) ./ (nn-1);
            
            % do linear interpolation between spc_bottom(1) and spc_bottom(nn)
            k     = -1;
            for ii = 1:nn
               k              = k + 1;
               spc_bottom(ii) = spc_bottom(1) + delta_spc.*k;
            end % end for ii loop
            
            % Calculate the clutter-to-noise ratio
            clutter_power           = sum(spc(lt_clt_index:rt_clt_index) - spc_bottom);
            noise_power             = nos_mean * length(spc);
            % sometimes, the area is negative (the peak value is not at DC)
            cnr_lin                 = clutter_power / noise_power;
            if(cnr_lin > 0)
               prof_clt_cnr_dB(c,1)    = (10.*log10(cnr_lin));
            end % end if(clutter_power <= 0)
            
            %%  Remove clutter from spectrum
            %   (Box #5 in Fig. 10)
            
            % replace the 3 spectral points with the interpolated values
            % (note that values at lt_clt_index and rt_clt_index are unchanged)
            
            % do a linear interpolation
            k     = -1;
            for ii = lt_clt_index:rt_clt_index
               k                 = k + 1;
               spc_dcl_lin(ii)   = spc(lt_clt_index) + delta_spc.*k;
            end % end for ii loop
            
            %%  Find the single peak moment for this decluttered spectrum
            %   (Box #6 in Fig. 10)
            
            % Verify there are enough points to estimate moments
            f  = spc_dcl_lin > nos_threshold;
            if(sum(f) >= valid_pts_thres)
               
               % put NaN into all other locations
               spc_working       = spc_dcl_lin;
               spc_working(~f)   = spc_working(~f) .* NaN;
               
               % find the indices to the single peak with largest magnitude
               [~, index_limit] = func_find_single_peak(spc_working);
               
               % save attributes of this residual peak
               % how many points are in this peak?
               prof_clt_residual_lt_index(c,1)        = index_limit(1);
               prof_clt_residual_rt_index(c,1)        = index_limit(2);
               prof_clt_residual_num_valid_pts(c,1)   = index_limit(2) - index_limit(1) + 1;
               
               % if there are enough peaks, calculate some moments
               % (Box #7 in Fig. 10)
               
               if(prof_clt_residual_num_valid_pts(c,1) >= valid_pts_thres)
                  
                  % Construct short arrays
                  spk_short = spc_working(index_limit(1):index_limit(2));
                  %Vd_short  = Vd(index_limit(1):index_limit(2));
                  
                  % calculate the snr
                  signal_power         = sum(spk_short - nos_mean);
                  noise_power          = nos_mean * length(spc);
                  prof_clt_residual_snr_dB(c,1)    = 10.*log10(signal_power / noise_power);
                  
                  % Get the peak velocity index
                  [~, max_index]  = max(spc_working);
                  prof_clt_residual_Vpeak_index(c,1)  = max_index;
                  prof_clt_residual_Vpeak(c,1)        = Vd(max_index);
                  
                  % From the peak value, estimate the power drop to nearest neighbor
                  % Check to make sure peak is not at either end of spectrum
                  % check left edge of Vd (else move index...)
                  if((max_index - 1) < 1)
                     max_index = 2;
                  end
                  % check right edge of Vd (else move index...)
                  if((max_index + 1) > p)
                     max_index = p - 1;
                  end
                  
                  % Calculate drop to left and right neighbors
                  drop_lt   = 10.*log10(spc_working(max_index)) - 10.*log10(spc_working(max_index - 1));
                  drop_rt   = 10.*log10(spc_working(max_index)) - 10.*log10(spc_working(max_index + 1));
                  
                  drop_1bin                        = max([drop_lt drop_rt]);
                  prof_clt_residual_drop_dB(c,1)   = drop_1bin;
                  
                  % Deterine if residual Vpeak is at interpolation edge
                  % (Box #10 in Fig. 10)
                  test_index  = prof_clt_residual_Vpeak_index(c,1);
                  
                  % This the peak at the left interpolation edge?
                  if(test_index == lt_clt_index)
                     prof_clt_residual_sacr_flag(c,1)    = -1;
                  end % end if test left index
                  
                  % This the peak at the right interpolation edge?
                  if(test_index == rt_clt_index)
                     prof_clt_residual_sacr_flag(c,1)    =  1;
                  end % end if test right index
                  
               end % end if(prof_clt_residual_num_valid_pts(c,1) >= valid_pts_thres)
               
            end % end if(sum(f) >= valid_pts_thres)
            
            %% Do logic to determine if declutter spectrum should be kept
            
            % Test #1 - (spc(DC) > neighbors) & (spc(DC) > nos_thres)?
            %   must be true because we are modifying spc.
            
            % Test #2 - Does residual have enough points for valid moment?
            %   If there are not enough points, NaN out spc_dcl_lin
            if(prof_clt_residual_num_valid_pts(c,1) < valid_pts_thres)
               spc_dcl_lin                   = spc_dcl_lin .* NaN;
            end % end if(prof_clt_residual_num_valid_pts(c,1) < valid_pts_thres)
            
            % Test #3 - Is residual Vpeak away from interpolation edges?
            %   If peak is at interpolation edges, NaN out spc_dcl_lin
            if(~isnan(prof_clt_residual_sacr_flag(c,1)))
               spc_dcl_lin                   = spc_dcl_lin .* NaN;
            end % end if(~isnan(prof_clt_residual_sacr_flag(c,1))
            
            % Test #4 - Is residual snr < cnr?
            %   If peak is at interpolation edges, NaN out spc_dcl_lin
            %if(prof_clt_residual_snr_dB(c,1) < prof_clt_cnr_dB(c,1))
            %   spc_dcl_lin                   = spc_dcl_lin .* NaN;
            %end % end if(~isnan(prof_clt_residual_sacr_flag(c,1))
            
         end % end if(drop_1bin > drop_dB_thres)
         
      end % end if(spc_dB(DC_index) > nos_threshold_dB) & greater than neighbors
      
      % Save this spectra to the output array   
      prof_spc_dcl_lin(c,:)   = spc_dcl_lin;
      
   end % end if(~isnan(nos_threshold))
   
end % end for c loop
            
            