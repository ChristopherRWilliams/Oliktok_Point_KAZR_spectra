function [snr_out, nos_out, zdb_out, Vmean_out, Vsig_out, ...
   Vskew_out, Vkurt_out, Vpeak_out, Ppeak_out, ...
   Vslope_lt_out, Vslope_rt_out, std_1bin_out, std_2bin_out, ...
   max_1bin_delta_out, max_2bin_delta_out] = ...
   func_find_moments(spk_short, Vd_short, nos_mean, Npts, ...
   ht_ref_dB, radar_const_dB)

% inputs:
% spk_short       - spectrum (linear units)
% Vd_short        - Doppler velocity bins
% Npts            - length of original spectrum
% ht_ref_dB       - height expressed in dB, 20*log10(ht)
% radar_const_dB  - calibration constant

% outputs

%% Estimate the moments using short arrays

if(sum(spk_short-nos_mean) > 0)

   % noise power
   noise_power = nos_mean * Npts;
   nos_out = 10*log10(noise_power);
   
   % snr
   snr_out = 10*log10(sum(spk_short - nos_mean) / noise_power);
      
   % zdb (dB) = 10*log10(sum(spk_short - nos_mean) + 20*log10(ht_ref) + radar_const_dB;
   % zdb (dB) = snr_out + nos_out +  20*log10(ht_ref) + radar_const_dB;
   zdb_out = (snr_out + nos_out) + ht_ref_dB + radar_const_dB;
   
   % Vmean
   Vmean_out = sum(spk_short .* Vd_short) / sum(spk_short);
   
   % Vsig
   Vsig_out = sqrt(sum((Vd_short - Vmean_out).^2 .* ...
      spk_short) / sum(spk_short));
   
   % Calculate the skewness
   % skewness is: sum((spk_short - Vmean).^3 .* S(v)) /sum(S(v))
   %              ----------------------------------------------
   %                 Vsig.^3
   
   skewness_top      = (sum((Vd_short - Vmean_out).^3 .* spk_short))/ sum(spk_short);
   skewness_bot      = Vsig_out.^3;
   Vskew_out         = skewness_top / skewness_bot;
   
   % Calculate the kurtosis
   % kurtosis is: sum((spk_short - Vmean).^4 .* S(v)) /sum(S(v))
   %              ----------------------------------------------
   %                 Vvar.^2
   
   kurtosis_top      = (sum((Vd_short - Vmean_out).^4 .* spk_short))/ sum(spk_short);
   kurtosis_bot      = Vsig_out.^4;
   Vkurt_out         = kurtosis_top / kurtosis_bot;
   
   % Find the velocity peak magnitude (Vpeak)
   [~, imax] = max(spk_short);
   
   Vpeak_out   = Vd_short(imax);
   Ppeak_out   = 10.*log10(spk_short(imax));
   nos_mean_dB = 10.*log10(nos_mean);
   
   % Calculate the left and right slopes
   % ===================================
   
   % look to the left
   Vrange_left    = Vd_short(imax) - Vd_short(1);   
   
   % look to the right
   Vrange_right   = Vd_short(end) - Vd_short(imax);
   
   Vslope_lt_out     = (Ppeak_out - nos_mean_dB) / Vrange_left;
   Vslope_rt_out     = (Ppeak_out - nos_mean_dB) / Vrange_right;
   
   % Calculate the 1-bin neighbor variance
   % =====================================
   spk_short_dB   = 10.*log10(spk_short);
   
   delta_1bin     = spk_short_dB(2:end) - spk_short_dB(1:end-1);
   std_1bin_out   = std(delta_1bin);
   
   delta_2bin     = spk_short_dB(3:end) - spk_short_dB(1:end-2);
   std_2bin_out   = std(delta_2bin);
   
   % calculate the 1bin drop from the peak
   % =====================================
   [~,index_peak]       = max(spk_short_dB);
   
   index_lt             = index_peak - 1;
   index_rt             = index_peak + 1;
   
   % look to the left and right
   if(index_lt >= 1)
      drop_lt           = spk_short_dB(index_peak) - spk_short_dB(index_lt);
   else
      drop_lt           = 0;
   end % end if(lt_index >= 1)

   if(index_rt <= length(spk_short_dB))
      drop_rt            = spk_short_dB(index_peak) - spk_short_dB(index_rt);
   else
      drop_rt           = 0;
   end % end if(rt_index <= length(spk_short_dB))
   
   if(drop_lt > drop_rt)
      max_1bin_delta_out   = drop_lt;
   else
      max_1bin_delta_out   = drop_rt;
   end % end if(drop_lt > drop_rt)
   
   % calculate the 2bin drop from the peak
   % =====================================
   
   index_lt             = index_peak - 2;
   index_rt             = index_peak + 2;
   
   % look to the left and right
   if(index_lt >= 1)
      drop_lt           = spk_short_dB(index_peak) - spk_short_dB(index_lt);
   else
      drop_lt           = 0;
   end % end if(lt_index >= 1)

   if(index_rt <= length(spk_short_dB))
      drop_rt            = spk_short_dB(index_peak) - spk_short_dB(index_rt);
   else
      drop_rt           = 0;
   end % end if(rt_index <= length(spk_short_dB))
   
   if(drop_lt > drop_rt)
      max_2bin_delta_out   = drop_lt;
   else
      max_2bin_delta_out   = drop_rt;
   end % end if(drop_lt > drop_rt)

end % end if(sum(spk_short - mean_noise) > 0)
