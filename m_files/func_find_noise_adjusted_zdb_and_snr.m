function [zdb_adj, snr_adj, nos_adj, new_nos_power] = func_find_noise_adjusted_zdb_and_snr(snr_input, nos_input, range_input, radar_const_dB)

% This routine finds the noise adjusted snr and reflectivity with PRC = 1.
% Estimate the reflectivity for a gain of PRC = 1;
% Z = 10*log10(PRC * r^2 * SNR_linear)
% Z = 10*log10(PRC) + 20*log10(r) + 10*log10(SNR_linear)
% Z = 10*log10(PRC) + 20*log10(r) + SNR_dB
%    if PRC = 1, then:
% Z = 0 + 20*log10(r) + SNR_dB
% The radar constant is given in dB, so, reflectivity is:
% Z = radar_const_dB + 20*log10(r) + SNR_dB;

% The SNR_dB needs to be adjusted by the mean noise observed at other hts:
% SNR_dB_adj = SNR_dB_old + noise_old(in dB) - noise_new(in dB);

% updated:05-August-2016

%% define the output variables

n                   = length(snr_input);
zdb_adj             = ones(1,n) .* NaN;
snr_adj             = ones(1,n) .* NaN;
nos_adj             = ones(1,n) .* NaN;

%% Define the range in log

range_20log10r      = 20*log10(range_input);

%% Find the new noise power

% Don't use the noise estimates below 1000 m
gg = ~isnan(nos_input) & (range_input > 1000);

% do the noise adjust if there are at least 10 range gates
if(sum(gg) >= 10)
   
   % find the average of the 10 smallest noise powers
   nos_sort = sort(nos_input(gg));
   
   new_nos_power = mean(nos_sort(1:10));
   
else
   new_nos_power = min(nos_input);
end % end if(sum(gg) >= 10)

%% Find the noise adjusted snr and Z (with PRC = 1)

% process each range
for j = 1:n
      
   % Calculate the nos adjustment
   nos_adj(j)  = nos_input(j) - new_nos_power;
   
   % adjust the SNR_dB
   % snr_adj(j) = snr_input(j) + (nos_input(j) - new_nos_power);
   snr_adj(j)  = snr_input(j) + nos_adj(j);

   % find the reflectivity with PRC = 1
   %zdb_adj(j)     = radar_const_dB + range_20log10r(j) + snr_adj(j);
   % The DOE radar constant is based on total power, not SNR.
   zdb_adj(j)     = radar_const_dB + range_20log10r(j) + snr_adj(j) + new_nos_power;

end % end for j loop
