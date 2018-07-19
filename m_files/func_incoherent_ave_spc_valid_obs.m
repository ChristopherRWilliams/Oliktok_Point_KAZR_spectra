function [spc_desired_output, cnt_spc_output]   = func_incoherent_ave_spc_valid_obs(spc_input, valid_obs_input, time_orig_minsec, time_desired_minsec)

% This routine averages spectra over a defined time window

% updated: 5-July-2017
% =======================================================================

%% Inputs

% spc_input    - input spectra, (m,p) dimension: (profile, npts)
% time_orig_minsec_input   - original time of each profile
% time_desired_minsec   - desired start time of each profile

%% Define the output variables

m        = length(time_desired_minsec);
[~,p]    = size(spc_input);

spc_desired_output   = ones(m,p) .* NaN;
cnt_spc_output       = ones(m,1) .* NaN;

%% Define the start and end times of each dwell

start_time  = time_desired_minsec;
delta_t     = time_desired_minsec(6) - time_desired_minsec(5);

end_time(1:m-1,1)   = start_time(2:m);
end_time(m,1)       = end_time(m-1) + delta_t;

%% Process each time interval

for prof_num = 1:m
   
   % find all profiles within the start & end times   
   f = (time_orig_minsec >= start_time(prof_num)) & ...
       (time_orig_minsec <  end_time(prof_num));

    if(sum(f) > 0)
       
       % Get the spectra into a common matrix
       spc_temp         = spc_input(f,:);
       
       valid_obs_temp   = valid_obs_input(f);
       
       % how many of these spectra are valid?
       % look at first element of each spectrum
       spc_first_bin    = spc_temp(:,1);
       g = ~isnan(spc_first_bin);

       % Also use the valid_obs flag, 1 = good moment, NaN = clutter
       k = ~isnan(valid_obs_temp);
       
       % keep the number of valid spectra
       cnt_spc_output(prof_num) = sum(g & k);

       if(cnt_spc_output(prof_num) > 1)
          
         % average the spectra within this interval
         spc_desired_output(prof_num,:) = mean(spc_temp(g & k,:));
       
       end % end if(cnt_spc_output(prof_num) > 1)
       
    end % end if(sum(f) > 0)

end % end for prof_num loop

   