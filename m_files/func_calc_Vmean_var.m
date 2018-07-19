function    [Vmean_var_output, cnt_Vmean_var_output] = ...
      func_calc_Vmean_var(Vmean_input, ...
      time_orig_minsec, time_desired_minsec)


% For each dwell interval, this routine shifts the spectra to a common
% Vmean, averages spectra, then shi over a defined time window

% updated: 29-July-2017
% =======================================================================

%% Inputs and Outputs

% Inputs
% Vmean_input              - mean velocity of each spectrum
% cnt_threshold            - number of needed samples
% time_orig_minsec_input   - original time of each profile
% time_desired_minsec      - desired start time of each profile

% Outputs
% spc_desired_output       - desired spectra
% cnt_spc_output           - number of spectra used in average

%% Define the output variables

m                       = length(time_desired_minsec);
Vmean_var_output        = ones(m,1) .* NaN;
cnt_Vmean_var_output    = ones(m,1) .* NaN;

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
      Vmean_temp    = Vmean_input(f,:);
      
      % process only valid spectra, which is also when Vmean_temp is valid
      g = ~isnan(Vmean_temp);
      
      % keep the number of valid spectra
      cnt_Vmean_var_output(prof_num) = sum(g);
      
      if(cnt_Vmean_var_output(prof_num) > 2)
         
         % Calculate the Vmean_var
         Vmean_var_output(prof_num)    = std(Vmean_temp(g)).^2;
         
      end % end if(cnt_Vmean_var_output(prof_num) > 2)
      
   end % end if(sum(f) > 0)
   
   
end % end for prof_num loop
