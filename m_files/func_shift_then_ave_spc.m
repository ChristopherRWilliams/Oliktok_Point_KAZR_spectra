function [spc_desired_output, cnt_spc_output] = ...
   func_shift_then_ave_spc(spc_input, Vmean_orig_input,...
   Vmean_desired_input, Vd_orig, time_orig_minsec, time_desired_minsec)


% For each dwell interval, this routine shifts the spectra to a common
% Vmean, averages spectra, then shi over a defined time window

% updated: 29-July-2017
% =======================================================================

%% Inputs and Outputs

% Inputs
% spc_input                - input spectra, (m,p) dimension: (profile, npts)
% Vmean_input              - mean velocity of each spectrum
% time_orig_minsec_input   - original time of each profile
% time_desired_minsec      - desired start time of each profile

% Outputs
% spc_desired_output       - desired spectra
% cnt_spc_output           - number of spectra used in average

%% Define the output variables

m        = length(time_desired_minsec);
[~,npts]    = size(spc_input);

spc_desired_output   = ones(m,npts) .* NaN;
cnt_spc_output       = ones(m,1) .* NaN;

%% Determine the index for DC

% find the DC location
%DC_index = npts/2 + 1;
% this is the same as
%[~, DC_index] = min(abs(Vd_orig));

dVd   = Vd_orig(6) - Vd_orig(5);

%% Define the start and end times of each dwell

start_time  = time_desired_minsec;
delta_t     = time_desired_minsec(6) - time_desired_minsec(5);

end_time(1:m-1,1)   = start_time(2:m);
end_time(m,1)       = end_time(m-1) + delta_t;

%% Process each time interval

for prof_num = 1:m
   
   % Only do the shift if there is a valid Vmean_desired_input
   if(~isnan(Vmean_desired_input(prof_num)))
      
      % find all profiles within the start & end times
      f = (time_orig_minsec >= start_time(prof_num)) & ...
         (time_orig_minsec <  end_time(prof_num));
      
      if(sum(f) > 0)
         
         % Get the spectra into a common matrix
         spc_temp      = spc_input(f,:);
         Vmean_temp    = Vmean_orig_input(f,:);
         
         % process only valid spectra, which is also when Vmean_temp is valid
         g = ~isnan(Vmean_temp);
         
         % keep the number of valid spectra
         cnt_spc_output(prof_num) = sum(g);
         
         if(cnt_spc_output(prof_num) > 0)
            
            % keep only the valid spectra
            spc_temp      = spc_temp(g,:);
            Vmean_temp    = Vmean_temp(g);
            
            % Shift each spectrum so that Vmean = Vmean_desired
            shift_spc  = spc_temp;
            for i = 1:sum(g)
               
               % Determine how many dVd to shift spectrum
               Vmean_shift_amount = round((Vmean_desired_input(prof_num) - Vmean_temp(i)) / dVd);
               
               % shift the spectrum
               shift_spc(i,:)   = circshift(spc_temp(i,:), Vmean_shift_amount);
               
            end % end for i loop
            
            % calculate the average shifted spectrum
            spc_desired_output(prof_num,:) = mean(shift_spc);
            
         end % end if(cnt_spc_output(prof_num) > 1)
         
      end % end if(sum(f) > 0)
      
   end % end if(~isnan(Vmean_desired_input(prof_num)))
   
end % end for prof_num loop
