function [output_time, output_data_1, output_data_2, output_data_3] = ...
   func_fill_time_gaps_with_NaN_profiles(input_time, input_data_1, input_data_2, input_data_3)
   
% This routine places a NaN profile when the time between profiles is too
% large.

% if there are time-gaps in profiles, the plot will have
% streaks. Put a profile of NaNs between any profiles that are
% too far apart

% Inputs
% input_time      - time vector for plotting
% input_data_1    - 2D matrix of data to be plotted
% input_data_2    - 2D matrix of data to be plotted
% input_data_3    - 2D matrix of data to be plotted

% Outputs
% output_time      - updated time vector for plotting
% output_data_1    - updated 2D matrix of data to be plotted
% output_data_2    - updated 2D matrix of data to be plotted
% output_data_3    - updated 2D matrix of data to be plotted

% updated: 17-May-2017
% ========================================================================

%% Define the output variables

% define the output variables
output_time    = input_time;
output_data_1  = input_data_1;
output_data_2  = input_data_2;
output_data_3  = input_data_3;

%% Determine if any NaN profiles need to be added

% Calculate the time between profiles (no need to estimate last profile)

[n]               = length(input_time);
dif_time          = input_time(2:n) - input_time(1:n-1);
median_dif_time   = median(dif_time);
time_thres        = 3*median(dif_time);

% are there any time steps greater than 3*expected time gap?
f = dif_time > time_thres;
if(sum(f) > 0)
   put_NaN_profile_flag = 1;
else
   put_NaN_profile_flag = 0;
end % end if(sum(f) > 0)

%% Add extra profiles if needed

if(put_NaN_profile_flag)
   
   % f is pointing to the profiles with gaps > time_thres
   % get the index of these gaps
   gap_profile_index = find(f);
   
   % Start from the end and go backwards -
   % put a profile of NaN into the datasets
   [num_NaN_profiles]   = sum(f);
   
   for i = num_NaN_profiles:-1:1
      
      % Get the index - this is the last good profile before
      % the gap
      gap_index   = gap_profile_index(i);
      
      % move the data from the gap to the end over one location
      % work only with the output variables
      output_time(gap_index+2:end+1)       = output_time(gap_index+1:end);
      output_data_1(gap_index+2:end+1,:)   = output_data_1(gap_index+1:end,:);
      output_data_2(gap_index+2:end+1,:)   = output_data_2(gap_index+1:end,:);
      output_data_3(gap_index+2:end+1,:)   = output_data_3(gap_index+1:end,:);
      
      % put a NaN in the i+1 location
      output_time(gap_index+1)           = output_time(gap_index) + median_dif_time;
      output_data_1(gap_index+1,:)       = output_data_1(gap_index,:) .* NaN;
      output_data_2(gap_index+1,:)       = output_data_3(gap_index,:) .* NaN;
      output_data_3(gap_index+1,:)       = output_data_3(gap_index,:) .* NaN;
      
   end % end for i loop
   
end % end if(put_NaN_profile_flag)
