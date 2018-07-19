function [profile_time_vec, time_minsec] = func_convert_netCDF_time_to_vector(base_time, time_offset)

% Inputs:
% base_time       - time of profile since January 1, 1970.

% Outputs:
% profile_time_vec   - vector format of profile time:
%                    - 1: year, 2: month, 3: dayofmonth, 4: hour
%                    - 5: minute, 6: second, 7: millisecond

% updated: 27-November-2013
% =======================================================================

%% Convert base_time from Epoch time (seconds since Jan 1, 1970)

% Get the number of days since January 1, 1970
% (base_time needs to be converted to double):
base_time_in_days_since_epoch = double(base_time) / (24*60*60);

% Get the number of days from January 1, 0000 to the Epoch date:
time_in_days_from_0000_to_epoch = datenum([1970 1 1 0 0 0]);

%The base_time expressed using the MATLAB convention is:
base_time_in_days_MATLAB = time_in_days_from_0000_to_epoch + ...
   base_time_in_days_since_epoch;

%To convert from the MATLAB daynumber to vector use:
base_time_vector = datevec(base_time_in_days_MATLAB);
% format is: [year, month, dayofmonth, hour, minute, second]

%Get the time for each profile:

%The time of the first profile is base_time + time_offset(1).

%Since base_time_in_days_MATLAB has the units of days, the
% time_offset needs to be converted into fractions of day:
% time_offset_day = time_offset ./ (24*60*60);

%Since the time_offset is in seconds, the first profile is at
% base_time_in_days_MATLAB + time_offset(1)

%The vector date for profile i is given by:
m = length(time_offset);
profile_time_day = zeros(m,1);
profile_time_vec = zeros(m,7);

for r = 1:m
   
   profile_time_day(r) = base_time_in_days_MATLAB + ...
      time_offset(r)/(24*60*60);
   
   % Get the time vector for each profile
   profile_time_vec(r,1:6) = datevec(profile_time_day(r));
   
   second_fraction         = profile_time_vec(r,6);
   second                  = floor(second_fraction);
   millisecond             = (second_fraction - second) * 1000;
   profile_time_vec(r,6)   = second;
   profile_time_vec(r,7)   = millisecond;
   
end % end for r loop

time_minsec = profile_time_vec(:,5) + profile_time_vec(:,6)/60 + ...
   (profile_time_vec(:,7)/1000)/60;
