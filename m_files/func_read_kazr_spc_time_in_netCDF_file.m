function [kaz_spc_time] = func_read_kazr_spc_time_in_netCDF_file(input_spc_filename)

% Program: func_read_kazr_spc_time_in_netCDF_file

% This routine only extracts the time stamp from the netCDF file

% updated: 17-June-2017
% ========================================================================

%% Open the netCDF file
ncid = netcdf.open(input_spc_filename,'NC_NOWRITE');

%% Get the number of variables in this data file

%[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
[~, nvars, ~, ~] = netcdf.inq(ncid);

%% Load base_time and time_offset variables into MATLAB

for r = 0:nvars-1
    
    %[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,r);
    [varname, ~, ~, ~] = netcdf.inqVar(ncid,r);

    % Extract variable: base_time
    if(length(varname) == 9)
       compare_string = varname == 'base_time';
       if(sum(compare_string) == 9)
          
          % get this variable
          base_time = netcdf.getVar(ncid,r);
       end % end if(sum(compare_string) == 7)
    end % end if(length(varname) == 7)

    % Extract variable: time_offset
    if(length(varname) == 11)
       compare_string = varname == 'time_offset';
       if(sum(compare_string) == 11)
          
          % get this variable
          time_offset = netcdf.getVar(ncid,r);
       end % end if(sum(compare_string) == 7)
    end % end if(length(varname) == 7)
        
end % end for r loop

%% Convert base_time from Epoch time (seconds since Jan 1, 1970)

[kaz_spc_time, ~] = func_convert_netCDF_time_to_vector(base_time, time_offset);

%% close the opened file

netcdf.close(ncid)

