%% Program: func_read_nc_file

% updated: 05-Oct-2017

% This routine reads the omega velocity netCDF data file.

% Open the netCDF file
ncid = netcdf.open(netCDF_input_filename,'NC_NOWRITE');

% Get the number of variables in this data file
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get the information about each dimension
%for r = 0:ndims-1
%    [dimName dimLen] = netcdf.inqDim(ncid,r);
%
%    disp(['Dimension Name: ',dimName,', Length: ',num2str(dimLen)]);
%end % end for r loop

% Load each Variable into MATLAB
for r = 0:nvars-1
   
   [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,r);
   
   %disp(['Variable Name: ',varname,', Dim IDs: ',num2str(dimids),...
   %    ', Number of Attributes: ',num2str(numatts)]);
   
   %disp(['Processing variable: ',varname,', ...']);
   
   % Get the data
   
   % Test to make sure this is the variable named 'Spectra'
   
   % set the flag to true
   data = netcdf.getVar(ncid,r);
   
   % convert the data into a double
   data = double(data);
   
   
   % Put the data into the variable name
   eval([varname '= data;']);
   
   
end % end for r loop

% close the opened file
netcdf.close(ncid)

