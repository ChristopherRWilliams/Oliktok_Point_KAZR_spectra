%% Program: func_read_kazr_netCDF_file_no_spectra

% updated: 28-December-2016

% Open the netCDF file
ncid = netcdf.open(input_spc_filename,'NC_NOWRITE');

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
    get_data = 1;
    
    % Check to see if this variable is 'spectra'
    % if yes, then set get_data to 0
    if(length(varname) == 7)
       compare_string = varname == 'spectra';
       if(sum(compare_string) == 7)
          get_data = 0;
       end % end if(sum(compare_string) == 7)
       %compare_string = varname == 'Spectra';
       %if(sum(compare_string) == 7)
       %   get_data = 0;
       %end % end if(sum(compare_string) == 7)
    end % end if(length(varname) == 7)

    % Rename the 'range' variable on reading
    
    if(length(varname) == 5)
       compare_string = varname == 'range';
       if(sum(compare_string) == 5)

         data = netcdf.getVar(ncid,r);
       
         % convert the data into a double
         kaz_range = double(data);
       
         % Rotate the matrices so that time is in rows and range is in columns.
         n_dim = length(size(data));  % check to see that it is a 2-dimension variable
       
         % switch the 2-dim matrices
         if((n_dim > 1) && (n_dim < 3))
             % Do the permutation if the second dimension is not equal to 1
            [m, n] = size(kaz_range);
            if(n > 1)
               kaz_range = permute(kaz_range,[2 1]);
            end % end if(n > 1)
         end % end if((n_dim > 1) && (n_dim < 3))
       end % end if(sum(compare_string) == 5)
    end % end if(length(varname) == 5)
       
    % Rename the 'range_bounds' variable on reading
    if(length(varname) == 12)
       compare_string = varname == 'range_bounds';
       if(sum(compare_string) == 12)

         data = netcdf.getVar(ncid,r);
       
         % convert the data into a double
         kaz_range_bounds = double(data);
       
         % Rotate the matrices so that time is in rows and range is in columns.
         n_dim = length(size(data));  % check to see that it is a 2-dimension variable
       
         % switch the 2-dim matrices
         if((n_dim > 1) && (n_dim < 3))
             % Do the permutation if the second dimension is not equal to 1
            [m, n] = size(kaz_range_bounds);
            if(n > 1)
               kaz_range_bounds = permute(kaz_range_bounds,[2 1]);
            end % end if(n > 1)
         end % end if((n_dim > 1) && (n_dim < 3))
       end % end if(sum(compare_string) == 12)
    end % end if(length(varname) == 12)    
    
    %% process each of the other variables
    
    if(get_data)
       data = netcdf.getVar(ncid,r);
       
       % convert the data into a double
       data = double(data);
       
       % Rotate the matrices so that time is in rows and range is in columns.
       n_dim = length(size(data));  % check to see that it is a 2-dimension variable
       
       % switch the 2-dim matrices
       if((n_dim > 1) && (n_dim < 3))
          % Do the permutation if the second dimension is not equal to 1
          [m, n] = size(data);
          if(n > 1)
             data = permute(data,[2 1]);
          end % end if(n > 1)
       end % end if((n_dim > 1) && (n_dim < 3))
       
       % switch the 3-dim matrices
       %if(n > 2)
       %   data = permute(data,[3 2 1]);
       %end % end if(n > 12)
       
       % Put the data into the variable name
       eval([varname '= data;']);
       
    end % end if(get_data)
    % Get the attributes
    %for c = 0:numatts-1
    %    attname = netcdf.inqAttName(ncid,r,c);
    %    % Get value of attribute.
    %    attval = netcdf.getAtt(ncid,r,attname);
    %
    %    disp(['Attribute: ',attname,...
    %        ', ',attval]);
    %    
    %end % end for c loop

    %disp(' ');    
    
end % end for r loop

% Get the global attributes
for c = 0:ngatts-1
   
   % Get the attribute name
   global_att_name = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),c);
   
   % Get the attribute value
   gattval           = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),global_att_name);
   
   % Save the value in the new variable named by global_att_name
   eval([global_att_name '= gattval;']);
   
end % end for c loop

% close the opened file
netcdf.close(ncid)

