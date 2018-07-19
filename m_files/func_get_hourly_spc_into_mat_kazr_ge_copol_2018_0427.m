function [kaz_ge_orig_spc_lin, kaz_ge_orig_time,kaz_ge_Vd,...
   kaz_ge_range, kaz_ge_radar_const_dB, kaz_ge_Z_near_field_correction,...
   kaz_ge_Nspc] = func_get_hourly_spc_into_mat_kazr_ge_copol_2018_0427(...
   start_year, month_loop, dayofmonth_loop, hour_loop, ...
   raw_netCDF_spc_local_directory, suffix)

% updated = '27-April-2018'

% This routine constructs hourly files of KAZR spectra.
% **********************************************************************


%% Pre-define the output matrices

kaz_ge_orig_spc_lin              = ones(1,1,1) .* NaN;
kaz_ge_orig_time                 = ones(1,7) .* NaN;
kaz_ge_Vd                        = NaN;
kaz_ge_range                     = NaN;
kaz_ge_radar_const_dB            = NaN;
kaz_ge_Z_near_field_correction   = NaN;
kaz_ge_Nspc                      = NaN;

%% Display a message stating which hour of data is being processed

year_str       = num2str(start_year);
month_num_str  = num2str_2digits(month_loop);
dayofmonth_str = num2str_2digits(dayofmonth_loop);
hour_str       = num2str_2digits(hour_loop);

disp(['...starting to process day: ',year_str,'-',month_num_str,'-',dayofmonth_str,', hour: ',hour_str,' ...']);

%% Process one hour of data

start_time_test   = hour_loop;
end_time_test     = hour_loop + 1;

% Get the names of all spc files in the local directory
wildcard_name           = [raw_netCDF_spc_local_directory,suffix];
files_to_process        = dir(wildcard_name);

% the filenames to process are accessed using:
% files_to_process(i).name

num_files_to_process    = length(files_to_process);

if(num_files_to_process > 0)
   
   % reset prof_cnt
   prof_cnt    = 0;
   
   % process each file
   for file_loop = 1:num_files_to_process
      
      % get the next filename to process
      filename_no_dir      = files_to_process(file_loop).name;
      input_spc_filename   = [raw_netCDF_spc_local_directory,filename_no_dir];
      
      %disp(['Processing file ',num2str(file_loop),' out of ',num2str(num_files_to_process),': ',filename_no_dir,' ...']);
      
      % Get the time in this file
      [kaz_ge_time]    = func_read_kazr_spc_time_in_netCDF_file(input_spc_filename);
      
      kaz_ge_time_hr   = kaz_ge_time(:,4) + kaz_ge_time(:,5)./60 + ...
         kaz_ge_time(:,6)./(60*60) + (kaz_ge_time(:,7)./1000)./(60*60);
      
      % Find how many profiles are within the time window
      f = (kaz_ge_time_hr >= start_time_test) & (kaz_ge_time_hr < end_time_test);
      
      % find the files within this day
      g  = (kaz_ge_time(:,3) == dayofmonth_loop);
      
      if(sum(f&g) > 0)
         %disp('paused...')
         %pause
         prof_cnt   = prof_cnt + sum(f);
         
         % get file size information
         % Get all of the variables (except spectra) in this file
         func_read_kazr_netCDF_file_no_spectra
         
         % Define the spectra attributes
         max_nhts    = length(kaz_range);
         max_npts    = length(velocity_bins);
         
      end % end if(sum(f) > 0)
      
   end % end for file_loop
   
   if(prof_cnt > 0)
      
      %% Define the number of profiles to be saved...
      
      kaz_ge_orig_spc_lin   = ones(prof_cnt,max_nhts,max_npts) .* NaN;
      kaz_ge_orig_time      = ones(prof_cnt,7) .* NaN;
      %kazr_spc_orig_Vd        = ones(1,max_npts);
      % kaz_ge_orig_range     = ones(1,max_nhts);
      
      %% Read each file again, and extract the spectra of valid profiles
      
      % reset batch_prof_num
      batch_prof_num    = 0;
      
      % process each file
      for file_loop = 1:num_files_to_process
         
         % get the next filename to process
         filename_no_dir      = files_to_process(file_loop).name;
         input_spc_filename   = [raw_netCDF_spc_local_directory,filename_no_dir];
         
         %disp(' ')
         %disp(['Processing file ',num2str(file_loop),' out of ',num2str(num_files_to_process),'... ',input_spc_filename,' ...']);
         %disp(['Getting spectra, file ',num2str(file_loop),' out of ',num2str(num_files_to_process),': ',filename_no_dir,' ...']);
         
         % Get the time in this file
         [kaz_ge_time]    = func_read_kazr_spc_time_in_netCDF_file(input_spc_filename);
         
         kaz_ge_time_hr   = kaz_ge_time(:,4) + kaz_ge_time(:,5)./60 + ...
            kaz_ge_time(:,6)./(60*60) + (kaz_ge_time(:,7)./1000)./(60*60);
         
         % Find how many profiles are within the time window
         f = (kaz_ge_time_hr >= start_time_test) & (kaz_ge_time_hr < end_time_test);
         
         % find the files within this day
         g  = (kaz_ge_time(:,3) == dayofmonth_loop);
         
         if(sum(f&g) > 0)
            
            %disp('paused...')
            %pause
            
            % Get all of the variables (except spectra) in this file
            func_read_kazr_netCDF_file_no_spectra
            
            %% Get the Spectra Information from the netCDF file
            
            % open the netCDF file, again (it was closed after getting the
            % moments)
            
            % Open the netCDF file
            ncid = netcdf.open(input_spc_filename,'NC_NOWRITE');
            
            % Get the number of variables in this data file
            [~, nvars, ~, ~] = netcdf.inq(ncid);
            
            % Find the location of the variable named 'spectra'
            variable_num = -1;
            for rr = 0:nvars-1
               
               %[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,rr);
               [varname, ~, ~, ~] = netcdf.inqVar(ncid,rr);
               
               % Is this variable named 'Spectra'?
               %variable_length = length(varname);
               if(length(varname) == 7)
                  compare_string = varname == 'spectra';
                  if(sum(compare_string) == 7)
                     variable_num = rr;
                  end % end if(sum(compare_string) == 7)
               end % end if(length(varname) == 7)
               
            end % end for rr loop
            
            % The variable_num points to the variable named 'Spectra'
            
            %[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,variable_num);
            
            %[dimname, dimlen] = netcdf.inqDim(ncid,dimids(1));
            
            % Get the attributes
            varid       = netcdf.inqVarID(ncid,'spectra');
            
            spectra_scale_factor    = netcdf.getAtt(ncid,varid,'scale_factor');
            spectra_add_offset      = netcdf.getAtt(ncid,varid,'add_offset');
            %spectra_comment         = netcdf.getAtt(ncid,varid,'comment');
            
            % convert the single to double precision
            spectra_scale_factor    = double(spectra_scale_factor);
            spectra_add_offset      = double(spectra_add_offset);
            
            % Get the radar calibration constant, stored as a global
            % attribute
            n                       = length(cal_constant);
            kaz_ge_radar_const_dB   = str2double(cal_constant(1:n-4));
            
            % Get the number of spectral averages
            Nspc                 = str2double(num_spectral_averages);
            
            % The dimensions of the variable Spectra are: [n 256], where n is
            % the number of spectra saved in the compressed data file. Spectra
            % without large enough SNR are not saved.
            
            % The locator_mask(r,c) has the index number for the spectra that
            % corresponds to the pixel (r,c).
            
            % Display a message
            %                   disp(['spectra scale factor: ',num2str(spectra_scale_factor)]);
            %                   disp(['spectra offset      : ',num2str(spectra_add_offset)]);
            %                   disp(['spectra radar const : ',num2str(kaz_radar_const_dB)]);
            %                   disp(['Num spectra averaged: ',num2str(Nspc)]);
            %
            %% Calculate the near-field reflectivity correction
            
            % The far field is approximately
            % r_far = 2*D^2/lambda;
            
            % Get the radar operatering freq, stored as a global attribute
            n = length(radar_operating_frequency);
            radar_operating_frequency_Hz = str2double(radar_operating_frequency(1:n-4));
            
            % estimate the radar operating wavelength
            lambda      = (3*10^8)/(radar_operating_frequency_Hz);
            
            % Get the antenna diameter, stored as a global attribute
            if(exist('antenna_diameter','var') == 1)
               n = length(antenna_diameter);
               antenna_diameter_m = str2double(antenna_diameter(1:n-1));
            else
               antenna_diameter_m = 2;
            end % end (exist('antenna_diameter','var') == 1)
            
            % calculate the far field:
            r_far       = 2*antenna_diameter_m^2/lambda;
            %disp(['The operating frequency is: ',num2str(radar_operating_frequency_Hz./(10^9)),' GHz']);
            %disp(['The far field is at range: ',num2str(r_far),' m']);
            
            % The universal near-field reflectivity correction is:
            % numerator = (5.26*10^-5) *(r0/r_far)^2.50;
            % denominator = 0.0117 + (r0/r_far).^2.50;
            % F(r0)/F0 = numerator / denominator;
            
            num_gates               = length(kaz_range);
            r0_r_far_ratio          = ones(1,num_gates) .* NaN;
            Z_near_field_correction = ones(1,num_gates) .* NaN;
            
            for c = 1:num_gates
               r0_r_far_ratio(c) = kaz_range(c) / r_far;
               if(kaz_range(c) < r_far)
                  Z_near_field_correction(c) = 10.*log10(((5.26*10^-5) + (r0_r_far_ratio(c))^2.50) / (0.0117 + (r0_r_far_ratio(c)).^2.50));
               else
                  Z_near_field_correction(c) = 0;
               end % end if(KAZR_range < r_far)
            end % end for c loop
            
            % rename variables
            kaz_ge_Z_near_field_correction  = Z_near_field_correction;
            
            
            % get some constants for extracting the spectra
            npts         = length(velocity_bins);
            count_value  = [npts 1];
            
            %% Process each valid profile
            
            % process the profiles that have good times
            % get the indices of good profiles
            prof_index  = find(f);
            
            % for each valid profile
            for r = 1:length(prof_index)
               
               % increment the batch_prof_num
               batch_prof_num    = batch_prof_num + 1;
               
               % get the profile number for this opened file
               prof_num = prof_index(r);
               
               if(rem(prof_num,200) == 0)
                  disp(['.....getting spectra, profile #: ',num2str(prof_num),'...',...
                     num2str_2digits(kaz_ge_time(prof_num,4)),':',...
                     num2str_2digits(kaz_ge_time(prof_num,5)),':',...
                     num2str_2digits(kaz_ge_time(prof_num,6))]);
               end % end if((rem(r,10) == 0)
               
               % Get the time for this profile
               kaz_ge_orig_time(batch_prof_num,:)   = kaz_ge_time(prof_num,:);
               
               % Get a profile of spectra and expand
               for c = 1:max_nhts
                  
                  if(locator_mask(prof_num,c) >= 0)
                     %range_KAZR(c) = range(c);
                     
                     % Define the starting point:
                     start_value = [0 locator_mask(prof_num,c)];
                     
                     %data = netcdf.getVar(ncid,variable_num,start_value,count_value);
                     data           = netcdf.getVar(ncid,variable_num,start_value,count_value,'int16');
                     
                     % scale and offset the spectra
                     spc_orig_dB       = (double(data) .* spectra_scale_factor) + spectra_add_offset;
                     spc_orig_lin      = 10.^(spc_orig_dB./10);
                     
                     % Flip the spectra, so + is downward (Doppler)
                     % And expand to 3*Vd_orig
                     % ***Process original spectra***
                     
                     % keep value in location #1, -Nyquist velocity
                     % flip all other points.
                     flip_spc_orig_lin      = spc_orig_lin;
                     for k = 2:length(spc_orig_lin)
                        flip_spc_orig_lin(k) = spc_orig_lin(length(spc_orig_lin)+2 - k);
                     end % end for c loop
                     
                     %% Expand the spectra to correct for aliasing
                     %spc_expand     = [flip_spc_orig_lin ; flip_spc_orig_lin; flip_spc_orig_lin];
                     
                     % save in the array of profiles
                     kaz_ge_orig_spc_lin(batch_prof_num,c,:)    = flip_spc_orig_lin;
                     %spc_kazr_orig_lin(prof_num,c,:)    = spc_expand;
                     
                     %disp('paused, line 394...')
                     %pause
                     
                  end % end if(locator_mask(prof_num,c) >= 0)
               end % end for c loop
               
            end % end for r loop
            
            %% close the opened file
            
            netcdf.close(ncid)
            
         end % end if(sum(f) > 0)
         
      end % end for file_loop
      
      %% save the spectra in a mat file
      
      % rename some variables
      kaz_ge_Vd       = velocity_bins;
      kaz_ge_range    = kaz_range;
      kaz_ge_Nspc     = Nspc;
      
      % define the output variable name
      %spc_directory        = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\mat_hourly_orig_spc_files\';
      %output_mat_spc_filename  = [spc_directory,'oli_orig_spc_kazr_ge_copol_',year_str,month_num_str,dayofmonth_str,'_',hour_str,'.mat'];
      
      %save(output_mat_spc_filename,'kaz_ge_orig_spc_lin',...
      %   'kaz_ge_orig_time','kaz_ge_Vd',...
      %   'kaz_ge_range','kaz_ge_radar_const_dB',...
      %   'kaz_ge_Z_near_field_correction','kaz_ge_Nspc','-v7.3');
      
   end % end if(prof_cnt > 0)
   
end % end if(num_files_to_process > 0)
