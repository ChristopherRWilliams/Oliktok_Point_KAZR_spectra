% Program = 'main_save_daily_ge_15sec_all_mom_as_netCDF_2018_0718'
% updated = '18-July-2018'

% This routine reads KAZR spectra and estimates moments.

% **********************************************************************

%% Define the drive letters

site_prefix             = 'oli';
version_suffex          = '_v15';

%% Define the days to process

start_year  = 2016;
start_month =    6;
end_month   =    6;
start_day   =   19;
end_day     =   19;

for month_loop = start_month:end_month
   
   for dayofmonth_loop = start_day:end_day
      
      
      %% Define the file to be loaded
      
      year_str       = num2str(start_year);
      month_num_str  = num2str_2digits(month_loop);
      dayofmonth_str = num2str_2digits(dayofmonth_loop);

      input_directory      = '..\mat_daily_15sec_ave_moments\';
      %input_directory      = 'L:\oli\kazr\multi_peak_processing_version_2017_1021\ge_copol\mat_daily_15sec_ave_moments\';
      filename_root        = [site_prefix,'_kazr_ge_15sec_shift_mom_daily_'];
      input_filename       = [input_directory,filename_root,year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
      
      good_input_filename  = exist(input_filename,'file');
      
      %% Load in the daily files
      if(good_input_filename == 2)
         
         
         %% load in the kazr moments
         
         disp(['loading mat file: ',input_filename,' ...']);
         
         % the input file is named: input_netCDF_filename
         %func_load_all_netCDF_variables
         load(input_filename)
         % the variables are:
         %   kaz_ge_range                         608x1                     4864  double
         %   kaz_ge_separate_Ppeak               5760x608x4            112066560  double
         %   kaz_ge_separate_Vkurt               5760x608x4            112066560  double
         %   kaz_ge_separate_Vleft               5760x608x4            112066560  double
         %   kaz_ge_separate_Vmean               5760x608x4            112066560  double
         %   kaz_ge_separate_Vpeak               5760x608x4            112066560  double
         %   kaz_ge_separate_Vright              5760x608x4            112066560  double
         %   kaz_ge_separate_Vsig                5760x608x4            112066560  double
         %   kaz_ge_separate_Vskew               5760x608x4            112066560  double
         %   kaz_ge_separate_Vslope_lt           5760x608x4            112066560  double
         %   kaz_ge_separate_Vslope_rt           5760x608x4            112066560  double
         %   kaz_ge_separate_index_lt            5760x608x4            112066560  double
         %   kaz_ge_separate_index_rt            5760x608x4            112066560  double
         %   kaz_ge_separate_max_1bin_delta      5760x608x4            112066560  double
         %   kaz_ge_separate_max_2bin_delta      5760x608x4            112066560  double
         %   kaz_ge_separate_nos                 5760x608x4            112066560  double
         %   kaz_ge_separate_snr                 5760x608x4            112066560  double
         %   kaz_ge_separate_std_1bin            5760x608x4            112066560  double
         %   kaz_ge_separate_std_2bin            5760x608x4            112066560  double
         %   kaz_ge_separate_zdb_arm             5760x608x4            112066560  double
         %   kaz_ge_separate_zdb_crw             5760x608x4            112066560  double
         %   kaz_ge_single_Ppeak                 5760x608               28016640  double
         %   kaz_ge_single_Vkurt                 5760x608               28016640  double
         %   kaz_ge_single_Vleft                 5760x608               28016640  double
         %   kaz_ge_single_Vmean                 5760x608               28016640  double
         %   kaz_ge_single_Vmean_var             5760x608               28016640  double
         %   kaz_ge_single_Vpeak                 5760x608               28016640  double
         %   kaz_ge_single_Vright                5760x608               28016640  double
         %   kaz_ge_single_Vsig                  5760x608               28016640  double
         %   kaz_ge_single_Vskew                 5760x608               28016640  double
         %   kaz_ge_single_Vslope_lt             5760x608               28016640  double
         %   kaz_ge_single_Vslope_rt             5760x608               28016640  double
         %   kaz_ge_single_cnt_Vmean_var         5760x608               28016640  double
         %   kaz_ge_single_index_lt              5760x608               28016640  double
         %   kaz_ge_single_index_rt              5760x608               28016640  double
         %   kaz_ge_single_max_1bin_delta        5760x608               28016640  double
         %   kaz_ge_single_max_2bin_delta        5760x608               28016640  double
         %   kaz_ge_single_nos                   5760x608               28016640  double
         %   kaz_ge_single_nos_adj               5760x608               28016640  double
         %   kaz_ge_single_nos_max               5760x608               28016640  double
         %   kaz_ge_single_nos_mean              5760x608               28016640  double
         %   kaz_ge_single_nos_profile           5760x1                    46080  double
         %   kaz_ge_single_nos_std               5760x608               28016640  double
         %   kaz_ge_single_snr                   5760x608               28016640  double
         %   kaz_ge_single_snr_adj               5760x608               28016640  double
         %   kaz_ge_single_std_1bin              5760x608               28016640  double
         %   kaz_ge_single_std_2bin              5760x608               28016640  double
         %   kaz_ge_single_zdb_arm               5760x608               28016640  double
         %   kaz_ge_single_zdb_arm_near          5760x608               28016640  double
         %   kaz_ge_single_zdb_crw               5760x608               28016640  double
         %   kaz_ge_sub_Ppeak                    5760x608x4            112066560  double
         %   kaz_ge_sub_Vkurt                    5760x608x4            112066560  double
         %   kaz_ge_sub_Vleft                    5760x608x4            112066560  double
         %   kaz_ge_sub_Vmean                    5760x608x4            112066560  double
         %   kaz_ge_sub_Vpeak                    5760x608x4            112066560  double
         %   kaz_ge_sub_Vright                   5760x608x4            112066560  double
         %   kaz_ge_sub_Vsig                     5760x608x4            112066560  double
         %   kaz_ge_sub_Vskew                    5760x608x4            112066560  double
         %   kaz_ge_sub_Vslope_lt                5760x608x4            112066560  double
         %   kaz_ge_sub_Vslope_rt                5760x608x4            112066560  double
         %   kaz_ge_sub_index_lt                 5760x608x4            112066560  double
         %   kaz_ge_sub_index_rt                 5760x608x4            112066560  double
         %   kaz_ge_sub_max_1bin_delta           5760x608x4            112066560  double
         %   kaz_ge_sub_max_2bin_delta           5760x608x4            112066560  double
         %   kaz_ge_sub_nos                      5760x608x4            112066560  double
         %   kaz_ge_sub_snr                      5760x608x4            112066560  double
         %   kaz_ge_sub_std_1bin                 5760x608x4            112066560  double
         %   kaz_ge_sub_std_2bin                 5760x608x4            112066560  double
         %   kaz_ge_sub_zdb_arm                  5760x608x4            112066560  double
         %   kaz_ge_sub_zdb_crw                  5760x608x4            112066560  double
         %   kaz_ge_time_stamp                   5760x7                   322560  double
                             
         %% Rename the variables for output in netCDF format
                  
         kaz_ge_single_year         = kaz_ge_time_stamp(:,1);
         kaz_ge_single_month        = kaz_ge_time_stamp(:,2);
         kaz_ge_single_dayofmonth   = kaz_ge_time_stamp(:,3);
         kaz_ge_single_hour         = kaz_ge_time_stamp(:,4);
         kaz_ge_single_minute       = kaz_ge_time_stamp(:,5);
         kaz_ge_single_second       = kaz_ge_time_stamp(:,6);
         kaz_ge_single_millisecond  = kaz_ge_time_stamp(:,7);
         
         % Get the ARM base_time 
         ARM_time                   = datenum([kaz_ge_single_year kaz_ge_single_month kaz_ge_single_dayofmonth kaz_ge_single_hour kaz_ge_single_minute kaz_ge_single_second]);
         base_time                  = (ARM_time(1) - datenum([1970 1 1 0 0 0]))*86400;
         time_offset                = (ARM_time - ARM_time(1)).*86400;
         ARM_midnight               = datenum([kaz_ge_single_year(1) kaz_ge_single_month(1) kaz_ge_single_dayofmonth(1) 0 0 0]);
         time_from_midnight         = (ARM_time - ARM_midnight).*86400;

         % get the range variables 
         kaz_ge_single_range        = kaz_ge_range';
         dht                        = kaz_ge_range(6) - kaz_ge_range(5);
         kaz_ge_single_range_bounds(1,:) = kaz_ge_range' - dht/2;
         kaz_ge_single_range_bounds(2,:) = kaz_ge_range' + dht/2;
                           
         kaz_ge_single_HS_nos_thres = kaz_ge_single_nos_max;
         kaz_ge_single_zdb          = kaz_ge_single_zdb_crw;

         % rename the reflectivity
         kaz_ge_sub_zdb             = kaz_ge_sub_zdb_crw;
         kaz_ge_separate_zdb        = kaz_ge_separate_zdb_crw;

         %% Do a 3x3 pixel filter
         
         %disp('pause')
         %pause
         % Define the cnt_threshold
         %cnt_threshold        = 3;
         %cnt_threshold_edge   = 2;
         
         %[kaz_ge_single_Vmeany]            = func_remove_isolated_pixels_3x3_threshold(...
         %   kaz_ge_single_Vmean, cnt_threshold, cnt_threshold_edge);
         %[kaz_ge_single_Vmean]            = func_remove_isolated_pixels_5x5(kaz_ge_single_Vmean);
         [kaz_ge_single_3x3_filter]        = func_remove_isolated_pixels_3x3(kaz_ge_single_Vmean);
         
         %% Set reflectivity in lower gates to NaN
         
%         bad_gate_range                            = 1:2;
%         kaz_ge_single_zdb(:,bad_gate_range)       = kaz_ge_single_zdb(:,bad_gate_range) .* NaN;
%         kaz_ge_sub_zdb(:,bad_gate_range,:)        = kaz_ge_sub_zdb(:,bad_gate_range,:) .* NaN;
%         kaz_ge_separate_zdb(:,bad_gate_range,:)   = kaz_ge_separate_zdb(:,bad_gate_range,:) .* NaN;         
         
         %% NaN out all other variables using the 3x3 filtered variable
         
         [m,n,p]    = size(kaz_ge_sub_Vmean);
         
         for r = 1:m
            
            if(rem(r,1000) == 0)
               disp(['NaNing out variables, row: ',num2str(r),', out of ',num2str(m),' rows...']);
            end % end if(rem(r,60) == 0)
            
            % find the range gates with NaNs
            f                                = isnan(kaz_ge_single_Vmean(r,:));
            
            % Process the single variables
            if(sum(f) > 0)
               
               kaz_ge_single_Ppeak(r,f)           = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vkurt(r,f)           = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vleft(r,f)           = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vright(r,f)          = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vsig(r,f)            = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vskew(r,f)           = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vslope_lt(r,f)       = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vslope_rt(r,f)       = ones(1,sum(f)) .* NaN;
               kaz_ge_single_index_lt(r,f)        = ones(1,sum(f)) .* NaN;
               kaz_ge_single_index_rt(r,f)        = ones(1,sum(f)) .* NaN;
               %kaz_ge_single_max_1bin_delta(r,f)  = ones(1,sum(f)) .* NaN;
               %kaz_ge_single_max_2bin_delta(r,f)  = ones(1,sum(f)) .* NaN;
               kaz_ge_single_HS_nos_thres(r,f)    = ones(1,sum(f)) .* NaN;
               kaz_ge_single_nos_mean(r,f)        = ones(1,sum(f)) .* NaN;
               kaz_ge_single_nos_std(r,f)         = ones(1,sum(f)) .* NaN;
               kaz_ge_single_snr(r,f)             = ones(1,sum(f)) .* NaN;
               %kaz_ge_single_std_1bin(r,f)        = ones(1,sum(f)) .* NaN;
               %kaz_ge_single_std_2bin(r,f)        = ones(1,sum(f)) .* NaN;
               kaz_ge_single_zdb(r,f)             = ones(1,sum(f)) .* NaN;
               kaz_ge_single_Vmean_var(r,f)       = ones(1,sum(f)) .* NaN;
               kaz_ge_single_cnt_Vmean_var(r,f)   = ones(1,sum(f)) .* NaN;
               
               % Process the sub & separate variables
               
               for k = 1:p
                  
                  kaz_ge_sub_Ppeak(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_Vkurt(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_Vleft(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_Vmean(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_Vpeak(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_Vright(r,f,k)              = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_Vsig(r,f,k)                = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_Vskew(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_snr(r,f,k)                 = ones(1,sum(f),1) .* NaN;
                  kaz_ge_sub_zdb(r,f,k)                 = ones(1,sum(f),1) .* NaN;

                  kaz_ge_separate_Ppeak(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_Vkurt(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_Vleft(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_Vmean(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_Vpeak(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_Vright(r,f,k)              = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_Vsig(r,f,k)                = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_Vskew(r,f,k)               = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_snr(r,f,k)                 = ones(1,sum(f),1) .* NaN;
                  kaz_ge_separate_zdb(r,f,k)                 = ones(1,sum(f),1) .* NaN;
                  
               end % end for k loop
               
            end % end if(sum(f) > 0)
            
         end % end for r loop
         
         %disp('paused...')
         %pause
         
         %% Do the actual netCDF writing...
         
         % define the number of peaks to save
         max_num_peaks = 2;
         
         % Set the netCDF compression level.
         % The deflate_level determines the amount of netCDF compression. Value 
         % between 0 and 9, where 0 is no compression and 9 is the most compression.
         deflate_level  = 9;

         % The filename will be:
         netCDF_filename_root    = [site_prefix,'kazrge_15secshiftmom_cwilliams_'];
         netCDF_directory        = '..\nc_daily_15sec_ave_moments\';
         
         %netCDF_directory        = 'K:\multi_peak_spc\ge_copol\nc_files\';
         
         netCDF_output_filename         = [netCDF_directory,netCDF_filename_root,year_str,month_num_str,dayofmonth_str,version_suffex,'.nc'];
         
         [done] = func_save_ge_all_mom_in_netCDF(netCDF_output_filename, ...
            deflate_level, ...
            kaz_ge_single_year,           kaz_ge_single_month,...
            kaz_ge_single_dayofmonth,     kaz_ge_single_hour,...
            kaz_ge_single_minute,         kaz_ge_single_second,...
            kaz_ge_single_millisecond,    ARM_time,...
            base_time,                    time_offset,...
            ARM_midnight,                 time_from_midnight,...
            kaz_ge_single_range,...
            kaz_ge_single_range_bounds,   kaz_ge_single_Ppeak,...
            kaz_ge_single_Vkurt,          kaz_ge_single_Vleft,...
            kaz_ge_single_Vmean,          kaz_ge_single_Vpeak,...
            kaz_ge_single_Vright,         kaz_ge_single_Vsig,...
            kaz_ge_single_Vskew,          kaz_ge_single_Vslope_lt,...
            kaz_ge_single_Vslope_rt,      kaz_ge_single_index_lt,...
            kaz_ge_single_index_rt,       kaz_ge_single_HS_nos_thres,...
            kaz_ge_single_nos_mean,       kaz_ge_single_nos_std,...
            kaz_ge_single_snr,            kaz_ge_single_zdb,...
            kaz_ge_single_Vmean_var,      kaz_ge_single_cnt_Vmean_var,...
            kaz_ge_sub_Ppeak,             kaz_ge_sub_Vkurt,...
            kaz_ge_sub_Vleft,             kaz_ge_sub_Vmean,...
            kaz_ge_sub_Vpeak,             kaz_ge_sub_Vright,...
            kaz_ge_sub_Vsig,              kaz_ge_sub_Vskew,...
            kaz_ge_sub_snr,               kaz_ge_sub_zdb,...
            kaz_ge_separate_Ppeak,        kaz_ge_separate_Vkurt,...
            kaz_ge_separate_Vleft,        kaz_ge_separate_Vmean,...
            kaz_ge_separate_Vpeak,        kaz_ge_separate_Vright,...
            kaz_ge_separate_Vsig,         kaz_ge_separate_Vskew,...
            kaz_ge_separate_snr,          kaz_ge_separate_zdb,...
            max_num_peaks,                kaz_ge_single_3x3_filter);

         %% Clear out the memory
         
         % The size of the matrices change as the modes change.
         % So, need to clear variables defined within code.
         
         clear   kaz_ge_single_year           kaz_ge_single_month
         clear   kaz_ge_single_dayofmonth     kaz_ge_single_hour
         clear   kaz_ge_single_minute         kaz_ge_single_second
         clear   kaz_ge_single_millisecond    kaz_ge_single_range
         clear   kaz_ge_single_range_bounds   kaz_ge_single_Ppeak
         clear   kaz_ge_single_Vkurt          kaz_ge_single_Vleft
         clear   kaz_ge_single_Vmean          kaz_ge_single_Vpeak
         clear   kaz_ge_single_Vright         kaz_ge_single_Vsig
         clear   kaz_ge_single_Vskew          kaz_ge_single_Vslope_lt
         clear   kaz_ge_single_Vslope_rt      kaz_ge_single_index_lt
         clear   kaz_ge_single_index_rt       kaz_ge_single_HS_nos_thres
         clear   kaz_ge_single_nos_mean       kaz_ge_single_nos_std
         clear   kaz_ge_single_snr            kaz_ge_single_zdb
         clear   kaz_ge_single_Vmean_var      kaz_ge_single_cnt_Vmean_var
         clear   kaz_ge_sub_Ppeak             kaz_ge_sub_Vkurt
         clear   kaz_ge_sub_Vleft             kaz_ge_sub_Vmean
         clear   kaz_ge_sub_Vpeak             kaz_ge_sub_Vright
         clear   kaz_ge_sub_Vsig              kaz_ge_sub_Vskew
         clear   kaz_ge_sub_snr               kaz_ge_sub_zdb
         clear   kaz_ge_separate_Ppeak        kaz_ge_separate_Vkurt
         clear   kaz_ge_separate_Vleft        kaz_ge_separate_Vmean
         clear   kaz_ge_separate_Vpeak        kaz_ge_separate_Vright
         clear   kaz_ge_separate_Vsig         kaz_ge_separate_Vskew
         clear   kaz_ge_separate_snr          kaz_ge_separate_zdb
         clear   kaz_ge_single_3x3_filter

         
      end % end if(good_filename == 2)
      
   end % end for day_loop
end % end for month_loop

