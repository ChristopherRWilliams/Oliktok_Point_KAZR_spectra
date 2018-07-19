function [done] = func_save_ge_all_mom_in_netCDF(netCDF_output_filename, ...
   deflate_level,...
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
   kaz_ge_separate_snr,          kaz_ge_separate_zdb, ...
   max_num_peaks,                kaz_ge_single_3x3_filter)


% This routine saves the single peak moments in netCDF format

% updated: 09-February-2018
% ======================================================================

%% steps:

% 0. The variables to be saved are:
% 1. Define the bad data flag
% 2. Get the Dimensions of the arrays
% 3. Set all NaNs to bad_flag
% 4. Open the netCDF file
% 5. Define the dimensions in the output file
% 6. Define Variable names
% 7. Write the data to the output file
% 8. Add attributes to each variable
% 9. close the data file

%% 0. The variables to be saved are:

%   kaz_ge_single_year                42593x1                 340744  double
%   kaz_ge_single_month               42593x1                 340744  double
%   kaz_ge_single_dayofmonth          42593x1                 340744  double
%   kaz_ge_single_hour                42593x1                 340744  double
%   kaz_ge_single_minute              42593x1                 340744  double
%   kaz_ge_single_second              42593x1                 340744  double
%   kaz_ge_single_millisecond         42593x1                 340744  double

%   kaz_ge_single_range                   1x607                 4856  double
%   kaz_ge_single_range_bounds            2x607                 9712  double

%   kaz_ge_single_3x3_filter          42593x607            206831608  double
%   kaz_ge_single_HS_nos_thres        42593x607            206831608  double
%   kaz_ge_single_Ppeak               42593x607            206831608  double
%   kaz_ge_single_Vkurt               42593x607            206831608  double
%   kaz_ge_single_Vleft               42593x607            206831608  double
%   kaz_ge_single_Vmean               42593x607            206831608  double
%   kaz_ge_single_Vpeak               42593x607            206831608  double
%   kaz_ge_single_Vright              42593x607            206831608  double
%   kaz_ge_single_Vsig                42593x607            206831608  double
%   kaz_ge_single_Vskew               42593x607            206831608  double
%   kaz_ge_single_Vslope_lt           42593x607            206831608  double
%   kaz_ge_single_Vslope_rt           42593x607            206831608  double
%   kaz_ge_single_index_lt            42593x607            206831608  double
%   kaz_ge_single_index_rt            42593x607            206831608  double
%   kaz_ge_single_max_1bin_delta      42593x607            206831608  double
%   kaz_ge_single_max_2bin_delta      42593x607            206831608  double
%   kaz_ge_single_nos_mean            42593x607            206831608  double
%   kaz_ge_single_nos_std             42593x607            206831608  double
%   kaz_ge_single_snr                 42593x607            206831608  double
%   kaz_ge_single_std_1bin            42593x607            206831608  double
%   kaz_ge_single_std_2bin            42593x607            206831608  double
%   kaz_ge_single_zdb                 42593x607            206831608  double
%   kaz_ge_single_Vmean_var             5760x608               28016640  double
%   kaz_ge_single_cnt_Vmean_var         5760x608               28016640  double

%   kaz_ge_separate_Ppeak               5760x608x4            112066560  double
%   kaz_ge_separate_Vkurt               5760x608x4            112066560  double
%   kaz_ge_separate_Vleft               5760x608x4            112066560  double
%   kaz_ge_separate_Vmean               5760x608x4            112066560  double
%   kaz_ge_separate_Vpeak               5760x608x4            112066560  double
%   kaz_ge_separate_Vright              5760x608x4            112066560  double
%   kaz_ge_separate_Vsig                5760x608x4            112066560  double
%   kaz_ge_separate_Vskew               5760x608x4            112066560  double
%   kaz_ge_separate_snr                 5760x608x4            112066560  double
%   kaz_ge_separate_zdb                 5760x608x4            112066560  double

%% 1. Define the bad data flag (bad_flag)

bad_flag = -99;

%% 2. Get the Dimensions of the arrays

[ntime, nrange] = size(kaz_ge_single_snr);

%% 3. Set all NaNs to bad_flag - single matrices

g = isnan(kaz_ge_single_HS_nos_thres);
if(sum(sum(g)) > 0)
   kaz_ge_single_HS_nos_thres(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Ppeak);
if(sum(sum(g)) > 0)
   kaz_ge_single_Ppeak(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vkurt);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vkurt(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vleft);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vleft(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vmean);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vmean(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vpeak);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vpeak(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vright);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vright(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vsig);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vsig(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vskew);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vskew(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vslope_lt);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vslope_lt(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vslope_rt);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vslope_rt(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_index_lt);
if(sum(sum(g)) > 0)
   kaz_ge_single_index_lt(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_index_rt);
if(sum(sum(g)) > 0)
   kaz_ge_single_index_rt(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_nos_mean);
if(sum(sum(g)) > 0)
   kaz_ge_single_nos_mean(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_nos_std);
if(sum(sum(g)) > 0)
   kaz_ge_single_nos_std(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_snr);
if(sum(sum(g)) > 0)
   kaz_ge_single_snr(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_zdb);
if(sum(sum(g)) > 0)
   kaz_ge_single_zdb(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_Vmean_var);
if(sum(sum(g)) > 0)
   kaz_ge_single_Vmean_var(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_cnt_Vmean_var);
if(sum(sum(g)) > 0)
   kaz_ge_single_cnt_Vmean_var(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

g = isnan(kaz_ge_single_3x3_filter);
if(sum(sum(g)) > 0)
   kaz_ge_single_3x3_filter(g) = ones(sum(sum(g)),1) .* bad_flag;
end % end if(sum(sum(g)) > 0);

%% 3. Set all NaNs to bad_flag - sub peaks

g = isnan(kaz_ge_sub_Ppeak);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Ppeak(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_Vkurt);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Vkurt(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_Vleft);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Vleft(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_Vmean);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Vmean(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_Vpeak);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Vpeak(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_Vright);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Vright(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_Vsig);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Vsig(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_Vskew);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_Vskew(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_snr);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_snr(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_sub_zdb);
if(sum(sum(sum(g))) > 0)
   kaz_ge_sub_zdb(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

%% 3. Set all NaNs to bad_flag - separate peaks

g = isnan(kaz_ge_separate_Ppeak);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Ppeak(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_Vkurt);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Vkurt(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_Vleft);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Vleft(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_Vmean);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Vmean(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_Vpeak);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Vpeak(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_Vright);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Vright(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_Vsig);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Vsig(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_Vskew);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_Vskew(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_snr);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_snr(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

g = isnan(kaz_ge_separate_zdb);
if(sum(sum(sum(g))) > 0)
   kaz_ge_separate_zdb(g) = ones(sum(sum(sum(g))),1) .* bad_flag;
end % end if(sum(sum(sum(g))) > 0);

%% 4. Open the netCDF file

%ncid = netcdf.create(netCDF_output_filename, 'NC_CLOBBER');
ncid = netcdf.create(netCDF_output_filename, 'NETCDF4');

%% 5. Define the dimensions in the output file

time_dim_ID          = netcdf.defDim(ncid,'time',ntime);
height_dim_ID        = netcdf.defDim(ncid,'range',nrange);
%npts_dim_ID     = netcdf.defDim(ncid,'spectrum_length',npts);
single_dim_ID        = netcdf.defDim(ncid,'single_value',1);
max_num_peak_dim_ID  = netcdf.defDim(ncid,'max_num_peak',max_num_peaks);

%% 6. Define Variable names

% single value variables
badflag_ID             = netcdf.defVar(ncid,'kaz_ge_single_bad_flag','float',single_dim_ID);

% Range dimension variables
range_ID               = netcdf.defVar(ncid,'kaz_ge_single_range','float',height_dim_ID);
lower_range_bound_ID   = netcdf.defVar(ncid,'kaz_ge_single_lower_range_bound','float',height_dim_ID);
upper_range_bound_ID   = netcdf.defVar(ncid,'kaz_ge_single_upper_range_bound','float',height_dim_ID);

% Time dimension variabiles
year_ID                = netcdf.defVar(ncid,'kaz_ge_single_year','float',time_dim_ID);
month_ID               = netcdf.defVar(ncid,'kaz_ge_single_month','float',time_dim_ID);
dayofmonth_ID          = netcdf.defVar(ncid,'kaz_ge_single_dayofmonth','float',time_dim_ID);
hour_ID                = netcdf.defVar(ncid,'kaz_ge_single_hour','float',time_dim_ID);
minute_ID              = netcdf.defVar(ncid,'kaz_ge_single_minute','float',time_dim_ID);
second_ID              = netcdf.defVar(ncid,'kaz_ge_single_second','float',time_dim_ID);
millisecond_ID         = netcdf.defVar(ncid,'kaz_ge_single_millisecond','float',time_dim_ID);

base_time_ID           = netcdf.defVar(ncid,'base_time','float',single_dim_ID);
time_offset_ID         = netcdf.defVar(ncid,'time_offset','float',time_dim_ID);
time_from_midnight_ID  = netcdf.defVar(ncid,'time_from_midnight','float',time_dim_ID);

% Time and Range dimension variables
kaz_ge_single_HS_nos_thres_ID       = netcdf.defVar(ncid,'kaz_ge_single_HS_nos_thres','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Ppeak_ID              = netcdf.defVar(ncid,'kaz_ge_single_Ppeak','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vkurt_ID              = netcdf.defVar(ncid,'kaz_ge_single_Vkurt','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vleft_ID              = netcdf.defVar(ncid,'kaz_ge_single_Vleft','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vmean_ID              = netcdf.defVar(ncid,'kaz_ge_single_Vmean','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vpeak_ID              = netcdf.defVar(ncid,'kaz_ge_single_Vpeak','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vright_ID             = netcdf.defVar(ncid,'kaz_ge_single_Vright','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vsig_ID               = netcdf.defVar(ncid,'kaz_ge_single_Vsig','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vskew_ID              = netcdf.defVar(ncid,'kaz_ge_single_Vskew','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vslope_lt_ID          = netcdf.defVar(ncid,'kaz_ge_single_Vslope_lt','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vslope_rt_ID          = netcdf.defVar(ncid,'kaz_ge_single_Vslope_rt','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_index_lt_ID           = netcdf.defVar(ncid,'kaz_ge_single_index_lt','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_index_rt_ID           = netcdf.defVar(ncid,'kaz_ge_single_index_rt','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_nos_mean_ID           = netcdf.defVar(ncid,'kaz_ge_single_nos_mean','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_nos_std_ID            = netcdf.defVar(ncid,'kaz_ge_single_nos_std','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_snr_ID                = netcdf.defVar(ncid,'kaz_ge_single_snr','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_zdb_ID                = netcdf.defVar(ncid,'kaz_ge_single_zdb','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_Vmean_var_ID          = netcdf.defVar(ncid,'kaz_ge_single_Vmean_var','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_cnt_Vmean_var_ID      = netcdf.defVar(ncid,'kaz_ge_single_cnt_Vmean_var','float',[time_dim_ID height_dim_ID]);
kaz_ge_single_3x3_filter_ID         = netcdf.defVar(ncid,'kaz_ge_single_3x3_filter','float',[time_dim_ID height_dim_ID]);

kaz_ge_sub_Ppeak_ID            = netcdf.defVar(ncid,'kaz_ge_sub_Ppeak','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_Vkurt_ID            = netcdf.defVar(ncid,'kaz_ge_sub_Vkurt','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_Vleft_ID            = netcdf.defVar(ncid,'kaz_ge_sub_Vleft','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_Vmean_ID            = netcdf.defVar(ncid,'kaz_ge_sub_Vmean','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_Vpeak_ID            = netcdf.defVar(ncid,'kaz_ge_sub_Vpeak','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_Vright_ID           = netcdf.defVar(ncid,'kaz_ge_sub_Vright','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_Vsig_ID             = netcdf.defVar(ncid,'kaz_ge_sub_Vsig','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_Vskew_ID            = netcdf.defVar(ncid,'kaz_ge_sub_Vskew','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_snr_ID              = netcdf.defVar(ncid,'kaz_ge_sub_snr','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_sub_zdb_ID              = netcdf.defVar(ncid,'kaz_ge_sub_zdb','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);

kaz_ge_separate_Ppeak_ID            = netcdf.defVar(ncid,'kaz_ge_separate_Ppeak','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_Vkurt_ID            = netcdf.defVar(ncid,'kaz_ge_separate_Vkurt','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_Vleft_ID            = netcdf.defVar(ncid,'kaz_ge_separate_Vleft','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_Vmean_ID            = netcdf.defVar(ncid,'kaz_ge_separate_Vmean','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_Vpeak_ID            = netcdf.defVar(ncid,'kaz_ge_separate_Vpeak','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_Vright_ID           = netcdf.defVar(ncid,'kaz_ge_separate_Vright','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_Vsig_ID             = netcdf.defVar(ncid,'kaz_ge_separate_Vsig','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_Vskew_ID            = netcdf.defVar(ncid,'kaz_ge_separate_Vskew','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_snr_ID              = netcdf.defVar(ncid,'kaz_ge_separate_snr','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);
kaz_ge_separate_zdb_ID              = netcdf.defVar(ncid,'kaz_ge_separate_zdb','float',[time_dim_ID height_dim_ID max_num_peak_dim_ID]);

%% Define the compression parameters for each netCDF variable

% The deflate_level determines the amount of netCDF compression. Value 
% between 0 and 9, where 0 is no compression and 9 is the most compression.
% deflate_level is an input parameter

netcdf.defVarDeflate(ncid, badflag_ID, true, true, deflate_level);

netcdf.defVarDeflate(ncid, range_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, lower_range_bound_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, upper_range_bound_ID, true, true, deflate_level);

netcdf.defVarDeflate(ncid, year_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, month_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, dayofmonth_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, hour_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, minute_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, second_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, millisecond_ID, true, true, deflate_level);

netcdf.defVarDeflate(ncid, base_time_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, time_offset_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, time_from_midnight_ID, true, true, deflate_level);

netcdf.defVarDeflate(ncid, kaz_ge_single_HS_nos_thres_ID, true, true, deflate_level);
netcdf.defVarDeflate(ncid, kaz_ge_single_Ppeak_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vkurt_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vleft_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vmean_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vpeak_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vright_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vsig_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vskew_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vslope_lt_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vslope_rt_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_index_lt_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_index_rt_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_nos_mean_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_nos_std_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_snr_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_zdb_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_Vmean_var_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_cnt_Vmean_var_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_single_3x3_filter_ID, true, true, deflate_level); 

netcdf.defVarDeflate(ncid, kaz_ge_sub_Ppeak_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_Vkurt_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_Vleft_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_Vmean_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_Vpeak_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_Vright_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_Vsig_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_Vskew_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_snr_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_sub_zdb_ID, true, true, deflate_level); 

netcdf.defVarDeflate(ncid, kaz_ge_separate_Ppeak_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_Vkurt_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_Vleft_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_Vmean_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_Vpeak_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_Vright_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_Vsig_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_Vskew_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_snr_ID, true, true, deflate_level); 
netcdf.defVarDeflate(ncid, kaz_ge_separate_zdb_ID, true, true, deflate_level); 

%% 7. Write the data to the output file

% Need to get out of define mode before writting data
netcdf.endDef(ncid);

% Write the single value variables
netcdf.putVar(ncid,badflag_ID,single(bad_flag));

% Write the Range dimension variables
netcdf.putVar(ncid,range_ID,single(kaz_ge_single_range));
netcdf.putVar(ncid,lower_range_bound_ID,single(kaz_ge_single_range_bounds(1,:)));
netcdf.putVar(ncid,upper_range_bound_ID,single(kaz_ge_single_range_bounds(2,:)));

% Write the Time dimension variabiles
netcdf.putVar(ncid,year_ID,single(kaz_ge_single_year));
netcdf.putVar(ncid,month_ID,single(kaz_ge_single_month));
netcdf.putVar(ncid,dayofmonth_ID,single(kaz_ge_single_dayofmonth));
netcdf.putVar(ncid,hour_ID,single(kaz_ge_single_hour));
netcdf.putVar(ncid,minute_ID,single(kaz_ge_single_minute));
netcdf.putVar(ncid,second_ID,single(kaz_ge_single_second));
netcdf.putVar(ncid,millisecond_ID,single(kaz_ge_single_millisecond));

netcdf.putVar(ncid,base_time_ID,single(base_time));
netcdf.putVar(ncid,time_offset_ID,single(time_offset));
netcdf.putVar(ncid,time_from_midnight_ID,single(time_from_midnight));

% Write the Time and Range dimension variables
netcdf.putVar(ncid,kaz_ge_single_HS_nos_thres_ID,  single(kaz_ge_single_HS_nos_thres));
netcdf.putVar(ncid,kaz_ge_single_Ppeak_ID,         single(kaz_ge_single_Ppeak));
netcdf.putVar(ncid,kaz_ge_single_Vkurt_ID,         single(kaz_ge_single_Vkurt));
netcdf.putVar(ncid,kaz_ge_single_Vleft_ID,         single(kaz_ge_single_Vleft));
netcdf.putVar(ncid,kaz_ge_single_Vmean_ID,         single(kaz_ge_single_Vmean));
netcdf.putVar(ncid,kaz_ge_single_Vpeak_ID,         single(kaz_ge_single_Vpeak));
netcdf.putVar(ncid,kaz_ge_single_Vright_ID,        single(kaz_ge_single_Vright));
netcdf.putVar(ncid,kaz_ge_single_Vsig_ID,          single(kaz_ge_single_Vsig));
netcdf.putVar(ncid,kaz_ge_single_Vskew_ID,         single(kaz_ge_single_Vskew));
netcdf.putVar(ncid,kaz_ge_single_Vslope_lt_ID,     single(kaz_ge_single_Vslope_lt));
netcdf.putVar(ncid,kaz_ge_single_Vslope_rt_ID,     single(kaz_ge_single_Vslope_rt));
netcdf.putVar(ncid,kaz_ge_single_index_lt_ID,      single(kaz_ge_single_index_lt));
netcdf.putVar(ncid,kaz_ge_single_index_rt_ID,      single(kaz_ge_single_index_rt));
netcdf.putVar(ncid,kaz_ge_single_nos_mean_ID,      single(kaz_ge_single_nos_mean));
netcdf.putVar(ncid,kaz_ge_single_nos_std_ID,       single(kaz_ge_single_nos_std));
netcdf.putVar(ncid,kaz_ge_single_snr_ID,           single(kaz_ge_single_snr));
netcdf.putVar(ncid,kaz_ge_single_zdb_ID,           single(kaz_ge_single_zdb));

netcdf.putVar(ncid,kaz_ge_single_Vmean_var_ID,  single(kaz_ge_single_Vmean_var));
netcdf.putVar(ncid,kaz_ge_single_cnt_Vmean_var_ID, single(kaz_ge_single_cnt_Vmean_var));

netcdf.putVar(ncid,kaz_ge_single_3x3_filter_ID, single(kaz_ge_single_3x3_filter));

netcdf.putVar(ncid,kaz_ge_sub_Ppeak_ID,      single(kaz_ge_sub_Ppeak(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_Vkurt_ID,      single(kaz_ge_sub_Vkurt(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_Vleft_ID,      single(kaz_ge_sub_Vleft(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_Vmean_ID,      single(kaz_ge_sub_Vmean(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_Vpeak_ID,      single(kaz_ge_sub_Vpeak(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_Vright_ID,     single(kaz_ge_sub_Vright(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_Vsig_ID,       single(kaz_ge_sub_Vsig(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_Vskew_ID,      single(kaz_ge_sub_Vskew(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_sub_snr_ID,        single(kaz_ge_sub_snr(:,:,1:max_num_peaks)));
%netcdf.putVar(ncid,kaz_ge_sub_zdb_ID,        single(permute(kaz_ge_sub_zdb,[2 1 3])));
netcdf.putVar(ncid,kaz_ge_sub_zdb_ID,        single(kaz_ge_sub_zdb(:,:,1:max_num_peaks)));

netcdf.putVar(ncid,kaz_ge_separate_Ppeak_ID,      single(kaz_ge_separate_Ppeak(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_Vkurt_ID,      single(kaz_ge_separate_Vkurt(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_Vleft_ID,      single(kaz_ge_separate_Vleft(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_Vmean_ID,      single(kaz_ge_separate_Vmean(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_Vpeak_ID,      single(kaz_ge_separate_Vpeak(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_Vright_ID,     single(kaz_ge_separate_Vright(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_Vsig_ID,       single(kaz_ge_separate_Vsig(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_Vskew_ID,      single(kaz_ge_separate_Vskew(:,:,1:max_num_peaks)));
netcdf.putVar(ncid,kaz_ge_separate_snr_ID,        single(kaz_ge_separate_snr(:,:,1:max_num_peaks)));
%netcdf.putVar(ncid,kaz_ge_separate_zdb_ID,        single(permute(kaz_ge_separate_zdb,[1 2 3])));
netcdf.putVar(ncid,kaz_ge_separate_zdb_ID,        single(kaz_ge_separate_zdb(:,:,1:max_num_peaks)));

%% 8. Add attributes to each variable

% Go back into the define mode
netcdf.reDef(ncid);

% Put attributes in the order in which the variables were written

netcdf.putAtt(ncid,0,'long_name','Bad Data Flag, -99');
netcdf.putAtt(ncid,0,'units','integer');

netcdf.putAtt(ncid,1,'long_name','Range above the radar');
netcdf.putAtt(ncid,1,'units','meter');

netcdf.putAtt(ncid,2,'long_name','Lower range bound above the radar');
netcdf.putAtt(ncid,2,'units','meter');
netcdf.putAtt(ncid,3,'long_name','Upper range bound above the radar');
netcdf.putAtt(ncid,3,'units','meter');

netcdf.putAtt(ncid,4,'long_name','Year');
netcdf.putAtt(ncid,4,'units','ADC');
netcdf.putAtt(ncid,5,'long_name','Month of Year');
netcdf.putAtt(ncid,5,'units','day');
netcdf.putAtt(ncid,6,'long_name','Day of Month');
netcdf.putAtt(ncid,6,'units','day');
netcdf.putAtt(ncid,7,'long_name','Hour of Day');
netcdf.putAtt(ncid,7,'units','hour (UTC)');
netcdf.putAtt(ncid,8,'long_name','Minute of Hour');
netcdf.putAtt(ncid,8,'units','minute');
netcdf.putAtt(ncid,9,'long_name','Second');
netcdf.putAtt(ncid,9,'units','second');
netcdf.putAtt(ncid,10,'long_name','Millisecond');
netcdf.putAtt(ncid,10,'units','millisecond');

netcdf.putAtt(ncid,11,'long_name','Base time in Epoch');
netcdf.putAtt(ncid,11,'units','seconds since 1970-1-1 0:00:00 0:00');
netcdf.putAtt(ncid,11,'String',datestr(ARM_time(1),'yyyy-mm-dd HH:MM:SS'));

netcdf.putAtt(ncid,12,'long_name','Time offset from base time');
netcdf.putAtt(ncid,12,'units',['seconds since ',datestr(ARM_time(1),'yyyy-mm-dd HH:MM:SS'),' 0:00']);

netcdf.putAtt(ncid,13,'long_name','Time offset from midnight');
netcdf.putAtt(ncid,13,'units',['seconds since ',datestr(ARM_midnight(1),'yyyy-mm-dd HH:MM:SS'),' 0:00']);

% netcdf.putVar(ncid,kaz_ge_single_HS_nos_thres_ID,single(kaz_ge_single_HS_nos_thres'));
% netcdf.putVar(ncid,kaz_ge_single_Ppeak_ID,      single(kaz_ge_single_Ppeak'));
% netcdf.putVar(ncid,kaz_ge_single_Vkurt_ID,      single(kaz_ge_single_Vkurt'));
% netcdf.putVar(ncid,kaz_ge_single_Vleft_ID,      single(kaz_ge_single_Vleft'));
netcdf.putAtt(ncid,14,'long_name','Hildebrand & Sekhon (1974) noise threshold');
netcdf.putAtt(ncid,14,'units','dBm');
netcdf.putAtt(ncid,15,'long_name','Power in spectral bin at peak of spectrum');
netcdf.putAtt(ncid,15,'units','dBm');
netcdf.putAtt(ncid,16,'long_name','Power spectrum weighted velocity Kurtosis (positive values are downward)');
netcdf.putAtt(ncid,16,'units','(m/s)^4');
netcdf.putAtt(ncid,17,'long_name','Power spectrum left edge velocity (upward side of spectrum)');
netcdf.putAtt(ncid,17,'units','m/s');

% netcdf.putVar(ncid,kaz_ge_single_Vmean_ID,      single(kaz_ge_single_Vmean'));
% netcdf.putVar(ncid,kaz_ge_single_Vpeak_ID,      single(kaz_ge_single_Vpeak'));
% netcdf.putVar(ncid,kaz_ge_single_Vright_ID,     single(kaz_ge_single_Vright'));
% netcdf.putVar(ncid,kaz_ge_single_Vsig_ID,       single(kaz_ge_single_Vsig'));
% netcdf.putVar(ncid,kaz_ge_single_Vskew_ID,      single(kaz_ge_single_Vskew'));
netcdf.putAtt(ncid,18,'long_name','Power spectrum weighted mean velocity (positive values are downward)');
netcdf.putAtt(ncid,18,'units','m/s');
netcdf.putAtt(ncid,19,'long_name','Velocity of spectral bin at peak of spectrum (positive values are downward)');
netcdf.putAtt(ncid,19,'units','m/s');
netcdf.putAtt(ncid,20,'long_name','Power spectrum right edge velocity (downward side of spectrum)');
netcdf.putAtt(ncid,20,'units','m/s');
netcdf.putAtt(ncid,21,'long_name','Power spectrum weighted velocity standard deviation');
netcdf.putAtt(ncid,21,'units','m/s');
netcdf.putAtt(ncid,22,'long_name','Power spectrum weighted velocity skewness (positive values are downward)');
netcdf.putAtt(ncid,22,'units','(m/s)^3');

% netcdf.putVar(ncid,kaz_ge_single_Vslope_lt_ID,  single(kaz_ge_single_Vslope_lt'));
% netcdf.putVar(ncid,kaz_ge_single_Vslope_rt_ID,  single(kaz_ge_single_Vslope_rt'));
% netcdf.putVar(ncid,kaz_ge_single_index_lt_ID,   single(kaz_ge_single_index_lt'));
% netcdf.putVar(ncid,kaz_ge_single_index_rt_ID,   single(kaz_ge_single_index_rt'));
netcdf.putAtt(ncid,23,'long_name','Slope from spectrum peak to left edge (upward side of spectrum)');
netcdf.putAtt(ncid,23,'units','dBm/(m/s)');
netcdf.putAtt(ncid,24,'long_name','Slope from spectrum peak to right edge (downward side of spectrum)');
netcdf.putAtt(ncid,24,'units','dBm/(m/s)');
netcdf.putAtt(ncid,25,'long_name','Spectral bin index at left edge of spectrum (upward side of spectrum), this is Vd(index)');
netcdf.putAtt(ncid,25,'units','index');
netcdf.putAtt(ncid,26,'long_name','Spectral bin index at right edge of spectrum (downward side of spectrum), this is Vd(index)');
netcdf.putAtt(ncid,26,'units','index');

netcdf.putAtt(ncid,27,'long_name','mean noise value');
netcdf.putAtt(ncid,27,'units','dBm');
netcdf.putAtt(ncid,28,'long_name','noise standard deviation');
netcdf.putAtt(ncid,28,'units','dBm');

% netcdf.putVar(ncid,kaz_ge_single_snr_ID,        single(kaz_ge_single_snr'));
% netcdf.putVar(ncid,kaz_ge_single_zdb_ID,        single(kaz_ge_single_zdb'));
netcdf.putAtt(ncid,29,'long_name','Power spectrum Signal-to-Noise Ratio');
netcdf.putAtt(ncid,29,'units','dB');
netcdf.putAtt(ncid,30,'long_name','Measured reflectivity factor');
netcdf.putAtt(ncid,30,'units','dBZ');

% netcdf.putVar(ncid,kaz_ge_single_Vmean_var_ID,  single(kaz_ge_single_Vmean_var'));
% netcdf.putVar(ncid,kaz_ge_single_cnt_Vmean_var_ID, single(kaz_ge_single_cnt_Vmean_var'));
netcdf.putAtt(ncid,31,'long_name','Velocity variance of individual Vmean values during dwell interval');
netcdf.putAtt(ncid,31,'units','m^2/s^2');
netcdf.putAtt(ncid,32,'long_name','Count of valid Vmean values during dwell interval used in Vmean_var (velocity variance) calculation');
netcdf.putAtt(ncid,32,'units','integer');

netcdf.putAtt(ncid,33,'long_name','3x3 filter flag. Flag = 0 means less than 3 neighors. Flag = 1 means 3 or more neighbors.');
netcdf.putAtt(ncid,33,'units','integer');

% netcdf.putVar(ncid,kaz_ge_sub_Ppeak_ID,      single(kaz_ge_sub_Ppeak'));
% netcdf.putVar(ncid,kaz_ge_sub_Vkurt_ID,      single(kaz_ge_sub_Vkurt'));
% netcdf.putVar(ncid,kaz_ge_sub_Vleft_ID,      single(kaz_ge_sub_Vleft'));
netcdf.putAtt(ncid,34,'long_name','Sub peak (defined as peak within single peak): Power in spectral bin at peak of spectrum');
netcdf.putAtt(ncid,34,'units','dBm');
netcdf.putAtt(ncid,35,'long_name','Sub peak (defined as peak within single peak): Power spectrum weighted velocity Kurtosis (positive values are downward)');
netcdf.putAtt(ncid,35,'units','(m/s)^4');
netcdf.putAtt(ncid,36,'long_name','Sub peak (defined as peak within single peak): Power spectrum left edge velocity (upward side of spectrum)');
netcdf.putAtt(ncid,36,'units','m/s');

% netcdf.putVar(ncid,kaz_ge_sub_Vmean_ID,      single(kaz_ge_sub_Vmean'));
% netcdf.putVar(ncid,kaz_ge_sub_Vpeak_ID,      single(kaz_ge_sub_Vpeak'));
% netcdf.putVar(ncid,kaz_ge_sub_Vright_ID,     single(kaz_ge_sub_Vright'));
% netcdf.putVar(ncid,kaz_ge_sub_Vsig_ID,       single(kaz_ge_sub_Vsig'));
% netcdf.putVar(ncid,kaz_ge_sub_Vskew_ID,      single(kaz_ge_sub_Vskew'));
netcdf.putAtt(ncid,37,'long_name','Sub peak (defined as peak within single peak): Power spectrum weighted mean velocity (positive values are downward)');
netcdf.putAtt(ncid,37,'units','m/s');
netcdf.putAtt(ncid,38,'long_name','Sub peak (defined as peak within single peak): Velocity of spectral bin at peak of spectrum (positive values are downward)');
netcdf.putAtt(ncid,38,'units','m/s');
netcdf.putAtt(ncid,39,'long_name','Sub peak (defined as peak within single peak): Power spectrum right edge velocity (downward side of spectrum)');
netcdf.putAtt(ncid,39,'units','m/s');
netcdf.putAtt(ncid,40,'long_name','Sub peak (defined as peak within single peak): Power spectrum weighted velocity standard deviation');
netcdf.putAtt(ncid,40,'units','m/s');
netcdf.putAtt(ncid,41,'long_name','Sub peak (defined as peak within single peak): Power spectrum weighted velocity skewness (positive values are downward)');
netcdf.putAtt(ncid,41,'units','(m/s)^3');

% netcdf.putVar(ncid,kaz_ge_sub_snr_ID,        single(kaz_ge_sub_snr'));
% netcdf.putVar(ncid,kaz_ge_sub_zdb_ID,        single(kaz_ge_sub_zdb'));
netcdf.putAtt(ncid,42,'long_name','Sub peak (defined as peak within single peak): Power spectrum Signal-to-Noise Ratio');
netcdf.putAtt(ncid,42,'units','dB');
netcdf.putAtt(ncid,43,'long_name','Sub peak (defined as peak within single peak): Measured reflectivity factor');
netcdf.putAtt(ncid,43,'units','dBZ');

% netcdf.putVar(ncid,kaz_ge_separate_Ppeak_ID,      single(kaz_ge_separate_Ppeak'));
% netcdf.putVar(ncid,kaz_ge_separate_Vkurt_ID,      single(kaz_ge_separate_Vkurt'));
% netcdf.putVar(ncid,kaz_ge_separate_Vleft_ID,      single(kaz_ge_separate_Vleft'));
netcdf.putAtt(ncid,44,'long_name','Separate peak (defined as noise between single peak and this peak): Power in spectral bin at peak of spectrum');
netcdf.putAtt(ncid,44,'units','dBm');
netcdf.putAtt(ncid,45,'long_name','Separate peak (defined as noise between single peak and this peak): Power spectrum weighted velocity Kurtosis (positive values are downward)');
netcdf.putAtt(ncid,45,'units','(m/s)^4');
netcdf.putAtt(ncid,46,'long_name','Separate peak (defined as noise between single peak and this peak): Power spectrum left edge velocity (upward side of spectrum)');
netcdf.putAtt(ncid,46,'units','m/s');

% netcdf.putVar(ncid,kaz_ge_separate_Vmean_ID,      single(kaz_ge_separate_Vmean'));
% netcdf.putVar(ncid,kaz_ge_separate_Vpeak_ID,      single(kaz_ge_separate_Vpeak'));
% netcdf.putVar(ncid,kaz_ge_separate_Vright_ID,     single(kaz_ge_separate_Vright'));
% netcdf.putVar(ncid,kaz_ge_separate_Vsig_ID,       single(kaz_ge_separate_Vsig'));
% netcdf.putVar(ncid,kaz_ge_separate_Vskew_ID,      single(kaz_ge_separate_Vskew'));
netcdf.putAtt(ncid,47,'long_name','Separate peak (defined as noise between single peak and this peak): Power spectrum weighted mean velocity (positive values are downward)');
netcdf.putAtt(ncid,47,'units','m/s');
netcdf.putAtt(ncid,48,'long_name','Separate peak (defined as noise between single peak and this peak): Velocity of spectral bin at peak of spectrum (positive values are downward)');
netcdf.putAtt(ncid,48,'units','m/s');
netcdf.putAtt(ncid,49,'long_name','Separate peak (defined as noise between single peak and this peak): Power spectrum right edge velocity (downward side of spectrum)');
netcdf.putAtt(ncid,49,'units','m/s');
netcdf.putAtt(ncid,50,'long_name','Separate peak (defined as noise between single peak and this peak): Power spectrum weighted velocity standard deviation');
netcdf.putAtt(ncid,50,'units','m/s');
netcdf.putAtt(ncid,51,'long_name','Separate peak (defined as noise between single peak and this peak): Power spectrum weighted velocity skewness (positive values are downward)');
netcdf.putAtt(ncid,51,'units','(m/s)^3');

% netcdf.putVar(ncid,kaz_ge_separate_snr_ID,        single(kaz_ge_separate_snr'));
% netcdf.putVar(ncid,kaz_ge_separate_zdb_ID,        single(kaz_ge_separate_zdb'));
netcdf.putAtt(ncid,52,'long_name','Separate peak (defined as noise between single peak and this peak): Power spectrum Signal-to-Noise Ratio');
netcdf.putAtt(ncid,52,'units','dB');
netcdf.putAtt(ncid,53,'long_name','Separate peak (defined as noise between single peak and this peak): Measured reflectivity factor');
netcdf.putAtt(ncid,53,'units','dBZ');

%% 9. close the data file

netcdf.close(ncid)

%% 10. Save a value

done = 1;

