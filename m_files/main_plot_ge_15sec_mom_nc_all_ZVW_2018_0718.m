% Program = 'main_plot_ge_15sec_mom_nc_all_ZVW_2017_1021.m'
% updated = '04-December-2017'

% This routine reads KAZR spectra and estimates moments.

% **********************************************************************

%%    1. Define some variables

site_prefix             = 'oli';
version_suffex          = '_v15';

max_vel                 = 8;  % define maximum velocity in color plot
max_Vsig                = 2;  % define max vel STD in color plot

hourly_flag             = 1;

%%    2. Define the hours to process

% Days to process
% 114, 24 April
% 115, 25 April
% 117, 27 April
% 131, 11 May
% 140, 20 May
% 143, 23 May
% 144, 24 May
%days_to_process  = [114 115 115 117 131 131 140 143 144];
%hour_start_range = [  9   8  18   5   0  17   6  21  20];
%hour_end_range   = [ 12  11  20  15   0  23  16  23  22];

start_year  = 2016;
start_month =    6;
end_month   =    6;
start_day   =   19;
end_day     =   19;
%start_hour  =    4;
%end_hour    =    4;

%%    3. Process each day

for month_loop = start_month:end_month
   
   for dayofmonth_loop = start_day:end_day
      
      
      %%    4. Define the file to be loaded
      
      year_str       = num2str(start_year);
      month_num_str  = num2str_2digits(month_loop);
      dayofmonth_str = num2str_2digits(dayofmonth_loop);
      
      input_directory         = '..\nc_daily_15sec_ave_moments\';
      filename_root           = [site_prefix,'kazrge_15secshiftmom_cwilliams_'];
      netCDF_input_filename   = [input_directory,filename_root,year_str,month_num_str,dayofmonth_str,version_suffex,'.nc'];
      
      good_input_filename  = exist(netCDF_input_filename,'file');
      
      output_daily_dir     = '..\images_from_nc_files\daily\';
      output_hourly_dir    = '..\images_from_nc_files\hourly\';
      
      %%    5. Load in this daily file
      
      if(good_input_filename == 2)
         
         % load moments
         %disp(['loading netCDF file: ',input_netCDF_filename,' ...']);
         disp(['loading nc file: ',netCDF_input_filename,' ...']);
         
         % the input file is named: input_netCDF_filename
         func_read_nc_file
         % the variables are:
         %   kaz_ge_separate_Ppeak                5760x607x2            55941120  double
         %   kaz_ge_separate_Vkurt                5760x607x2            55941120  double
         %   kaz_ge_separate_Vleft                5760x607x2            55941120  double
         %   kaz_ge_separate_Vmean                5760x607x2            55941120  double
         %   kaz_ge_separate_Vpeak                5760x607x2            55941120  double
         %   kaz_ge_separate_Vright               5760x607x2            55941120  double
         %   kaz_ge_separate_Vsig                 5760x607x2            55941120  double
         %   kaz_ge_separate_Vskew                5760x607x2            55941120  double
         %   kaz_ge_separate_snr                  5760x607x2            55941120  double
         %   kaz_ge_separate_zdb                  5760x607x2            55941120  double
         %   kaz_ge_single_3x3_filter             5760x607              27970560  double
         %   kaz_ge_single_HS_nos_thres           5760x607              27970560  double
         %   kaz_ge_single_Ppeak                  5760x607              27970560  double
         %   kaz_ge_single_Vkurt                  5760x607              27970560  double
         %   kaz_ge_single_Vleft                  5760x607              27970560  double
         %   kaz_ge_single_Vmean                  5760x607              27970560  double
         %   kaz_ge_single_Vmean_var              5760x607              27970560  double
         %   kaz_ge_single_Vpeak                  5760x607              27970560  double
         %   kaz_ge_single_Vright                 5760x607              27970560  double
         %   kaz_ge_single_Vsig                   5760x607              27970560  double
         %   kaz_ge_single_Vskew                  5760x607              27970560  double
         %   kaz_ge_single_Vslope_lt              5760x607              27970560  double
         %   kaz_ge_single_Vslope_rt              5760x607              27970560  double
         %   kaz_ge_single_bad_flag                  1x1                       8  double
         %   kaz_ge_single_cnt_Vmean_var          5760x607              27970560  double
         %   kaz_ge_single_dayofmonth             5760x1                   46080  double
         %   kaz_ge_single_hour                   5760x1                   46080  double
         %   kaz_ge_single_index_lt               5760x607              27970560  double
         %   kaz_ge_single_index_rt               5760x607              27970560  double
         %   kaz_ge_single_lower_range_bound       607x1                    4856  double
         %   kaz_ge_single_millisecond            5760x1                   46080  double
         %   kaz_ge_single_minute                 5760x1                   46080  double
         %   kaz_ge_single_month                  5760x1                   46080  double
         %   kaz_ge_single_nos_mean               5760x607              27970560  double
         %   kaz_ge_single_nos_std                5760x607              27970560  double
         %   kaz_ge_single_range                   607x1                    4856  double
         %   kaz_ge_single_second                 5760x1                   46080  double
         %   kaz_ge_single_snr                    5760x607              27970560  double
         %   kaz_ge_single_upper_range_bound       607x1                    4856  double
         %   kaz_ge_single_year                   5760x1                   46080  double
         %   kaz_ge_single_zdb                    5760x607              27970560  double
         %   kaz_ge_sub_Ppeak                     5760x607x2            55941120  double
         %   kaz_ge_sub_Vkurt                     5760x607x2            55941120  double
         %   kaz_ge_sub_Vleft                     5760x607x2            55941120  double
         %   kaz_ge_sub_Vmean                     5760x607x2            55941120  double
         %   kaz_ge_sub_Vpeak                     5760x607x2            55941120  double
         %   kaz_ge_sub_Vright                    5760x607x2            55941120  double
         %   kaz_ge_sub_Vsig                      5760x607x2            55941120  double
         %   kaz_ge_sub_Vskew                     5760x607x2            55941120  double
         %   kaz_ge_sub_snr                       5760x607x2            55941120  double
         %   kaz_ge_sub_zdb                       5760x607x2            55941120  double
  
         % this old version transposes matrices because of how nc file was written...
         %func_load_all_netCDF_variables
         % the variables are:
         
         %% Set all bad values to NaN
         
         kaz_ge_single_zdb(kaz_ge_single_zdb       == kaz_ge_single_bad_flag) = kaz_ge_single_zdb(kaz_ge_single_zdb == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vmean(kaz_ge_single_Vmean   == kaz_ge_single_bad_flag) = kaz_ge_single_Vmean(kaz_ge_single_Vmean == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vsig(kaz_ge_single_Vsig     == kaz_ge_single_bad_flag) = kaz_ge_single_Vsig(kaz_ge_single_Vsig == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vskew(kaz_ge_single_Vskew   == kaz_ge_single_bad_flag) = kaz_ge_single_Vskew(kaz_ge_single_Vskew == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_snr(kaz_ge_single_snr       == kaz_ge_single_bad_flag) = kaz_ge_single_snr(kaz_ge_single_snr == kaz_ge_single_bad_flag)  .* NaN;
         
         kaz_ge_single_Ppeak(kaz_ge_single_Ppeak   == kaz_ge_single_bad_flag) = kaz_ge_single_Ppeak(kaz_ge_single_Ppeak == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vkurt(kaz_ge_single_Vkurt   == kaz_ge_single_bad_flag) = kaz_ge_single_Vkurt(kaz_ge_single_Vkurt == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vleft(kaz_ge_single_Vleft   == kaz_ge_single_bad_flag) = kaz_ge_single_Vleft(kaz_ge_single_Vleft == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vpeak(kaz_ge_single_Vpeak   == kaz_ge_single_bad_flag) = kaz_ge_single_Vpeak(kaz_ge_single_Vpeak == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vright(kaz_ge_single_Vright == kaz_ge_single_bad_flag) = kaz_ge_single_Vright(kaz_ge_single_Vright == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vslope_lt(kaz_ge_single_Vslope_lt == kaz_ge_single_bad_flag) = kaz_ge_single_Vslope_lt(kaz_ge_single_Vslope_lt == kaz_ge_single_bad_flag)  .* NaN;
         kaz_ge_single_Vslope_rt(kaz_ge_single_Vslope_rt == kaz_ge_single_bad_flag) = kaz_ge_single_Vslope_rt(kaz_ge_single_Vslope_rt == kaz_ge_single_bad_flag)  .* NaN;
         
         %% calculate some time vectors
         
         %kaz_time_stamp(:,8)  = kaz_time_stamp(:,4) + kaz_time_stamp(:,5)./60 + kaz_time_stamp(:,6)./(60*60);
         kaz_ge_single_time_hr   = kaz_ge_single_hour + kaz_ge_single_minute./60 + ...
            kaz_ge_single_second./(60*60) + (kaz_ge_single_millisecond./1000)./(60*60);
         
         year_str          = num2str(kaz_ge_single_year(1,1));
         month_num_str     = num2str_2digits(kaz_ge_single_month(1,1));
         dayofmonth_str    = num2str_2digits(kaz_ge_single_dayofmonth(1,1));
         
         %%    6. Plot daily crw Z, Vmean and Skewness
         
         % calculate the dht
         plot_ht                       = kaz_ge_single_range ./1000;
         plot_dht                      = plot_ht(6) - plot_ht(5);
         
         % rename the time variable
         %plot_time                     = kaz_spc_ge_time_hr;
         
         % define some plot attributes
         max_ht      = 8;
         
         % make some strings
         %year_str       = num2str_2digits(kaz_spc_ge_time(1,1));
         %month_str      = num2str_2digits(kaz_spc_ge_time(1,2));
         %day_str        = num2str_2digits(kaz_spc_ge_time(1,3));
         
         start_hour  = 0;
         end_hour    = 24;
         
         disp(['processing daily images, max ht: ',num2str(max_ht),' km...']);
         
         f  = (kaz_ge_single_time_hr >= start_hour) & (kaz_ge_single_time_hr <= end_hour);
         if(sum(f) > 10)
            plot_zdb             = kaz_ge_single_zdb(f,:);
            plot_Vmean           = kaz_ge_single_Vmean(f,:);
            plot_Vskew           = kaz_ge_single_Vskew(f,:);
            plot_time            = kaz_ge_single_time_hr(f);
            
            % if there are time-gaps in profiles, the plot will have
            % streaks. Put a profile of NaNs between any profiles that are
            % too far apart
            [plot_time, plot_zdb, plot_Vmean, plot_Vskew] = ...
               func_fill_time_gaps_with_NaN_profiles(plot_time,...
               plot_zdb, plot_Vmean, plot_Vskew);
            
            figure
            colormap('jet(20)')
            
            % Top Panel
            subplot(3,1,1)
            pcolor(plot_time,(plot_ht-plot_dht/2),plot_zdb');
            hold on
            shading flat
            caxis([-50 30])
            colorbar('ticks',-50:10:30,'TickLabels',{'-50','-40','-30','-20','-10','0','10','20','30','40'})
            axis([start_hour end_hour 0 max_ht]);
            grid on
            %set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
            set(gca,'xtick',start_hour:1:end_hour,'xticklabel',{'0',' ',' ','3',' ',' ','6',' ',' ','9',' ',' ','12',' ',' ','15',' ',' ','18',' ',' ','21',' ',' ','24'});
            set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
            ylabel('Height [km]')
            title(['a. ',site_prefix,', GE, KAZR, 15sec, Ref. [dBZ], ',year_str,'-',month_num_str,'-',dayofmonth_str])
            
            % Middle Panel
            subplot(3,1,2)
            pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vmean');
            hold on
            shading flat
            %caxis([-1 2])
            %colorbar('ticks',-1:0.5:2,'TickLabels',{'-1',' ','0',' ','1',' ','2'})
            caxis([-2 max_vel])
            colorbar('ticks',-2:1:max_vel,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
            axis([start_hour end_hour 0 max_ht]);
            grid on
            %set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
            set(gca,'xtick',start_hour:1:end_hour,'xticklabel',{'0',' ',' ','3',' ',' ','6',' ',' ','9',' ',' ','12',' ',' ','15',' ',' ','18',' ',' ','21',' ',' ','24'});
            set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
            ylabel('Height [km]')
            title(['b. ',site_prefix,', GE, KAZR, 15sec, Mean Velocity [m/s] (+ = down)'])
            
            % Bottom Panel
            subplot(3,1,3)
            pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vskew');
            hold on
            shading flat
            caxis([-1 1])
            colorbar('ticks',-1:0.5:1,'TickLabels',{'-1.0','-0.5','0.0','0.5','1.0'})
            axis([start_hour end_hour 0 max_ht])
            grid on
            %set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
            set(gca,'xtick',start_hour:1:end_hour,'xticklabel',{'0',' ',' ','3',' ',' ','6',' ',' ','9',' ',' ','12',' ',' ','15',' ',' ','18',' ',' ','21',' ',' ','24'});
            set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
            ylabel('Height [km]')
            xlabel('Hour of day, UTC');
            title(['c. ',site_prefix,', GE, KAZR, 15sec, Velocity Skewness'])
            
            %filename = ['fig_kazr_ge_raw_ZVmeanVsig_2016_',month_str,day_str,'_',num2str_2digits(start_hour),'.tif'];
            filename = [output_daily_dir,'fig_day_',site_prefix,'kazrge_15sec_crw_ZVVskew_',year_str,month_num_str,dayofmonth_str,'.tif'];
            print('-dtiff',filename)
            close all
            
         end % end if(sum(f) > 10)
         
         
         %%    Plot daily crw Z, Vmean and Spectrum Width
         
         % calculate the dht
         plot_ht                       = kaz_ge_single_range ./1000;
         plot_dht                      = plot_ht(6) - plot_ht(5);
         
         % rename the time variable
         %plot_time                     = kaz_spc_ge_time_hr;
         
         % define some plot attributes
         max_ht      = 8;
         
         % make some strings
         %year_str       = num2str_2digits(kaz_spc_ge_time(1,1));
         %month_str      = num2str_2digits(kaz_spc_ge_time(1,2));
         %day_str        = num2str_2digits(kaz_spc_ge_time(1,3));
         
         start_hour  = 0;
         end_hour    = 24;
         
         %disp(['processing daily image, max ht: ',num2str(max_ht),' km...']);
         
         f  = (kaz_ge_single_time_hr >= start_hour) & (kaz_ge_single_time_hr <= end_hour);
         if(sum(f) > 10)
            plot_snr             = kaz_ge_single_snr(f,:);
            plot_zdb             = kaz_ge_single_zdb(f,:);
            plot_Vmean           = kaz_ge_single_Vmean(f,:);
            plot_Vsig            = kaz_ge_single_Vsig(f,:);
            %plot_snr             = kaz_ge_single_snr(f,:);
            plot_time            = kaz_ge_single_time_hr(f);
                        
            % if there are time-gaps in profiles, the plot will have
            % streaks. Put a profile of NaNs between any profiles that are
            % too far apart
            [plot_time, plot_zdb, plot_Vmean, plot_Vsig] = ...
               func_fill_time_gaps_with_NaN_profiles(plot_time,...
               plot_zdb, plot_Vmean, plot_Vsig);
            
            figure
            colormap('jet(20)')
            
            % Top Panel
            subplot(3,1,1)
            pcolor(plot_time,(plot_ht-plot_dht/2),plot_zdb');
            hold on
            shading flat
            caxis([-50 30])
            colorbar('ticks',-50:10:30,'TickLabels',{'-50','-40','-30','-20','-10','0','10','20','30','40'})
            axis([start_hour end_hour 0 max_ht]);
            grid on
            %set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
            set(gca,'xtick',start_hour:1:end_hour,'xticklabel',{'0',' ',' ','3',' ',' ','6',' ',' ','9',' ',' ','12',' ',' ','15',' ',' ','18',' ',' ','21',' ',' ','24'});
            set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
            ylabel('Height [km]')
            title(['a. ',site_prefix,', GE, KAZR, 15sec, Ref. [dBZ], ',year_str,'-',month_num_str,'-',dayofmonth_str])
            
            % Middle Panel
            subplot(3,1,2)
            pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vmean');
            hold on
            shading flat
            %caxis([-1 2])
            %colorbar('ticks',-1:0.5:2,'TickLabels',{'-1',' ','0',' ','1',' ','2'})
            caxis([-2 max_vel])
            colorbar('ticks',-2:1:max_vel,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
            %caxis([-2 8])
            %colorbar('ticks',-2:1:8,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
            axis([start_hour end_hour 0 max_ht]);
            grid on
            %set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
            set(gca,'xtick',start_hour:1:end_hour,'xticklabel',{'0',' ',' ','3',' ',' ','6',' ',' ','9',' ',' ','12',' ',' ','15',' ',' ','18',' ',' ','21',' ',' ','24'});
            set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
            ylabel('Height [km]')
            title(['b. ',site_prefix,', GE, KAZR, 15sec, Mean Velocity [m/s] (+ = down)'])
            
            % Bottom Panel
            subplot(3,1,3)
            pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vsig');
            hold on
            shading flat
            caxis([0 max_Vsig])
            colorbar('ticks',0:0.25:max_Vsig,'TickLabels',{'0.0',' ','0.5',' ','1.0',' ','1.5',' ','2.0'})
            axis([start_hour end_hour 0 max_ht])
            grid on
            %set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
            set(gca,'xtick',start_hour:1:end_hour,'xticklabel',{'0',' ',' ','3',' ',' ','6',' ',' ','9',' ',' ','12',' ',' ','15',' ',' ','18',' ',' ','21',' ',' ','24'});
            set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
            ylabel('Height [km]')
            xlabel('Hour of day, UTC');
            title(['c. ',site_prefix,', GE, KAZR, 15sec, Spectrum STD [m/s]'])
            
            %filename = ['fig_kazr_ge_raw_ZVmeanVsig_2016_',month_str,day_str,'_',num2str_2digits(start_hour),'.tif'];
            filename = [output_daily_dir,'fig_day_',site_prefix,'kazrge_15sec_crw_ZVVsig_',year_str,month_num_str,dayofmonth_str,'.tif'];
            print('-dtiff',filename)
            close all
            
         end % end if(sum(f) > 10)
         
         %%    6. Plot hourly crw Z, Vmean and Spectrum Width
         
         if(hourly_flag)
            
            % calculate the dht
            plot_ht                       = kaz_ge_single_range ./1000;
            plot_dht                      = plot_ht(6) - plot_ht(5);
            
            % rename the time variable
            %plot_time                     = kaz_spc_ge_time_hr;
            
            % define some plot attributes
            max_ht      = 8;
            
            % make some strings
            %year_str       = num2str_2digits(kaz_spc_ge_time(1,1));
            %month_str      = num2str_2digits(kaz_spc_ge_time(1,2));
            %day_str        = num2str_2digits(kaz_spc_ge_time(1,3));
            
            disp(['processing hourly images, max ht: ',num2str(max_ht),' km...']);
            %disp(['processing day: ',year_str,'-',month_num_str,'-',dayofmonth_str]);
            
            % process each 1-hr section
            for r = 0:1:23
               %for r = 12:12
               start_hour  = r;
               end_hour    = r + 1;
               
               %disp(['processing hour: ',num2str(start_hour),', kazr single moment, max ht: ',num2str(max_ht),' km...']);
               
               f  = (kaz_ge_single_time_hr >= start_hour) & (kaz_ge_single_time_hr <= end_hour);
               if(sum(f) > 10)
                  plot_zdb             = kaz_ge_single_zdb(f,:);
                  plot_Vmean           = kaz_ge_single_Vmean(f,:);
                  plot_Vskew           = kaz_ge_single_Vskew(f,:);
                  %plot_snr             = kaz_ge_single_snr(f,:);
                  plot_time            = kaz_ge_single_time_hr(f);
                  
                  % if there are time-gaps in profiles, the plot will have
                  % streaks. Put a profile of NaNs between any profiles that are
                  % too far apart
                  
                  % copy the original time for the snr gap routine
                  plot_time_org  = plot_time;
                  [plot_time, plot_zdb, plot_Vmean, plot_Vskew] = ...
                     func_fill_time_gaps_with_NaN_profiles(plot_time,...
                     plot_zdb, plot_Vmean, plot_Vskew);
                  
                  figure
                  colormap('jet(20)')
                  
                  % Top Panel
                  subplot(3,1,1)
                  pcolor(plot_time,(plot_ht-plot_dht/2),plot_zdb');
                  hold on
                  shading flat
                  caxis([-50 30])
                  colorbar('ticks',-50:10:30,'TickLabels',{'-50','-40','-30','-20','-10','0','10','20','30','40'})
                  axis([start_hour end_hour 0 max_ht]);
                  grid on
                  set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
                  set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
                  ylabel('Height [km]')
                  title(['a. ',site_prefix,', GE, KAZR, Ref. [dBZ], ',year_str,'-',month_num_str,'-',dayofmonth_str,', Hour: ',num2str_2digits(start_hour)])
                  
                  % Middle Panel
                  subplot(3,1,2)
                  pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vmean');
                  hold on
                  shading flat
                  %caxis([-1 2])
                  %colorbar('ticks',-1:0.5:2,'TickLabels',{'-1',' ','0',' ','1',' ','2'})
                  caxis([-2 max_vel])
                  colorbar('ticks',-2:1:max_vel,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
                  %caxis([-2 8])
                  %colorbar('ticks',-2:1:8,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
                  axis([start_hour end_hour 0 max_ht]);
                  grid on
                  set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
                  set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
                  ylabel('Height [km]')
                  title(['b. ',site_prefix,', GE, KAZR, Mean Velocity [m/s] (+ = down)'])
                  
                  % Bottom Panel
                  subplot(3,1,3)
                  pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vskew');
                  hold on
                  shading flat
                  caxis([-1 1])
                  colorbar('ticks',-1:0.5:1,'TickLabels',{'-1.0','-0.5','0.0','0.5','1.0'})
                  axis([start_hour end_hour 0 max_ht])
                  grid on
                  set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
                  set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
                  ylabel('Height [km]')
                  xlabel(['Minute of hour ',num2str_2digits(start_hour),' UTC']);
                  title(['c. ',site_prefix,', GE, KAZR, Velocity Skewness'])
                  
                  %filename = ['fig_kazr_ge_raw_ZVmeanVsig_2016_',month_str,day_str,'_',num2str_2digits(start_hour),'.tif'];
                  filename = [output_hourly_dir,'fig_hr_',site_prefix,'kazrge_15sec_crw_ZVVskew_',year_str,month_num_str,dayofmonth_str,num2str_2digits(start_hour),'.tif'];
                  print('-dtiff',filename)
                  close all
                  
               end % end if(sum(f) > 10)
               
            end % end for r loop
            
         end % end if(hourly_flag)
         
         %%    Plot hourly crw Z, Vmean and Spectrum Width - decluttered
         
         if(hourly_flag)
            
            % calculate the dht
            plot_ht                       = kaz_ge_single_range ./1000;
            plot_dht                      = plot_ht(6) - plot_ht(5);
            
            % rename the time variable
            %plot_time                     = kaz_spc_ge_time_hr;
            
            % define some plot attributes
            max_ht      = 8;
            
            % make some strings
            %year_str       = num2str_2digits(kaz_spc_ge_time(1,1));
            %month_str      = num2str_2digits(kaz_spc_ge_time(1,2));
            %day_str        = num2str_2digits(kaz_spc_ge_time(1,3));
            
            %disp(['processing day: ',year_str,'-',month_num_str,'-',dayofmonth_str,' clutter flag, hourly...']);
            
            % process each 1-hr section
            for r = 0:1:23
               %for r = 2:2
               start_hour  = r;
               end_hour    = r + 1;
               
               %disp(['processing hour: ',num2str(start_hour),', kazr single moment, max ht: ',num2str(max_ht),' km...']);
               
               f  = (kaz_ge_single_time_hr >= start_hour) & (kaz_ge_single_time_hr <= end_hour);
               if(sum(f) > 10)
                  plot_snr             = kaz_ge_single_snr(f,:);
                  plot_zdb             = kaz_ge_single_zdb(f,:);
                  plot_Vmean           = kaz_ge_single_Vmean(f,:);
                  plot_Vsig            = kaz_ge_single_Vsig(f,:);
                  %plot_snr             = kaz_ge_single_snr(f,:);
                  plot_time            = kaz_ge_single_time_hr(f);
                  
                  % if there are time-gaps in profiles, the plot will have
                  % streaks. Put a profile of NaNs between any profiles that are
                  % too far apart
                  
                  % copy the original time for the snr gap routine
                  plot_time_org  = plot_time;
                  [plot_time, plot_zdb, plot_Vmean, plot_Vsig] = ...
                     func_fill_time_gaps_with_NaN_profiles(plot_time,...
                     plot_zdb, plot_Vmean, plot_Vsig);
                  
                  figure
                  colormap('jet(20)')
                  
                  % Top Panel
                  subplot(3,1,1)
                  pcolor(plot_time,(plot_ht-plot_dht/2),plot_zdb');
                  hold on
                  shading flat
                  caxis([-50 30])
                  colorbar('ticks',-50:10:30,'TickLabels',{'-50','-40','-30','-20','-10','0','10','20','30','40'})
                  axis([start_hour end_hour 0 max_ht]);
                  grid on
                  set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
                  set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
                  ylabel('Height [km]')
                  title(['a. ',site_prefix,', GE, KAZR, 15sec, Ref. [dBZ], ',year_str,'-',month_num_str,'-',dayofmonth_str,', Hour: ',num2str_2digits(start_hour)])
                  
                  % Middle Panel
                  subplot(3,1,2)
                  pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vmean');
                  hold on
                  shading flat
                  %caxis([-1 2])
                  %colorbar('ticks',-1:0.5:2,'TickLabels',{'-1',' ','0',' ','1',' ','2'})
                  caxis([-2 max_vel])
                  colorbar('ticks',-2:1:max_vel,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
                  %caxis([-1 1])
                  %colorbar('ticks',-1:0.25:1,'TickLabels',{'-1.0',' ','-0.5',' ','0.0',' ','0.5',' ','1.0'})
                  %colorbar('ticks',-6:2:6,'TickLabels',{'-6','-4','-2','0','2','4','6'})
                  %colorbar('ticks',-8:1:6,'TickLabels',{'-8',' ','-6',' ','-4',' ','-2',' ','0',' ','2'})
                  %colorbar('ticks',-2:1:8,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
                  axis([start_hour end_hour 0 max_ht]);
                  grid on
                  set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
                  set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
                  ylabel('Height [km]')
                  title(['b. ',site_prefix,', GE, KAZR, 15sec, Mean Velocity [m/s] (+ = down)'])
                  
                  % Bottom Panel
                  subplot(3,1,3)
                  pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vsig');
                  hold on
                  shading flat
                  caxis([0 max_Vsig])
                  colorbar('ticks',0:0.25:max_Vsig,'TickLabels',{'0.0',' ','0.5',' ','1.0',' ','1.5',' ','2.0'})
                  axis([start_hour end_hour 0 max_ht])
                  grid on
                  set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
                  set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
                  ylabel('Height [km]')
                  xlabel(['Minute of hour ',num2str_2digits(start_hour),' UTC']);
                  title(['c. ',site_prefix,', GE, KAZR, 15sec, Spectrum STD [m/s]'])
                  
                  %filename = ['fig_kazr_ge_raw_ZVmeanVsig_2016_',month_str,day_str,'_',num2str_2digits(start_hour),'.tif'];
                  filename = [output_hourly_dir,'fig_hr_',site_prefix,'kazrge_15sec_crw_ZVVsig_',year_str,month_num_str,dayofmonth_str,num2str_2digits(start_hour),'.tif'];
                  print('-dtiff',filename)
                  close all
                  
               end % end if(sum(f) > 10)
               
            end % end for r loop
         end % end if(hourly_flag)
          
      end % end if(good_input_spc_filename == 2)
      
   end % end for day_loop
end % end for month_loop

%% turn off the parallel cores

%matlabpool close
