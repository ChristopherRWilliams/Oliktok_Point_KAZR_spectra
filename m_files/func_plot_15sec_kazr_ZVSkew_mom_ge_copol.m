function [done] = func_plot_15sec_kazr_ZVSkew_mom_ge_copol(input_filename, output_directory, output_filename_root, site_prefix, max_ht, hourly_flag)

% This routine generates plots of ZVW using 15sec observations.

% updated = '29-July-2017'
% ========================================================================

good_input_filename  = exist(input_filename,'file');

%%    5. Load in this daily file

if(good_input_filename == 2)
   
   % load moments
   %disp(['loading netCDF file: ',input_netCDF_filename,' ...']);
   disp(['loading mat file: ',input_filename,' ...']);
   
   % the input file is named: input_netCDF_filename
   %func_load_all_netCDF_variables
   load(input_filename)
   % the variables are:
   %   kaz_ge_range                       608x1                  4864  double
   %   kaz_ge_single_Ppeak               5760x608            28016640  double
   %   kaz_ge_single_Vkurt               5760x608            28016640  double
   %   kaz_ge_single_Vleft               5760x608            28016640  double
   %   kaz_ge_single_Vmean               5760x608            28016640  double
   %   kaz_ge_single_Vpeak               5760x608            28016640  double
   %   kaz_ge_single_Vright              5760x608            28016640  double
   %   kaz_ge_single_Vsig                5760x608            28016640  double
   %   kaz_ge_single_Vskew               5760x608            28016640  double
   %   kaz_ge_single_Vslope_lt           5760x608            28016640  double
   %   kaz_ge_single_Vslope_rt           5760x608            28016640  double
   %   kaz_ge_single_index_lt            5760x608            28016640  double
   %   kaz_ge_single_index_rt            5760x608            28016640  double
   %   kaz_ge_single_max_1bin_delta      5760x608            28016640  double
   %   kaz_ge_single_max_2bin_delta      5760x608            28016640  double
   %   kaz_ge_single_nos                 5760x608            28016640  double
   %   kaz_ge_single_nos_adj             5760x608            28016640  double
   %   kaz_ge_single_nos_max             5760x608            28016640  double
   %   kaz_ge_single_nos_mean            5760x608            28016640  double
   %   kaz_ge_single_nos_profile         5760x1                 46080  double
   %   kaz_ge_single_nos_std             5760x608            28016640  double
   %   kaz_ge_single_snr                 5760x608            28016640  double
   %   kaz_ge_single_snr_adj             5760x608            28016640  double
   %   kaz_ge_single_std_1bin            5760x608            28016640  double
   %   kaz_ge_single_std_2bin            5760x608            28016640  double
   %   kaz_ge_single_zdb_arm             5760x608            28016640  double
   %   kaz_ge_single_zdb_arm_near        5760x608            28016640  double
   %   kaz_ge_single_zdb_crw             5760x608            28016640  double
   %   kaz_ge_time_stamp                 5760x7                322560  double
   
   %% Set all bad values to
   
   %          kaz_ge_single_zdb(kaz_ge_single_zdb       == kaz_ge_single_bad_flag) = kaz_ge_single_zdb(kaz_ge_single_zdb == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vmean(kaz_ge_single_Vmean   == kaz_ge_single_bad_flag) = kaz_ge_single_Vmean(kaz_ge_single_Vmean == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vsig(kaz_ge_single_Vsig     == kaz_ge_single_bad_flag) = kaz_ge_single_Vsig(kaz_ge_single_Vsig == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_snr(kaz_ge_single_snr       == kaz_ge_single_bad_flag) = kaz_ge_single_snr(kaz_ge_single_snr == kaz_ge_single_bad_flag)  .* NaN;
   %
   %          kaz_ge_single_Ppeak(kaz_ge_single_Ppeak   == kaz_ge_single_bad_flag) = kaz_ge_single_Ppeak(kaz_ge_single_Ppeak == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vkurt(kaz_ge_single_Vkurt   == kaz_ge_single_bad_flag) = kaz_ge_single_Vkurt(kaz_ge_single_Vkurt == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vleft(kaz_ge_single_Vleft   == kaz_ge_single_bad_flag) = kaz_ge_single_Vleft(kaz_ge_single_Vleft == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vpeak(kaz_ge_single_Vpeak   == kaz_ge_single_bad_flag) = kaz_ge_single_Vpeak(kaz_ge_single_Vpeak == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vright(kaz_ge_single_Vright == kaz_ge_single_bad_flag) = kaz_ge_single_Vright(kaz_ge_single_Vright == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vslope_lt(kaz_ge_single_Vslope_lt == kaz_ge_single_bad_flag) = kaz_ge_single_Vslope_lt(kaz_ge_single_Vslope_lt == kaz_ge_single_bad_flag)  .* NaN;
   %          kaz_ge_single_Vslope_rt(kaz_ge_single_Vslope_rt == kaz_ge_single_bad_flag) = kaz_ge_single_Vslope_rt(kaz_ge_single_Vslope_rt == kaz_ge_single_bad_flag)  .* NaN;
   
   %% calculate some time vectors
   
   % This is for the netCDF file format
   %kaz_time_stamp(:,8)  = kaz_time_stamp(:,4) + kaz_time_stamp(:,5)./60 + kaz_time_stamp(:,6)./(60*60);
   %          kaz_ge_single_time_hr   = kaz_ge_single_hour + kaz_ge_single_minute./60 + ...
   %             kaz_ge_single_second./(60*60) + (kaz_ge_single_millisecond./1000)./(60*60);
   %
   %          year_str          = num2str(kaz_ge_single_year(1,1));
   %          month_num_str     = num2str_2digits(kaz_ge_single_month(1,1));
   %          dayofmonth_str    = num2str_2digits(kaz_ge_single_dayofmonth(1,1));
   
   % This is for the netCDF file format
   kaz_ge_single_time_hr   = kaz_ge_time_stamp(:,4) + kaz_ge_time_stamp(:,5)./60 + kaz_ge_time_stamp(:,6)./(60*60);
   
   year_str          = num2str(kaz_ge_time_stamp(1,1));
   month_num_str     = num2str_2digits(kaz_ge_time_stamp(1,2));
   dayofmonth_str    = num2str_2digits(kaz_ge_time_stamp(1,3));
   
   %%    6. Plot hourly crw Z, Vmean and Spectrum Width
   
   if(hourly_flag)
      
      % calculate the dht
      plot_ht                       = kaz_ge_range ./1000;
      plot_dht                      = plot_ht(6) - plot_ht(5);
      
      % rename the time variable
      %plot_time                     = kaz_spc_ge_time_hr;
      
      % define some plot attributes
      %max_ht      = 6;
      
      % make some strings
      %year_str       = num2str_2digits(kaz_spc_ge_time(1,1));
      %month_str      = num2str_2digits(kaz_spc_ge_time(1,2));
      %day_str        = num2str_2digits(kaz_spc_ge_time(1,3));
      
      %disp(['processing day: ',year_str,'-',month_num_str,'-',dayofmonth_str]);
      
      % process each 1-hr section
      for r = 0:1:23
         %for r = 12:12
         start_hour  = r;
         end_hour    = r + 1;
         
         disp(['processing hour: ',num2str(start_hour),', kazr single moment, max ht: ',num2str(max_ht),' km...']);
         
         f  = (kaz_ge_single_time_hr >= start_hour) & (kaz_ge_single_time_hr <= end_hour);
         if(sum(f) > 10)
            plot_zdb             = kaz_ge_single_zdb_crw(f,:);
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
            title(['a. ',site_prefix,', GE, KAZR, CRW, Ref. [dBZ], ',year_str,'-',month_num_str,'-',dayofmonth_str,', Hour: ',num2str_2digits(start_hour)])
            
            % Middle Panel
            subplot(3,1,2)
            pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vmean');
            hold on
            shading flat
            caxis([-2 8])
            %colorbar('ticks',-6:2:6,'TickLabels',{'-6','-4','-2','0','2','4','6'})
            %colorbar('ticks',-8:1:6,'TickLabels',{'-8',' ','-6',' ','-4',' ','-2',' ','0',' ','2'})
            colorbar('ticks',-2:1:8,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
            axis([start_hour end_hour 0 max_ht]);
            grid on
            set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
            set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
            ylabel('Height [km]')
            title(['b. ',site_prefix,', GE, KAZR, CRW, Mean Velocity [m/s] (+ = down)'])
            
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
            title(['c. ',site_prefix,', GE, KAZR, CRW, Velocity Skewness'])
            
            %filename = ['fig_kazr_ge_raw_ZVmeanVsig_2016_',month_str,day_str,'_',num2str_2digits(start_hour),'.tif'];
            filename = [output_directory,output_filename_root,year_str,'_',month_num_str,dayofmonth_str,'_',num2str_2digits(start_hour),'.tif'];
            print('-dtiff',filename)
            close all
            
         end % end if(sum(f) > 10)
         
      end % end for r loop
      
   end % end if(hourly_flag)
   
   
   %%    6. Plot daily crw Z, Vmean and Spectrum Width
   
   % calculate the dht
   plot_ht                       = kaz_ge_range ./1000;
   plot_dht                      = plot_ht(6) - plot_ht(5);
   
   % rename the time variable
   %plot_time                     = kaz_spc_ge_time_hr;
   
   % define some plot attributes
   %max_ht      = 6;
   
   % make some strings
   %year_str       = num2str_2digits(kaz_spc_ge_time(1,1));
   %month_str      = num2str_2digits(kaz_spc_ge_time(1,2));
   %day_str        = num2str_2digits(kaz_spc_ge_time(1,3));
   
   start_hour  = 0;
   end_hour    = 24;
   
   %disp(['processing day, max ht: ',num2str(max_ht),' km...']);
   
   f  = (kaz_ge_single_time_hr >= start_hour) & (kaz_ge_single_time_hr <= end_hour);
   if(sum(f) > 10)
      plot_zdb             = kaz_ge_single_zdb_crw(f,:);
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
      title(['a. ',site_prefix,', GE, KAZR, CRW, 15sec, Ref. [dBZ], ',year_str,'-',month_num_str,'-',dayofmonth_str])
      
      % Middle Panel
      subplot(3,1,2)
      pcolor(plot_time,(plot_ht-plot_dht/2),plot_Vmean');
      hold on
      shading flat
      caxis([-2 8])
      %colorbar('ticks',-6:2:6,'TickLabels',{'-6','-4','-2','0','2','4','6'})
      %colorbar('ticks',-8:1:6,'TickLabels',{'-8',' ','-6',' ','-4',' ','-2',' ','0',' ','2'})
      colorbar('ticks',-2:1:8,'TickLabels',{'-2',' ','0',' ','2',' ','4',' ','6',' ','8'})
      axis([start_hour end_hour 0 max_ht]);
      grid on
      %set(gca,'xtick',start_hour:5/60:end_hour,'xticklabel',{'0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60'});
      set(gca,'xtick',start_hour:1:end_hour,'xticklabel',{'0',' ',' ','3',' ',' ','6',' ',' ','9',' ',' ','12',' ',' ','15',' ',' ','18',' ',' ','21',' ',' ','24'});
      set(gca,'ytick',0:0.5:max_ht,'yticklabel',{'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' ','8',' ','9',' ','10'});
      ylabel('Height [km]')
      title(['b. ',site_prefix,', GE, KAZR, CRW, 15sec, Mean Velocity [m/s] (+ = down)'])
      
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
      title(['c. ',site_prefix,', GE, KAZR, CRW, 15sec, Velocity Skewness'])
      
      %filename = ['fig_kazr_ge_raw_ZVmeanVsig_2016_',month_str,day_str,'_',num2str_2digits(start_hour),'.tif'];
      filename = [output_directory,output_filename_root,'daily_',year_str,'_',month_num_str,dayofmonth_str,'.tif'];
      print('-dtiff',filename)
      close all
      
   end % end if(sum(f) > 10) 
   
end % end if(good_input_spc_filename == 2)

%% Set dummy flag

done = 1;
