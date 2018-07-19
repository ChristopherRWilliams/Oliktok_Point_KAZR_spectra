function [done] = func_make_daily_mat_files(mat_input_directory, filename_root, ...
   output_file_daily, output_file_daily_single)


wildcard_name           = [mat_input_directory,filename_root,'*.mat'];
files_to_process        = dir(wildcard_name);

% the filenames to process are accessed using:
% files_to_process(i).name

num_files_to_process    = length(files_to_process);

if(num_files_to_process > 0)
   
   % load in each file and get the matrix dimensions
   
   row_cnt  = 0;
   
   % load in these files
   for file_loop = 1:num_files_to_process
      
      filename_no_dir      = files_to_process(file_loop).name;
      input_mat_filename   = [mat_input_directory,filename_no_dir];
      
      % is this a good filename?
      good_input_mat_filename    = exist(input_mat_filename,'file');
      
      if(good_input_mat_filename == 2)
         
         % load in the time and one variable to get the dimensions
         load(input_mat_filename,'kaz_ge_time_stamp','kaz_ge_sub_snr');
         
         % Determine how many rows to keep for this file
         % =============================================
         
         % convert the time to day.fraction
         hour_of_day       = kaz_ge_time_stamp(:,4) + kaz_ge_time_stamp(:,5)./60 + kaz_ge_time_stamp(:,6)./(60*60) + (kaz_ge_time_stamp(:,7)./1000)./(60*60);
         mm                = length(hour_of_day);
         
         row_cnt           = row_cnt + mm;
         
         % Determine how many colums to keep for this file
         % =============================================
         
         [m,col_cnt,peak_cnt]    = size(kaz_ge_sub_snr);
         
      end % end if(good_input_mat_filename == 2)
   end % end for file_loop
   
   %% Make the output matrices...
   
   %disp(' ');
   %disp([' make the matrices: ',num2str(row_cnt),' vs. ',num2str(col_cnt)]);
   %disp(' ');
   
   buffer_org_time                           = ones(row_cnt,7) .* NaN;
   
   % define the dominant peak
   buffer_kaz_ge_single_Ppeak                   = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vkurt                   = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vleft                   = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vmean                   = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vmean_var               = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_cnt_Vmean_var           = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vpeak                   = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vright                  = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vsig                    = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vskew                   = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vslope_lt               = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_Vslope_rt               = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_index_lt                = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_index_rt                = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_max_1bin_delta          = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_max_2bin_delta          = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_nos                     = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_nos_adj                 = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_nos_max                 = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_nos_mean                = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_nos_profile             = ones(row_cnt,1) .* NaN;
   buffer_kaz_ge_single_nos_std                 = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_snr                     = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_snr_adj                 = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_std_1bin                = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_std_2bin                = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_zdb_arm                 = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_zdb_arm_near            = ones(row_cnt,col_cnt) .* NaN;
   buffer_kaz_ge_single_zdb_crw                 = ones(row_cnt,col_cnt) .* NaN;
   
   % define the sub-peak
   buffer_kaz_ge_sub_Ppeak                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vkurt                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vleft                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vmean                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vpeak                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vright                  = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vsig                    = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vskew                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vslope_lt               = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_Vslope_rt               = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_index_lt                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_index_rt                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_max_1bin_delta          = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_max_2bin_delta          = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_nos                     = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_snr                     = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_std_1bin                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_std_2bin                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_zdb_arm                 = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_sub_zdb_crw                 = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   
   % define the separate-peak
   buffer_kaz_ge_separate_Ppeak                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vkurt                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vleft                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vmean                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vpeak                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vright                  = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vsig                    = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vskew                   = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vslope_lt               = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_Vslope_rt               = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_index_lt                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_index_rt                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_max_1bin_delta          = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_max_2bin_delta          = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_nos                     = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_snr                     = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_std_1bin                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_std_2bin                = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_zdb_arm                 = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   buffer_kaz_ge_separate_zdb_crw                 = ones(row_cnt,col_cnt,peak_cnt) .* NaN;
   
   %% Load in the individual data files
   
   % load in these files
   end_index   = 0;
   
   for file_loop = 1:num_files_to_process
      
      filename_no_dir      = files_to_process(file_loop).name;
      input_mat_filename   = [mat_input_directory,filename_no_dir];
      
      % is this a good filename?
      good_input_mat_filename    = exist(input_mat_filename,'file');
      
      if(good_input_mat_filename == 2)
         
         if(rem(file_loop,6) == 0)
            disp(['processing file # ',num2str(file_loop),' out of ',num2str(num_files_to_process),' files...']);
            disp(['...loading file: ',input_mat_filename]);
         end % end if(rem(file_loop,50) == 0)
         
         % load this file
         load(input_mat_filename);
         
         % Determine which rows to place this data
         % =============================================
         
         % convert the time to day.fraction
         hour_of_day         = kaz_ge_time_stamp(:,4) + kaz_ge_time_stamp(:,5)./60 + kaz_ge_time_stamp(:,6)./(60*60) + (kaz_ge_time_stamp(:,7)./1000)./(60*60);
         mm                  = length(hour_of_day);
         
         if(mm > 0)
            
            start_index = end_index + 1;
            end_index   = start_index + mm - 1;
            
            [~,n_hts]   = size(kaz_ge_single_Ppeak);
            
            buffer_org_time(start_index:end_index,:)                       = kaz_ge_time_stamp(:,:);
            
            % get the dominant peak
            buffer_kaz_ge_single_Ppeak(start_index:end_index,1:n_hts)         = kaz_ge_single_Ppeak(:,:);
            buffer_kaz_ge_single_Vkurt(start_index:end_index,1:n_hts)         = kaz_ge_single_Vkurt(:,:);
            buffer_kaz_ge_single_Vleft(start_index:end_index,1:n_hts)         = kaz_ge_single_Vleft(:,:);
            buffer_kaz_ge_single_Vmean(start_index:end_index,1:n_hts)         = kaz_ge_single_Vmean(:,:);
            buffer_kaz_ge_single_Vmean_var(start_index:end_index,1:n_hts)     = kaz_ge_single_Vmean_var(:,:);
            buffer_kaz_ge_single_cnt_Vmean_var(start_index:end_index,1:n_hts) = kaz_ge_single_cnt_Vmean_var(:,:);
            buffer_kaz_ge_single_Vpeak(start_index:end_index,1:n_hts)         = kaz_ge_single_Vpeak(:,:);
            buffer_kaz_ge_single_Vright(start_index:end_index,1:n_hts)        = kaz_ge_single_Vright(:,:);
            buffer_kaz_ge_single_Vsig(start_index:end_index,1:n_hts)          = kaz_ge_single_Vsig(:,:);
            buffer_kaz_ge_single_Vskew(start_index:end_index,1:n_hts)         = kaz_ge_single_Vskew(:,:);
            buffer_kaz_ge_single_Vslope_lt(start_index:end_index,1:n_hts)     = kaz_ge_single_Vslope_lt(:,:);
            buffer_kaz_ge_single_Vslope_rt(start_index:end_index,1:n_hts)     = kaz_ge_single_Vslope_rt(:,:);
            buffer_kaz_ge_single_index_lt(start_index:end_index,1:n_hts)      = kaz_ge_single_index_lt(:,:);
            buffer_kaz_ge_single_index_rt(start_index:end_index,1:n_hts)      = kaz_ge_single_index_rt(:,:);
            buffer_kaz_ge_single_max_1bin_delta(start_index:end_index,1:n_hts) = kaz_ge_single_max_1bin_delta(:,:);
            buffer_kaz_ge_single_max_2bin_delta(start_index:end_index,1:n_hts) = kaz_ge_single_max_2bin_delta(:,:);
            buffer_kaz_ge_single_nos(start_index:end_index,1:n_hts)           = kaz_ge_single_nos(:,:);
            buffer_kaz_ge_single_nos_adj(start_index:end_index,1:n_hts)       = kaz_ge_single_nos_adj(:,:);
            buffer_kaz_ge_single_nos_max(start_index:end_index,1:n_hts)       = kaz_ge_single_nos_max(:,:);
            buffer_kaz_ge_single_nos_mean(start_index:end_index,1:n_hts)      = kaz_ge_single_nos_mean(:,:);
            buffer_kaz_ge_single_nos_profile(start_index:end_index,1)         = kaz_ge_single_nos_profile(:,1);
            buffer_kaz_ge_single_nos_std(start_index:end_index,1:n_hts)       = kaz_ge_single_nos_std(:,:);
            buffer_kaz_ge_single_snr(start_index:end_index,1:n_hts)           = kaz_ge_single_snr(:,:);
            buffer_kaz_ge_single_snr_adj(start_index:end_index,1:n_hts)       = kaz_ge_single_snr_adj(:,:);
            buffer_kaz_ge_single_std_1bin(start_index:end_index,1:n_hts)      = kaz_ge_single_std_1bin(:,:);
            buffer_kaz_ge_single_std_2bin(start_index:end_index,1:n_hts)      = kaz_ge_single_std_2bin(:,:);
            buffer_kaz_ge_single_zdb_arm(start_index:end_index,1:n_hts)       = kaz_ge_single_zdb_arm(:,:);
            buffer_kaz_ge_single_zdb_arm_near(start_index:end_index,1:n_hts)  = kaz_ge_single_zdb_arm_near(:,:);
            buffer_kaz_ge_single_zdb_crw(start_index:end_index,1:n_hts)       = kaz_ge_single_zdb_crw(:,:);
            
            % get the sub-peak
            buffer_kaz_ge_sub_Ppeak(start_index:end_index,1:n_hts,:)         = kaz_ge_sub_Ppeak(:,:,:);
            buffer_kaz_ge_sub_Vkurt(start_index:end_index,1:n_hts,:)         = kaz_ge_sub_Vkurt(:,:,:);
            buffer_kaz_ge_sub_Vleft(start_index:end_index,1:n_hts,:)         = kaz_ge_sub_Vleft(:,:,:);
            buffer_kaz_ge_sub_Vmean(start_index:end_index,1:n_hts,:)         = kaz_ge_sub_Vmean(:,:,:);
            buffer_kaz_ge_sub_Vpeak(start_index:end_index,1:n_hts,:)         = kaz_ge_sub_Vpeak(:,:,:);
            buffer_kaz_ge_sub_Vright(start_index:end_index,1:n_hts,:)        = kaz_ge_sub_Vright(:,:,:);
            buffer_kaz_ge_sub_Vsig(start_index:end_index,1:n_hts,:)          = kaz_ge_sub_Vsig(:,:,:);
            buffer_kaz_ge_sub_Vskew(start_index:end_index,1:n_hts,:)         = kaz_ge_sub_Vskew(:,:,:);
            buffer_kaz_ge_sub_Vslope_lt(start_index:end_index,1:n_hts,:)     = kaz_ge_sub_Vslope_lt(:,:,:);
            buffer_kaz_ge_sub_Vslope_rt(start_index:end_index,1:n_hts,:)     = kaz_ge_sub_Vslope_rt(:,:,:);
            buffer_kaz_ge_sub_index_lt(start_index:end_index,1:n_hts,:)      = kaz_ge_sub_index_lt(:,:,:);
            buffer_kaz_ge_sub_index_rt(start_index:end_index,1:n_hts,:)      = kaz_ge_sub_index_rt(:,:,:);
            buffer_kaz_ge_sub_max_1bin_delta(start_index:end_index,1:n_hts,:) = kaz_ge_sub_max_1bin_delta(:,:,:);
            buffer_kaz_ge_sub_max_2bin_delta(start_index:end_index,1:n_hts,:) = kaz_ge_sub_max_2bin_delta(:,:,:);
            buffer_kaz_ge_sub_nos(start_index:end_index,1:n_hts,:)           = kaz_ge_sub_nos(:,:,:);
            buffer_kaz_ge_sub_snr(start_index:end_index,1:n_hts,:)           = kaz_ge_sub_snr(:,:,:);
            buffer_kaz_ge_sub_std_1bin(start_index:end_index,1:n_hts,:)      = kaz_ge_sub_std_1bin(:,:,:);
            buffer_kaz_ge_sub_std_2bin(start_index:end_index,1:n_hts,:)      = kaz_ge_sub_std_2bin(:,:,:);
            buffer_kaz_ge_sub_zdb_arm(start_index:end_index,1:n_hts,:)       = kaz_ge_sub_zdb_arm(:,:,:);
            buffer_kaz_ge_sub_zdb_crw(start_index:end_index,1:n_hts,:)       = kaz_ge_sub_zdb_crw(:,:,:);
            
            % get the separate-peak
            buffer_kaz_ge_separate_Ppeak(start_index:end_index,1:n_hts,:)         = kaz_ge_separate_Ppeak(:,:,:);
            buffer_kaz_ge_separate_Vkurt(start_index:end_index,1:n_hts,:)         = kaz_ge_separate_Vkurt(:,:,:);
            buffer_kaz_ge_separate_Vleft(start_index:end_index,1:n_hts,:)         = kaz_ge_separate_Vleft(:,:,:);
            buffer_kaz_ge_separate_Vmean(start_index:end_index,1:n_hts,:)         = kaz_ge_separate_Vmean(:,:,:);
            buffer_kaz_ge_separate_Vpeak(start_index:end_index,1:n_hts,:)         = kaz_ge_separate_Vpeak(:,:,:);
            buffer_kaz_ge_separate_Vright(start_index:end_index,1:n_hts,:)        = kaz_ge_separate_Vright(:,:,:);
            buffer_kaz_ge_separate_Vsig(start_index:end_index,1:n_hts,:)          = kaz_ge_separate_Vsig(:,:,:);
            buffer_kaz_ge_separate_Vskew(start_index:end_index,1:n_hts,:)         = kaz_ge_separate_Vskew(:,:,:);
            buffer_kaz_ge_separate_Vslope_lt(start_index:end_index,1:n_hts,:)     = kaz_ge_separate_Vslope_lt(:,:,:);
            buffer_kaz_ge_separate_Vslope_rt(start_index:end_index,1:n_hts,:)     = kaz_ge_separate_Vslope_rt(:,:,:);
            buffer_kaz_ge_separate_index_lt(start_index:end_index,1:n_hts,:)      = kaz_ge_separate_index_lt(:,:,:);
            buffer_kaz_ge_separate_index_rt(start_index:end_index,1:n_hts,:)      = kaz_ge_separate_index_rt(:,:,:);
            buffer_kaz_ge_separate_max_1bin_delta(start_index:end_index,1:n_hts,:) = kaz_ge_separate_max_1bin_delta(:,:,:);
            buffer_kaz_ge_separate_max_2bin_delta(start_index:end_index,1:n_hts,:) = kaz_ge_separate_max_2bin_delta(:,:,:);
            buffer_kaz_ge_separate_nos(start_index:end_index,1:n_hts,:)           = kaz_ge_separate_nos(:,:,:);
            buffer_kaz_ge_separate_snr(start_index:end_index,1:n_hts,:)           = kaz_ge_separate_snr(:,:,:);
            buffer_kaz_ge_separate_std_1bin(start_index:end_index,1:n_hts,:)      = kaz_ge_separate_std_1bin(:,:,:);
            buffer_kaz_ge_separate_std_2bin(start_index:end_index,1:n_hts,:)      = kaz_ge_separate_std_2bin(:,:,:);
            buffer_kaz_ge_separate_zdb_arm(start_index:end_index,1:n_hts,:)       = kaz_ge_separate_zdb_arm(:,:,:);
            buffer_kaz_ge_separate_zdb_crw(start_index:end_index,1:n_hts,:)       = kaz_ge_separate_zdb_crw(:,:,:);
            
         end % end if(mm > 0)
         
      end % end if(good_input_mat_filename == 2)
   end % end for file_loop loop
   
   %% Rename the buffer variables to their original names
   
   %disp('paused, line 423...')
   %pause
   
   kaz_ge_time_stamp          = buffer_org_time;
   
   % get the dominant peak
   kaz_ge_single_Ppeak         = buffer_kaz_ge_single_Ppeak;
   kaz_ge_single_Vkurt         = buffer_kaz_ge_single_Vkurt;
   kaz_ge_single_Vleft         = buffer_kaz_ge_single_Vleft;
   kaz_ge_single_Vmean         = buffer_kaz_ge_single_Vmean;
   kaz_ge_single_Vmean_var     = buffer_kaz_ge_single_Vmean_var;
   kaz_ge_single_cnt_Vmean_var = buffer_kaz_ge_single_cnt_Vmean_var;
   kaz_ge_single_Vpeak         = buffer_kaz_ge_single_Vpeak;
   kaz_ge_single_Vright        = buffer_kaz_ge_single_Vright;
   kaz_ge_single_Vsig          = buffer_kaz_ge_single_Vsig;
   kaz_ge_single_Vskew         = buffer_kaz_ge_single_Vskew;
   kaz_ge_single_Vslope_lt     = buffer_kaz_ge_single_Vslope_lt;
   kaz_ge_single_Vslope_rt     = buffer_kaz_ge_single_Vslope_rt;
   kaz_ge_single_index_lt      = buffer_kaz_ge_single_index_lt;
   kaz_ge_single_index_rt      = buffer_kaz_ge_single_index_rt;
   kaz_ge_single_max_1bin_delta = buffer_kaz_ge_single_max_1bin_delta;
   kaz_ge_single_max_2bin_delta = buffer_kaz_ge_single_max_2bin_delta;
   kaz_ge_single_nos           = buffer_kaz_ge_single_nos;
   kaz_ge_single_nos_adj       = buffer_kaz_ge_single_nos_adj;
   kaz_ge_single_nos_max       = buffer_kaz_ge_single_nos_max;
   kaz_ge_single_nos_mean      = buffer_kaz_ge_single_nos_mean;
   kaz_ge_single_nos_profile   = buffer_kaz_ge_single_nos_profile;
   kaz_ge_single_nos_std       = buffer_kaz_ge_single_nos_std;
   kaz_ge_single_snr           = buffer_kaz_ge_single_snr;
   kaz_ge_single_snr_adj       = buffer_kaz_ge_single_snr_adj;
   kaz_ge_single_std_1bin      = buffer_kaz_ge_single_std_1bin;
   kaz_ge_single_std_2bin      = buffer_kaz_ge_single_std_2bin;
   kaz_ge_single_zdb_arm       = buffer_kaz_ge_single_zdb_arm;
   kaz_ge_single_zdb_arm_near  = buffer_kaz_ge_single_zdb_arm_near;
   kaz_ge_single_zdb_crw       = buffer_kaz_ge_single_zdb_crw;
   
   % get the sub-peak
   kaz_ge_sub_Ppeak         = buffer_kaz_ge_sub_Ppeak;
   kaz_ge_sub_Vkurt         = buffer_kaz_ge_sub_Vkurt;
   kaz_ge_sub_Vleft         = buffer_kaz_ge_sub_Vleft;
   kaz_ge_sub_Vmean         = buffer_kaz_ge_sub_Vmean;
   kaz_ge_sub_Vpeak         = buffer_kaz_ge_sub_Vpeak;
   kaz_ge_sub_Vright        = buffer_kaz_ge_sub_Vright;
   kaz_ge_sub_Vsig          = buffer_kaz_ge_sub_Vsig;
   kaz_ge_sub_Vskew         = buffer_kaz_ge_sub_Vskew;
   kaz_ge_sub_Vslope_lt     = buffer_kaz_ge_sub_Vslope_lt;
   kaz_ge_sub_Vslope_rt     = buffer_kaz_ge_sub_Vslope_rt;
   kaz_ge_sub_index_lt      = buffer_kaz_ge_sub_index_lt;
   kaz_ge_sub_index_rt      = buffer_kaz_ge_sub_index_rt;
   kaz_ge_sub_max_1bin_delta = buffer_kaz_ge_sub_max_1bin_delta;
   kaz_ge_sub_max_2bin_delta = buffer_kaz_ge_sub_max_2bin_delta;
   kaz_ge_sub_nos           = buffer_kaz_ge_sub_nos;
   kaz_ge_sub_snr           = buffer_kaz_ge_sub_snr;
   kaz_ge_sub_std_1bin      = buffer_kaz_ge_sub_std_1bin;
   kaz_ge_sub_std_2bin      = buffer_kaz_ge_sub_std_2bin;
   kaz_ge_sub_zdb_arm       = buffer_kaz_ge_sub_zdb_arm;
   kaz_ge_sub_zdb_crw       = buffer_kaz_ge_sub_zdb_crw;
   
   % get the separate-peak
   kaz_ge_separate_Ppeak         = buffer_kaz_ge_separate_Ppeak;
   kaz_ge_separate_Vkurt         = buffer_kaz_ge_separate_Vkurt;
   kaz_ge_separate_Vleft         = buffer_kaz_ge_separate_Vleft;
   kaz_ge_separate_Vmean         = buffer_kaz_ge_separate_Vmean;
   kaz_ge_separate_Vpeak         = buffer_kaz_ge_separate_Vpeak;
   kaz_ge_separate_Vright        = buffer_kaz_ge_separate_Vright;
   kaz_ge_separate_Vsig          = buffer_kaz_ge_separate_Vsig;
   kaz_ge_separate_Vskew         = buffer_kaz_ge_separate_Vskew;
   kaz_ge_separate_Vslope_lt     = buffer_kaz_ge_separate_Vslope_lt;
   kaz_ge_separate_Vslope_rt     = buffer_kaz_ge_separate_Vslope_rt;
   kaz_ge_separate_index_lt      = buffer_kaz_ge_separate_index_lt;
   kaz_ge_separate_index_rt      = buffer_kaz_ge_separate_index_rt;
   kaz_ge_separate_max_1bin_delta = buffer_kaz_ge_separate_max_1bin_delta;
   kaz_ge_separate_max_2bin_delta = buffer_kaz_ge_separate_max_2bin_delta;
   kaz_ge_separate_nos           = buffer_kaz_ge_separate_nos;
   kaz_ge_separate_snr           = buffer_kaz_ge_separate_snr;
   kaz_ge_separate_std_1bin      = buffer_kaz_ge_separate_std_1bin;
   kaz_ge_separate_std_2bin      = buffer_kaz_ge_separate_std_2bin;
   kaz_ge_separate_zdb_arm       = buffer_kaz_ge_separate_zdb_arm;
   kaz_ge_separate_zdb_crw       = buffer_kaz_ge_separate_zdb_crw;
   
   
   %% Save all values to a matlab file
   
   %output_file_directory   = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\mat_daily_60sec_ave_moments\';
   %output_file             = [output_file_directory,filename_root,'daily_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
   
   disp(['Saving the file: ',output_file_daily,' ...'])
   
   save(output_file_daily,'kaz_ge_time_stamp',...
      'kaz_ge_range',...
      'kaz_ge_single_Ppeak','kaz_ge_single_Vkurt','kaz_ge_single_Vleft',...
      'kaz_ge_single_Vmean','kaz_ge_single_Vpeak','kaz_ge_single_Vright',...
      'kaz_ge_single_Vmean_var','kaz_ge_single_cnt_Vmean_var',...
      'kaz_ge_single_Vsig','kaz_ge_single_Vskew','kaz_ge_single_Vslope_lt',...
      'kaz_ge_single_Vslope_rt','kaz_ge_single_index_lt','kaz_ge_single_index_rt',...
      'kaz_ge_single_max_1bin_delta','kaz_ge_single_max_2bin_delta',...
      'kaz_ge_single_nos','kaz_ge_single_nos_adj','kaz_ge_single_nos_max',...
      'kaz_ge_single_nos_mean','kaz_ge_single_nos_profile','kaz_ge_single_nos_std',...
      'kaz_ge_single_snr','kaz_ge_single_snr_adj','kaz_ge_single_zdb_arm',...
      'kaz_ge_single_std_1bin','kaz_ge_single_std_2bin',...
      'kaz_ge_single_zdb_arm_near','kaz_ge_single_zdb_crw',...
      'kaz_ge_sub_Ppeak','kaz_ge_sub_Vkurt','kaz_ge_sub_Vleft',...
      'kaz_ge_sub_Vmean','kaz_ge_sub_Vpeak','kaz_ge_sub_Vright',...
      'kaz_ge_sub_Vsig','kaz_ge_sub_Vskew','kaz_ge_sub_Vslope_lt',...
      'kaz_ge_sub_Vslope_rt','kaz_ge_sub_index_lt','kaz_ge_sub_index_rt',...
      'kaz_ge_sub_max_1bin_delta','kaz_ge_sub_max_2bin_delta',...
      'kaz_ge_sub_nos','kaz_ge_sub_snr','kaz_ge_sub_zdb_arm',...
      'kaz_ge_sub_std_1bin','kaz_ge_sub_std_2bin',...
      'kaz_ge_sub_zdb_crw',...
      'kaz_ge_separate_Ppeak','kaz_ge_separate_Vkurt','kaz_ge_separate_Vleft',...
      'kaz_ge_separate_Vmean','kaz_ge_separate_Vpeak','kaz_ge_separate_Vright',...
      'kaz_ge_separate_Vsig','kaz_ge_separate_Vskew','kaz_ge_separate_Vslope_lt',...
      'kaz_ge_separate_Vslope_rt','kaz_ge_separate_index_lt','kaz_ge_separate_index_rt',...
      'kaz_ge_separate_max_1bin_delta','kaz_ge_separate_max_2bin_delta',...
      'kaz_ge_separate_nos','kaz_ge_separate_snr','kaz_ge_separate_zdb_arm',...
      'kaz_ge_separate_std_1bin','kaz_ge_separate_std_2bin',...
      'kaz_ge_separate_zdb_crw')
   
   %% Save the single peak values to a matlab file
   
   %output_file_directory   = 'C:\Projects\Oliktok_Point\KAZR\multi_peak_temporal_ave_spc\mat_daily_60sec_ave_moments\';
   %output_file             = [output_file_directory,filename_root,'daily_single_',year_str,month_num_str,dayofmonth_str,version_suffex,'.mat'];
   
   disp(['Saving the file: ',output_file_daily_single,' ...'])
   
   save(output_file_daily_single,'kaz_ge_time_stamp',...
      'kaz_ge_range',...
      'kaz_ge_single_Ppeak','kaz_ge_single_Vkurt','kaz_ge_single_Vleft',...
      'kaz_ge_single_Vmean','kaz_ge_single_Vpeak','kaz_ge_single_Vright',...
      'kaz_ge_single_Vmean_var','kaz_ge_single_cnt_Vmean_var',...
      'kaz_ge_single_Vsig','kaz_ge_single_Vskew','kaz_ge_single_Vslope_lt',...
      'kaz_ge_single_Vslope_rt','kaz_ge_single_index_lt','kaz_ge_single_index_rt',...
      'kaz_ge_single_max_1bin_delta','kaz_ge_single_max_2bin_delta',...
      'kaz_ge_single_nos','kaz_ge_single_nos_adj','kaz_ge_single_nos_max',...
      'kaz_ge_single_nos_mean','kaz_ge_single_nos_profile','kaz_ge_single_nos_std',...
      'kaz_ge_single_snr','kaz_ge_single_snr_adj','kaz_ge_single_zdb_arm',...
      'kaz_ge_single_std_1bin','kaz_ge_single_std_2bin',...
      'kaz_ge_single_zdb_arm_near','kaz_ge_single_zdb_crw')
   
   %% Clear out the memory
   
   %clear kaz_spc* buffer*
   
end % end if(num_files_to_process > 0)

done = 1;
