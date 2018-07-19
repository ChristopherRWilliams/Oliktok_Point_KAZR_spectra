function [output_nos_mean, output_nos_std, output_nos_thres, output_spc_lin] = ...
   func_replace_Nstd_nos_with_NaN(input_spc_lin, input_Nspc, input_Nstd)


% This routine comresses spectra based on noise statistics.

% updated: 27-September-2017
% ========================================================================

%%    Define the variable names 

% Input variables:
% input_spc_lin            - 3D matrix of spectra
% input_Nspc               - 2D matrix of number of spectra in average
% input_Nstd               - number of std to set to NaN

% Outputs: 
% output_nos_mean          - noise mean
% output_nos_std           - noise std
% output_nos_thres         - noise threshold
% output_spc_lin           - 3D matrix of spectra with NaNs

%% Predefine output variables

[m,n,p]                    = size(input_spc_lin);

output_nos_mean            = zeros(m,n);
output_nos_std             = zeros(m,n);
output_nos_thres           = zeros(m,n);

output_spc_lin             = input_spc_lin;

initial_seed_length        = p/8;

%% Process each profile

for r = 1:m
   
   %% display a message
   if(rem(r,200) == 0)
      disp(['...putting NaNs into spc...processing profile #',num2str(r),' out of ',num2str(m),' profiles...']);
   end % end if(rem(prof_num,60) == 0)

   %% process each range gate
   
   for c = 1:n
      
      %% Get the spectrum to process
      
      spc   = squeeze(input_spc_lin(r,c,:));
      Nspc  = input_Nspc(r,c);      

      % process this spc      
      f = ~isnan(spc);
      if(sum(f) > 0)
                  
         % Perform the HS noise search
         [nos_mean, nos_max] = func_find_mean_HS_noise(spc, Nspc, initial_seed_length);
         
         % find the noise standard deviation
         g  = spc < nos_max;
         if(sum(g) > 3)
            nos_std  = std(spc(g));
         else
            nos_std  = NaN;
         end % end if(sum(f) > 3)
         
         % save the noise stats
         output_nos_mean(r,c)    = nos_mean;
         output_nos_std(r,c)     = nos_std;
         output_nos_thres(r,c)   = nos_max;

         % Determine # of std threshold
         test_ratio = floor((nos_max - nos_mean) / nos_std);         
         if(test_ratio < input_Nstd)            
            ratio    = test_ratio;
         else
            ratio    = input_Nstd;            
         end % end if(test_ratio < input_Nstd)
         
         % Set the output spc to NaN
         g  = spc < (nos_mean + ratio*nos_std);
         if(sum(g) > 0)            
            output_spc_lin(r,c,g) = ones(sum(g),1) .* NaN;
         end % end if(sum(g) > 0)
         
      end % end if(sum(f) > 0)
   end % end for c loop
end % end for r loop
