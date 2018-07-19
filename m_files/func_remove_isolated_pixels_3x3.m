function [output_data] = func_remove_isolated_pixels_3x3(input_data)

% This routine removes pixels that don't have orthogonal neighbors

% Inputs:
% input_data      - 2D matrix with NaN's indicated missing data
% cnt_threshold   - number of neighbors needed to keep the center pixel

% Outputs:
% output_data     - same data as input_data, except isolated pixels are
%                    replaced with NaNs

% updated: 29-September-2017
% =======================================================================

% Define the cnt_threshold
cnt_threshold        = 3;
cnt_threshold_edge   = 2;

%% Get the size of the input matrix

[m,n]          = size(input_data);
output_data    = input_data;

%% Process each pixel that has 4 othoginal neighbors

% process each row, expect for the first and last rows
for c = 2:n-1
   
   %disp(['doing isolated pixel filter, processing row: ',num2str(c),',...']);
   
   % process each profile, expect for the first and last columns
   for r = 2:m-1
         
      % Process only pixels that are not NaN
      if(~isnan(input_data(r,c)))
         
         % Get the 8 neighbors within 3x3 
         neighbors    = [input_data(r-1,c-1) input_data(r,c-1) ...
            input_data(r+1,c-1) input_data(r-1,c) input_data(r+1,c) ...
            input_data(r-1,c+1) input_data(r,c+1) input_data(r+1,c+1)];

         % NaN center pixel if not enough neighbors
         f = ~isnan(neighbors);
         if(sum(f) < cnt_threshold)
            output_data(r,c) = NaN;
         end % end if(sum(f) < cnt_threshold)
      end % end if(~isnan(input_data(r,c)));
      
   end % end for r loop
end % end for c loop

%% Process the bottom row (3 othoginal neighbors)

% process the first row
c = 1;

%disp(['doing isolated pixel filter, processing first row: ',num2str(c),',...']);

% process each profile
for r = 2:m-1
   
   % Process only pixels that are not NaN
   if(~isnan(input_data(r,c)))
      
      % Get the 5 neighbors within 2x3 
      neighbors    = [input_data(r-1,c) input_data(r+1,c) ...
         input_data(r-1,c+1) input_data(r,c+1) input_data(r+1,c+1)];
      
      % NaN center pixel if not enough neighbors
      f = ~isnan(neighbors);
      if(sum(f) < cnt_threshold_edge)
         output_data(r,c) = NaN;
      end % end if(sum(f) < cnt_threshold_edge)
   end % end if(~isnan(input_data(r,c)));
   
end % end for r loop

%% Process the top row (3 othoginal neighbors)

% process the first row
c = n;

%disp(['doing isolated pixel filter, processing top row: ',num2str(c),',...']);

% process each profile
for r = 2:m-1
   
   % Process only pixels that are not NaN
   if(~isnan(input_data(r,c)))
      
      % Get the 5 neighbors within 2x3 
      neighbors    = [input_data(r-1,c-1) input_data(r,c-1) ...
         input_data(r+1,c-1) input_data(r-1,c) input_data(r+1,c)];
      
      % NaN center pixel if not enough neighbors
      f = ~isnan(neighbors);
      if(sum(f) < cnt_threshold_edge)
         output_data(r,c) = NaN;
      end % end if(sum(f) < cnt_threshold_edge)
   end % end if(~isnan(input_data(r,c)));
   
end % end for r loop


%% Process the first profile (3 othoginal neighbors)

% process each row, expect for the first and last rows
for c = 2:n-1
   
   %disp(['doing isolated pixel filter, processing ht: ',num2str(plot_ht(c)),' m,...']);
   
   % process each profile, expect for the first and last columns
   r = 1;
   
   % Process only pixels that are not NaN
   if(~isnan(input_data(r,c)))
      
      % Get the 5 neighbors within 2x3 
      neighbors    = [input_data(r,c-1) input_data(r+1,c-1) ...
         input_data(r+1,c) input_data(r,c+1) input_data(r+1,c+1)];
      
      % NaN center pixel if not enough neighbors
      f = ~isnan(neighbors);
      if(sum(f) < cnt_threshold_edge)
         output_data(r,c) = NaN;
      end % end if(sum(f) < cnt_threshold_edge)
   end % end if(~isnan(input_data(r,c)));
   
end % end for c loop


%% Process the last profile (3 othoginal neighbors)

% process each row, expect for the first and last rows
for c = 2:n-1
   
   %disp(['doing isolated pixel filter, processing ht: ',num2str(plot_ht(c)),' m,...']);
   
   % process each profile, expect for the first and last columns
   r = m;
   
   % Process only pixels that are not NaN
   if(~isnan(input_data(r,c)))
      
      % Get the 5 neighbors within 2x3 
      neighbors    = [input_data(r-1,c-1) input_data(r,c-1) ...
         input_data(r-1,c) input_data(r-1,c+1) input_data(r,c+1)];
      
      % NaN center pixel if not enough neighbors
      f = ~isnan(neighbors);
      if(sum(f) < cnt_threshold_edge)
         output_data(r,c) = NaN;
      end % end if(sum(f) < cnt_threshold_edge)
   end % end if(~isnan(input_data(r,c)));
   
end % end for c loop

%% Process the four corners of the matrix (2 orthoginal neighbors)

r  = 1;
c  = 1;
% location (1,1)
if(~isnan(input_data(r,c)))

   % Get the 3 neighbors within 2x2 
   neighbors    = [input_data(r+1,c) input_data(r+1,c+1) input_data(r,c+1)];
   
   % NaN center pixel if not enough neighbors
   f = ~isnan(neighbors);
   if(sum(f) < cnt_threshold_edge)
      output_data(r,c) = NaN;
   end % end if(sum(f) < cnt_threshold_edge)
end % end if(~isnan(input_data(r,c)));

r  = 1;
c  = n;
% location (1,1)
if(~isnan(input_data(r,c)))

   % Get the 3 neighbors within 2x2 
   neighbors    = [input_data(r+1,c) input_data(r+1,c-1) input_data(r,c-1)];
   
   % NaN center pixel if not enough neighbors
   f = ~isnan(neighbors);
   if(sum(f) < cnt_threshold_edge)
      output_data(r,c) = NaN;
   end % end if(sum(f) < cnt_threshold_edge)
end % end if(~isnan(input_data(r,c)));

r  = m;
c  = 1;
% location (1,1)
if(~isnan(input_data(r,c)))

   % Get the 3 neighbors within 2x2 
   neighbors    = [input_data(r-1,c) input_data(r-1,c+1) input_data(r,c+1)];
   
   % NaN center pixel if not enough neighbors
   f = ~isnan(neighbors);
   if(sum(f) < cnt_threshold_edge)
      output_data(r,c) = NaN;
   end % end if(sum(f) < cnt_threshold_edge)
end % end if(~isnan(input_data(r,c)));

r  = m;
c  = n;
% location (1,1)
if(~isnan(input_data(r,c)))

   % Get the 3 neighbors within 2x2 
   neighbors    = [input_data(r-1,c) input_data(r-1,c-1) input_data(r,c-1)];
   
   % NaN center pixel if not enough neighbors
   f = ~isnan(neighbors);
   if(sum(f) < cnt_threshold_edge)
      output_data(r,c) = NaN;
   end % end if(sum(f) < cnt_threshold_edge)
end % end if(~isnan(input_data(r,c)));

