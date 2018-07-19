function [spk_single, index_limit] = func_find_single_peak(spk)

% Find a single peak containing the largest magnitude spectrum value.

% Input Variables
% spk             - input spectrum, bad values contain NaN

% Output Variables
% spk_clean     the same as input spk except outside points are NaN'd
% index_limit   End index values for the found peak
%               index_limit(1) -first velocity bin used in peak
%               index_limit(2) -last veocity bin used in peak
% Note: spk(index_limit(1):index_limit(2)) are all valid spectrum values 

% This routine finds the single peak by finding the maximum value and
%   marching down both sides of the peak until a NaN is found.

% Steps:
% 1. Check to see that there are valid points in spectra
% 2. Find the maximum value in the spectra.
% 3. Move to the left of the max value
% 4. Move to the left of the max value
% 5. NaN out all values outside of the indices

% Updated: 27-June-2018
% ************************************************************************

%% 1. Check to see that there are valid points in spectra
f = ~isnan(spk);

if(sum(f) < 2)
    spk_single = spk;
    index_limit = [NaN NaN];
    return
end % end if(sum(f) < 2)

%% 2. Find the maximum value in the spectra.

[~, i_max] = max(spk);

%% 3. Move to the left of the max value

current_index = i_max;
next_index    = current_index;
get_next      = 1;

while(get_next)
    % decrement the index
    next_index = next_index - 1;
    
    if(next_index < 1)
        next_index    = 1;
        get_next      = 0;
    end % end if(next_index < 1)
    
    % Get the value of the next spectral point
    next_value = spk(next_index);
    
    % Is next_value below the noise threshold?
    if(isnan(next_value))
        get_next = 0;
        % current_index is back one value...
    
    else  % the next value is value, is it decreasing or increasing?

       % move the current_index       
       current_index = next_index;
    
    end % end if(isnan(next_value))
    
end % end while(get_next)

% Place the current_index into the output variable
index_limit(1) = current_index;

%% 4. Move to the left of the max value

current_index = i_max;
next_index    = current_index;
get_next      = 1;

while(get_next)
    % increment the index
    next_index = next_index + 1;
    
    if(next_index > length(spk))
        next_index    = length(spk);
        get_next      = 0;
    end % end if(next_index > length(spk))
    
    % Get the value of the next spectral point
    next_value = spk(next_index);
    
    % Is next_value below the noise threshold?
    if(isnan(next_value))
        get_next = 0;
        % current_index is back one value
        
    else
       
       % move the current_index 
       current_index = next_index;
    
    end % end if(isnan(next_value))

end % end while(get_next)

% Place the current_index into the output variable
index_limit(2) = current_index;


%% 5. NaN out all values outside of the indices
% *********************************************
spk_single = spk;

% NaN out values to the left
if(index_limit(1) > 1)
   end_index = index_limit(1)-1;

   spk_single(1:end_index)         = spk_single(1:end_index) .* NaN;
end % end if(index_limit(1) > 1)

% look to the right of the peak
if(index_limit(2) < length(spk))
   start_index = index_limit(2)+1;
   
   spk_single(start_index:length(spk_single)) = ...
      spk_single(start_index:length(spk_single)) .* NaN;
end % end if(index_limit(1) > 1)

