function [spk_single, index_limit] = func_find_single_peak_Vmean_prior_valley(spk,valley_thres)

%function [spk_single, index_limit] = find_single_peak_Vmean_prior(spk, peak_index)

% Find the largest magnitude peak in spectra.

% Input Variables
% spk             - input spectrum, no weighting
% valley_thres    - depth of valley threshold

% Output Variables
% spk_clean     the same as input spk except outside points are NaN'd
% index_limit   End index values for the found peak
%               index_limit(1) -first velocity bin used in peak
%               index_limit(2) -last veocity bin used in peak

% This routine finds the single peak by finding the maximum value and
%   marching down both sides of the peak until a NaN is found.
% The new index is used to filter the spectrum

% Steps:
% 1. Check to see that there are valid points in spectra
% 2. Find the maximum value in the spectra.
% 3. Move to the left of the max value
% 4. Move to the left of the max value
% 5. NaN out all values outside of the indices

% Updated: 24-November-2014
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
current_value = spk(current_index);
next_index    = current_index;
get_next      = 1;

while(get_next)
    % decrement the index
    next_index = next_index - 1;
    
    if(next_index < 1)
        next_index    = 1;
        current_index = 1;
        get_next      = 0;
    end % end if(next_index < 1)
    
    % Get the value of the next spectral point
    next_value = spk(next_index);
    %disp(['next_index: ',num2str(next_index),', next_value: ',num2str(next_value),', current_value: ',num2str(current_value)]);
    %pause
    
    % Is next_value below the noise threshold?
    if(isnan(next_value))
        get_next = 0;
        % set the pointer to this noise value
        % want to know the noise values
        current_index   = next_index + 1;
    
    else  % the next value is value, is it decreasing or increasing?
       
       % Is the next_value less than the current_value?
       if(next_value < current_value)
          current_value = next_value;
          current_index = next_index;
       else
          % Is the next value going up the next peak?
          %delta_value = next_value - current_value;
          delta_value = 10.*log10(next_value) - 10.*log10(current_value);
          if(delta_value >= valley_thres)
             get_next = 0;
          end % end if(delta_value >= valley_thres)
       end % end if(next_value < current_value),
    
    end % end if(next_value < nos_max_dB)
    
end % end while(get_next)

% Place the current_index into the output variable
index_limit(1) = current_index;

%% 4. Move to the left of the max value

current_index = i_max;
current_value = spk(current_index);
next_index    = current_index;
get_next      = 1;

while(get_next)
    % decrement the index
    next_index = next_index + 1;
    
    if(next_index > length(spk))
        next_index    = length(spk);
        current_index = length(spk);
        get_next      = 0;
    end % end if(next_index < 1)
    
    % Get the value of the next spectral point
    next_value = spk(next_index);
    %disp(['next_index: ',num2str(next_index),', next_value: ',num2str(10*log10(next_value)),', current_value: ',num2str(10.*log10(current_value))]);
    %pause
    
    % Is next_value below the noise threshold?
    if(isnan(next_value))
        get_next = 0;
        % set the pointer to this noise value
        % want to know the noise values
        current_index   = next_index - 1;
        
    else
       
       % Is the next_value less than the current_value?
       if(next_value < current_value)
          current_value = next_value;
          current_index = next_index;
       else
          % Is the next value going up the next peak?
          delta_value = 10.*log10(next_value) - 10.*log10(current_value);
          if(delta_value >= valley_thres)
             get_next = 0;
          end % end if(delta_value >= valley_thres)
       end % end if(next_value < current_value),
    
    end % end if(next_value < nos_max_dB)

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
else
   %end_index = 1;
end % end if(index_limit(1) > 1)

% look to the right of the peak
if(index_limit(2) < length(spk))
   start_index = index_limit(2)+1;
   
   spk_single(start_index:length(spk_single)) = ...
      spk_single(start_index:length(spk_single)) .* NaN;
else
   %start_index = length(spk);
end % end if(index_limit(1) > 1)

