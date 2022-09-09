function pmf = stats_marg(image, range)
% Input     : image (Original Image)
%           : range (range and step value of histogram calculation, e.g. 0:255)
% Output    : pmf (Probability Mass Function)
    loc_pmf = hist(image,range);
    loc_pmf = loc_pmf/sum(loc_pmf);
    
    % update value
    pmf = loc_pmf;
end