% Usage: function [a, b, fit, its, in_idx] = ransacLineFitting(x, y, distThresh, agreeThresh, maxIterations)
%
% Estimate a line equation modeling x,y data using RANSAC
%
% Takes coordinates <x> and <y> corresponding to line samples
% and calculates a slope <a> and intercept <b> of a line. 
% 
% Additional Inputs:
% 'distThresh'    Defines how far points can be from the line and still be
%                 considered an inlier 
% 'agreeThresh'   The minimum percentage of inliers that must exist in
%                 the whole set for termination of RANSAC
% 'maxIterations' The maximum number of iterations allowed by RANSAC
%
% Author: Eric Psota
% Date: October 16, 2015

function [a, b, fit, its, in_idx] = GEOMransacLineFitting(x, y, distThresh, agreeThresh, maxIterations)

fit = 0;
for its_temp = 1:maxIterations
    
    rand_idx = randperm(length(x));
    rand_idx = rand_idx(1:2);
    
    a_s = (y(rand_idx(1)) - y(rand_idx(2)))/(x(rand_idx(1)) - x(rand_idx(2)));
    b_s = y(rand_idx(1)) - a_s*x(rand_idx(1));
    
    y_s = a_s*x + b_s;
    
    dist = abs(a_s*x - y + b_s)/sqrt(a_s^2 + 1);
    
    in_idx_temp = find(dist < distThresh);
    num_inliers = length(in_idx_temp);
    num_outliers = length(x) - num_inliers;
    
    fit_temp = 100*num_inliers / (num_inliers + num_outliers);
    
    if fit_temp >= agreeThresh
        fit = fit_temp;
        in_idx = in_idx_temp;
        its = its_temp;
        x_in = x(in_idx);
        x_in = [x_in(:),ones(length(x_in),1)];
        y_in = y(in_idx);
        y_in = y_in(:);
        params = (x_in'*x_in)\x_in'*y_in;
        a = params(1);
        b = params(2);
        break;
    elseif fit_temp > fit
        fit = fit_temp;
        in_idx = in_idx_temp;
        x_in = x(in_idx);
        x_in = [x_in(:),ones(length(x_in),1)];
        y_in = y(in_idx);
        y_in = y_in(:);
        params = (x_in'*x_in)\x_in'*y_in;
        a = params(1);
        b = params(2);
    else
        its = its_temp;
    end
end