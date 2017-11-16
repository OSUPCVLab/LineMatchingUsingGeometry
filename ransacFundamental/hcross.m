% HCROSS - Homogeneous cross product, result normalised to s = 1.
%
% Function to form cross product between two points, or lines,
% in homogeneous coodinates.  The result is normalised to lie
% in the scale = 1 plane.
% 
% Usage: c = hcross(a,b)


function c = hcross(a,b)
c = cross(a,b);
c = c/c(3);
