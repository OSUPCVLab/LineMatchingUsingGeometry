% RANSACFITHOMOGRAPHY - fits 2D homography using RANSAC
%
% Usage:   [H, inliers] = ransacfithomography(x1, x2, t)
%
% Arguments:
%          x1  - 2xN or 3xN set of homogeneous points.  If the data is
%                2xN it is assumed the homogeneous scale factor is 1.
%          x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not. 
%                Note that point coordinates are normalised to that their
%                mean distance from the origin is sqrt(2).  The value of
%                t should be set relative to this, say in the range 
%                0.001 - 0.01  
%
% Note that it is assumed that the matching of x1 and x2 are putative and it
% is expected that a percentage of matches will be wrong.
%
% Returns:
%          H       - The 3x3 homography such that x2 = H*x1.
%          inliers - An array of indices of the elements of x1, x2 that were
%                    the inliers for the best model.
%
% See Also: ransac, homography2d, homography1d

% Copyright (c) 2004-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% February 2004 - original version
% July     2004 - error in denormalising corrected (thanks to Andrew Stein)
% August   2005 - homogdist2d modified to fit new ransac specification.

function [H, inliers] = ransacfithomography2(x1, x2, t)

    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    [rows,npts] = size(x1);
    
    % Based on line end-points
    if rows~=4 && rows~=6
        error('x1 and x2 must have 4 or 6 rows');
    end
    
    
    
    if npts < 4
        error('Must have at least 4 lines to fit homography');
    end
    
    if rows == 4    % Pad data with homogeneous scale factor of 1
        
        x1 = [x1(1:2,:); ones(1,npts);x1(3:4,:); ones(1,npts)];
        x2 = [x2(1:2,:); ones(1,npts);x2(3:4,:); ones(1,npts)];
    end
        
    % Normalise each set of points so that the origin is at centroid and
    % mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    % scale parameter is 1.  Note that 'homography2d' will also call
    % 'normalise2dpts' but the code in 'ransac' that calls the distance
    % function will not - so it is best that we normalise beforehand.
    
    % This is based on line segments end-points:
    [x1, T1,T1s] = normalise2dpts(x1);
    [x2, T2,T2s] = normalise2dpts(x2);
%     T1=eye(3);  T2=eye(3);
    T12=[T1;T2];
    s = 4;  % Minimum No of points needed to fit a homography.
    
    fittingfn = @homography2d;
    distfn    = @homogdist2d;
    degenfn   = @isdegenerate;
    % x1 and x2 are 'stacked' to create a 6xN array for ransac
    [H, inliers] = ransac([x1; x2],T12, fittingfn, distfn, degenfn, s, t);
    
    % Now do a final least squares fit on the data points considered to
    % be inliers.
    H = homography2d(x1(:,inliers), x2(:,inliers));
    
    % Denormalise
    H = T2\H*T1; 

%     H = T1'*T1s'*H*(inv(T2))'*(inv(T2s))';
%----------------------------------------------------------------------
% Function to evaluate the symmetric transfer error of a homography with
% respect to a set of matched points as needed by RANSAC.

function [inliers, H] = homogdist2d(H, x, t,T12)
    
    % This is based on line end-points
    T1 = T12(1:3,:);      
    T2 = T12(4:6,:);
    H  = T2\H*T1;
    
    x1 = x(1:6,:);   % Extract x1 and x2 from x
    x2 = x(7:12,:);   
    
    x1 = [ T1\x1(1:3,:); T1\x1(4:6,:)];
    x2 = [ T2\x2(1:3,:); T2\x2(4:6,:)];
    
%     x1_l=cross(x1(1:3,:),x1(4:6,:));
%     x2_l=cross(x2(1:3,:),x2(4:6,:));
    
    % Calculate, in both directions, the transfered points    

    Hx1    = [H*x1(1:3,:); H*x1(4:6,:)];
    invHx2 = [(H)\x2(1:3,:); (H)\x2(4:6,:)] ;   
    
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 

    projL_into1=cross(invHx2(1:3,:),invHx2(4:6,:));
    projL_into2=cross(Hx1(1:3,:),Hx1(4:6,:));
    
    % line segments length
    l1_seg=sqrt((x1(1,:)-x1(4,:)).^2+(x1(2,:)-x1(5,:)).^2);
    l2_seg=sqrt((x2(1,:)-x2(4,:)).^2+(x2(2,:)-x2(5,:)).^2);
    
    l1_seg_p=sqrt((invHx2(1,:)-invHx2(4,:)).^2+(invHx2(2,:)-invHx2(5,:)).^2);
    l2_seg_p=sqrt((Hx1(1,:)-   Hx1(4,:)).^2+(Hx1(2,:)-Hx1(5,:)).^2);
    
    maxL1=max([l1_seg;l1_seg_p]);
    maxL2=max([l2_seg;l2_seg_p]);
    
    % difference in x,y coordinates (not orthognal distance)
    dpt1=sum(sqrt((x1([1 2 4 5],:)- invHx2([1 2 4 5],:)).^2),1);
    dpt2=sum(sqrt((x2([1 2 4 5],:)-    Hx1([1 2 4 5],:)).^2),1);
    
    dpts=dpt1+dpt2;

    % length of projected lines (sqrt(mx.^2 +my.^2)
    Len_into1=sqrt(sum(projL_into1(1:2,:).^2,1));
    Len_into2=sqrt(sum(projL_into2(1:2,:).^2,1));
    
    % Another form of distance error estimation
    
    d11=abs(Distance2D_pt2line(x1(1:3,:),x1(4:6,:),invHx2(1:3,:),1));
    d12=abs(Distance2D_pt2line(x1(1:3,:),x1(4:6,:),invHx2(4:6,:),1));
    E_im1=d11+d12;
    
    d21=abs(Distance2D_pt2line(x2(1:3,:),x2(4:6,:),Hx1(1:3,:),1));
    d22=abs(Distance2D_pt2line(x2(1:3,:),x2(4:6,:),Hx1(4:6,:),1));
    E_im2=d21+d22;
    
    d2=(d11+d12) + (d21+d22);
    % Normalise so that the homogeneous scale parameter for all coordinates
    % is 1.
    
%     x1     = hnormalise(x1);
%     x2     = hnormalise(x2);    
%     x1_l   = hnormalise(x1_l);
%     x2_l   = hnormalise(x2_l); 
% 
%     Hx1    = hnormalise(Hx1);
%     invHx2 = hnormalise(invHx2); 
    
    
    
%     sqrt(sum(x1.^2,1))
%     sqrt(sum(x2.^2,1))
%     for i=1:size(x,2)
%         A1=[x1(1,i) x1(2,i) 1;x1(4,i) x1(5,i) 1];
%         B1=l1_seg(i)/(3*(Len_into1(i)))*[1 0.5;0.5 1];
%         
%         A2=[x2(1,i) x2(2,i) 1;x2(4,i) x2(5,i) 1];
%         B2=l2_seg(i)/(3*(Len_into2(i)))*[1 0.5;0.5 1];
%         
%         E1=projL_into1(:,i)'*(A1'*B1*A1)*projL_into1(:,i);
%         E2=projL_into2(:,i)'*(A2'*B2*A2)*projL_into2(:,i);
% 
%         d2(i) = abs(E1)  + abs(E2);
%         
%     end
% 
%     d2
%     inliers = find(abs(d2) < t);    
%--------------------------------------    
%         d2 = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2)
        inliers = find(abs(d2) < t & dpts<100);    
% --------------------------------
% Function to determine if a set of 4 pairs of matched  points give rise
% to a degeneracy in the calculation of a homography as needed by RANSAC.
% This involves testing whether any 3 of the 4 points in each set is
% colinear. 
     
function r = isdegenerate(x)

%     x1 = x(1:3,:);    % Extract line1 and line2 from x (based on line segment equation)
%     x2 = x(4:6,:);    

    x1 = x(1:6,:);     % Extract line1 and line2 from x (based on line end-points) 
    x2 = x(7:12,:);    

    r = ...
    isparallel(x1(:,1),x1(:,2),x1(:,3)) | ...
    isparallel(x1(:,1),x1(:,2),x1(:,4)) | ...
    isparallel(x1(:,1),x1(:,3),x1(:,4)) | ...
    isparallel(x1(:,2),x1(:,3),x1(:,4)) | ...
    isparallel(x2(:,1),x2(:,2),x2(:,3)) | ...
    isparallel(x2(:,1),x2(:,2),x2(:,4)) | ...
    isparallel(x2(:,1),x2(:,3),x2(:,4)) | ...
    isparallel(x2(:,2),x2(:,3),x2(:,4));
   