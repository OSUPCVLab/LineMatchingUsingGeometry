% Demonstration of feature matching via simple correlation, and then using
% RANSAC to estimate the fundamental matrix and at the same time identify
% (mostly) inlying matches
%
% Usage:  testfund              - Demonstrates fundamental matrix calculation
%                                 on two default images
%         testfund(im1,im2)     - Computes fundamental matrix on two supplied images
%
% Edit code as necessary to tweak parameters

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

% February 2004
% August   2005 Octave compatibility

function [F X Points Fn] = testfund(im1,im2,points)
    
    if nargin == 0
	im1 = imread('im02.jpg');
	im2 = imread('im03.jpg');
    end

%     v = version; Octave=v(1)<'5';  % Crude Octave test        
%     thresh = 500;   % Harris corner threshold
%     nonmaxrad = 3;  % Non-maximal suppression radius
%     dmax = 50;      % Maximum search distance for matching
%     w = 11;         % Window size for correlation matching
%     
%     % Find Harris corners in image1 and image2
%     [cim1, r1, c1] = harris(im1, 1, thresh, 3);
%     show(im1,1), hold on, plot(c1,r1,'r+');
% 
%     [cim2, r2, c2] = harris(im2, 1, thresh, 3);
%     show(im2,2), hold on, plot(c2,r2,'r+');
%     drawnow
% 
%     correlation = 1;  % Change this between 1 or 0 to switch between the two
%                       % matching functions below
%     
%     if correlation  % Use normalised correlation matching
% 	[m1,m2] = matchbycorrelation(im1, [r1';c1'], im2, [r2';c2'], w, dmax);
% 	
%     else            % Use monogenic phase matching
% 	nscale = 1;
% 	minWaveLength = 10;
% 	mult = 4;
% 	sigmaOnf = .2;
% 	[m1,m2] = matchbymonogenicphase(im1, [r1';c1'], im2, [r2';c2'], w, dmax,...
% 					nscale, minWaveLength, mult, sigmaOnf);
%     end   
    
%     points=dlmread('demo_ASIFT_Win\matchings.txt');
%     points=points(2:end,:);
    
    m1=points(1:2,:);
    m2=points(3:4,:);
    
    % Display putative matches
    figure, subplot(2,1,1);imshow(im1), title('Putative matches in image 1');hold on;
                plot(m1(1,:),m1(2,:),'r+');hold off;
            subplot(2,1,2);imshow(im2), title('Putative matches in image 2');hold on;
                plot(m2(1,:),m2(2,:),'r+');hold off;
           


    % Assemble homogeneous feature coordinates for fitting of the
    % fundamental matrix, note that [x,y] corresponds to [col, row]
    x1 = [m1(1,:); m1(2,:); ones(1,length(m1))];
    x2 = [m2(1,:); m2(2,:); ones(1,length(m1))];    
    
    t = .001;  % Distance threshold for deciding outliers
    
    % Change the commenting on the lines below to switch between the use
    % of 7 or 8 point fundamental matrix solutions, or affine fundamental
    % matrix solution.
%   [F, inliers] = ransacfitfundmatrix7(x1, x2, t, 1);    
    [F, inliers Fn] = ransacfitfundmatrix(x1, x2, t, 1);
%   [F, inliers] = ransacfitaffinefund(x1, x2, t, 1);    
    
    fprintf('Number of inliers was %d (%d%%) \n', ...
	    length(inliers),round(100*length(inliers)/length(m1)))
    fprintf('Number of putative matches was %d \n', length(m1))        
    
    % Display both images overlayed with inlying matched feature points
    

   
    im3 = appendimages(im1,im2);
    figure,imshow(im3);hold on;
    cols1 = size(im1,2);
    for i = inliers
        line([m1(1,i) m2(1,i)+cols1],[m1(2,i) m2(2,i)], 'Color', 'w');
    end
    hold off;
    
    
%     response = input('Step through each epipolar line [y/n]?\n','s');
%     if response == 'n'
% 	return
%     end    
    
    % Step through each matched pair of points and display the
    % corresponding epipolar lines on the two images.
    
    l2 = F*x1;    % Epipolar lines in image2
    l1 = F'*x2;   % Epipolar lines in image1
    
    % Solve for epipoles
    [U,D,V] = svd(F);
    e1 = hnormalise(V(:,3));
    e2 = hnormalise(U(:,3));
    
%     figure,
%     id=1;
%     for n = inliers
%         subplot(2,1,1), imshow(im1), hold on,
%         plot(x1(1,n),x1(2,n),'r+');
%         hline(l1(:,n)); plot(e1(1), e1(2), 'g*');hold off;
%         title(['inlier #  ',num2str(id),'  ',num2str(length(inliers))])
%         subplot(2,1,2), imshow(im2), hold on;
%         plot(x2(1,n),x2(2,n),'r+');
%         hline(l2(:,n)); plot(e2(1), e2(2), 'g*');hold off;
% 
%         pause;
%         id=id+1;
%     end
    
    X=[x1(:,inliers);x2(:,inliers)];
    
    Points=[x1;x2];
    

    fprintf('                                         \n');
    
    
    