% HOMOGRAPHY2D - computes 2D homography
%
% Usage:   H = homography2d(x1, x2)
%          H = homography2d(x)
%
% Arguments:
%          x1  - 3xN set of homogeneous points
%          x2  - 3xN set of homogeneous points such that x1<->x2
%         
%           x  - If a single argument is supplied it is assumed that it
%                is in the form x = [x1; x2]
% Returns:
%          H - the 3x3 homography such that x2 = H*x1
%
% This code follows the normalised direct linear transformation 
% algorithm given by Hartley and Zisserman "Multiple View Geometry in
% Computer Vision" p92.
%

% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% May 2003  - Original version.
% Feb 2004  - Single argument allowed for to enable use with RANSAC.
% Feb 2005  - SVD changed to 'Economy' decomposition (thanks to Paul O'Leary)

function [H,l1,l2] = NewHomography2d_1(x1,x2,T12,type)
    

    
% %     % Attempt to normalise each set of points so that the origin 
% %     % is at centroid and mean distance from origin is sqrt(2).
%     [x1, T1,T1s] = normalise2dpts(x1);
%     [x2, T2,T2s] = normalise2dpts(x2);
if type =='P'
    if size(x1,1)==6
        x1=reshape(x1,3,2*size(x1,2));
        x2=reshape(x2,3,2*size(x2,2));
%         x1=[mean(x1([1 4],:)); mean(x1([2 5],:));x1(3,:)]
%         x2=[mean(x2([1 4],:)); mean(x2([2 5],:));x2(3,:)]
    end
    Npts = size(x1,2);
    A = zeros(3*Npts,9);
    
    O = [0 0 0];
    for n = 1:Npts
        X = x1(:,n)';
        x = x2(1,n); y = x2(2,n); w = x2(3,n);
        A(3*n-2,:) = [  O  -w*X  y*X];
        A(3*n-1,:) = [ w*X   O  -x*X];
    	A(3*n  ,:) = [-y*X  x*X   O ];
    end
    l1=[]; l2=[];
end

    
%     % 2. This section of building H is based on end-points of line segemnts correspondences
% if type ==2
%     Npts = size(x1,2);
%     A = zeros(2*Npts,9);
%     for n = 1:Npts
%         L1=( x2(2,n)- x2(5,n));
%         L2=( x2(4,n)- x2(1,n));
%         L3=( x2(1,n)*x2(5,n)-x2(4,n)*x2(2,n));
%         
%         X1s = x1(1:3,n)'; X1e = x1(4:6,n)'; 
%         
% %         A(2*n-1,:) = [ L1*X1s(1)  L2*X1s(1)  L3*X1s(1) L1*X1s(2)  L2*X1s(2)  L3*X1s(2) L1  L2  L3];
% %         A(2*n  ,:) = [ L1*X1e(1)  L2*X1e(1)  L3*X1e(1) L1*X1e(2)  L2*X1e(2)  L3*X1e(2) L1  L2  L3];
%         A(2*n-1,:) = [ L1*X1s  L2*X1s  L3*X1s];
%         A(2*n  ,:) = [ L1*X1e  L2*X1e  L3*X1e];
%     end
% end

if type =='L'
    
    Npts = size(x1,2);
    
    l1=cross(x1(1:3,:),x1(4:6,:));
    l2=cross(x2(1:3,:),x2(4:6,:));
    
    l11=l1./repmat(l1(end,:),3,1);
    l22=l2./repmat(l2(end,:),3,1);
    
    A = zeros(3*Npts,9);
    
    O = [0 0 0];
    for n = 1:Npts
        X = l11(:,n)';
        x = l22(1,n); y = l22(2,n); w = l22(3,n);
        A(3*n-2,:) = [  O  -w*X  y*X];
        A(3*n-1,:) = [ w*X   O  -x*X];
    	A(3*n  ,:) = [-y*X  x*X   O ];

    end

end


[U,D,V] = svd(A,0); % 'Economy' decomposition for speed
H = reshape(V(:,9),3,3)';

if type =='L'; H=inv(H');end


end


   