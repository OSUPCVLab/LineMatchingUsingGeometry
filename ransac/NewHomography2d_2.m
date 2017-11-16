function [H] = NewHomography2d_2(x1,x2,T12,type)
    

    
% %     % Attempt to normalise each set of points so that the origin 
% %     % is at centroid and mean distance from origin is sqrt(2).

if (strcmp('LP', type) || strcmp('P', type))
    if size(x1,1)==6
        xp1=reshape(x1,3,2*size(x1,2));
        xp2=reshape(x2,3,2*size(x2,2));
        
%         % This if you want to use mid-point of line segments instead of
%         % both end-points (in case of the orientation is unknown)
%         xp1=[mean(x1([1 4],:));mean(x1([2 5],:));mean(x1([3 6],:))];
%         xp2=[mean(x2([1 4],:));mean(x2([2 5],:));mean(x2([3 6],:))];
    else
        xp1=x1;
        xp2=x2;
        

    end
    Npts = size(xp1,2);
    A = zeros(3*Npts,9);
    
    O = [0 0 0];
    for n = 1:Npts
        X = xp1(:,n)';
        x = xp2(1,n); y = xp2(2,n); w = xp2(3,n);
        A(3*n-2,:) = [  O  -w*X  y*X];
        A(3*n-1,:) = [ w*X   O  -x*X];
    	A(3*n  ,:) = [-y*X  x*X   O ];
    end
   
    [U,D,V] = svd(A,0); % 'Economy' decomposition for speed
    H1 = reshape(V(:,9),3,3)';
end

 %%  This section of building H is based on end-points of line segemnts correspondences
% if (strcmp('LP', type) || strcmp('L', type))
%     Npts = size(x1,2);
%     A = zeros(2*Npts,9);
%     for n = 1:Npts
%         L1=( x1(2,n)- x1(5,n));
%         L2=( x1(4,n)- x1(1,n));
%         L3=( x1(1,n)*x1(5,n)-x1(4,n)*x1(2,n));
%         
%         X1s = x2(1:3,n)'; X1e = x2(4:6,n)'; 
%         
%         A(2*n-1,:) = [ L1*X1s(1)  L2*X1s(1)  L3*X1s(1) L1*X1s(2)  L2*X1s(2)  L3*X1s(2) L1  L2  L3];
%         A(2*n  ,:) = [ L1*X1e(1)  L2*X1e(1)  L3*X1e(1) L1*X1e(2)  L2*X1e(2)  L3*X1e(2) L1  L2  L3];
% 
%     end
%     [U,D,V] = svd(A,0); % 'Economy' decomposition for speed
%     H2 = reshape(V(:,9),3,3)';
%     H2 = inv(H2');
% end
%%
if (strcmp('LP', type) || strcmp('L', type))
    
    % DO I need to de-normalize the points ???????? (NO)
%     T1 = T12(1:3,:);      
%     T2 = T12(4:6,:);
%     x1 = [ T1\x1(1:3,:); T1\x1(4:6,:)];
%     x2 = [ T2\x2(1:3,:); T2\x2(4:6,:)];
    
    Npts = size(x1,2);
    
    l1=cross(x1(4:6,:),x1(1:3,:));
    l2=cross(x2(4:6,:),x2(1:3,:));
    
   
    A = zeros(3*Npts,9);
    
    O = [0 0 0];
    for n = 1:Npts
        X = l1(:,n)';
        x = l2(1,n); y = l2(2,n); w = l2(3,n);
        A(3*n-2,:) = [  O  -w*X  y*X];
        A(3*n-1,:) = [ w*X   O  -x*X];
    	A(3*n  ,:) = [-y*X  x*X   O ];

    end
    [U,D,V] = svd(A,0); % 'Economy' decomposition for speed
    H2 = reshape(V(:,9),3,3)';
    H2 = inv(H2');
end
%%

if (strcmp('P', type)); H=H1;end
if (strcmp('L', type)); H=H2;end
if (strcmp('LP', type)); H=[H1;H2];end
end


   