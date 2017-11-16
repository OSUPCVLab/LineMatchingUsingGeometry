function [H] = NewHomography2d_2L(l1,l2,type)
    
if (strcmp('LP', type) || strcmp('L', type))
    
    Npts = size(l1,2);
    
  
%     l1=l1./repmat(l1(end,:),3,1);
%     l2=l2./repmat(l2(end,:),3,1);
    
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
    H = reshape(V(:,9),3,3)';
    H = inv(H');
end

end


   