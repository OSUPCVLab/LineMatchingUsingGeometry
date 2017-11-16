
Ip = 0
I1 = double(imread(imstr(1,:)));                         	% loads image 1
[h1 w1 d1] = size(I1); 

for seq = 1:numimg-1
    
    I2 = double(imread(imstr(seq+1,:)));                       	% load images 2
    [h2 w2 d2] = size(I2);  
    
    for loop = seq:-1:1
        if seq == loop
            H = transrec(:,:,loop);
        else
            H = transrec(:,:,loop)*H;
        end
    end   
    H = H/H(3,3);
    
    % construct transformation matrix (T) 
    T = H; 

    % warps incoming corners to determine the size of the output image (in to out) 
    cp = T*[ 1 1 w2 w2 ; 1 h2 1 h2 ; 1 1 1 1 ];
    Xpr = min( [ cp(1,:) 0 ] ) : max( [cp(1,:) w1] ); % min x : max x 
    Ypr = min( [ cp(2,:) 0 ] ) : max( [cp(2,:) h1] ); % min y : max y 
    [Xp,Yp] = ndgrid(Xpr,Ypr); 
    [wp hp] = size(Xp); % = size(Yp)  

    % do backwards transform (from out to in) 
    X = T \ [ Xp(:) Yp(:) ones(wp*hp,1) ]'; % warp  
    for i = 1:size(X,2)
        X(:,i) = X(:,i)/X(3,i);
    end

    % re-sample pixel values with bilinear interpolation 
    clear Ip; 
    xI = reshape( X(1,:),wp,hp)'; 
    yI = reshape( X(2,:),wp,hp)'; 
    Ip(:,:,1) = interp2(I2(:,:,1), xI, yI, '*bilinear'); % red 
%     figure; imshow(Ip,[])

    % offset and copy original image into the warped image 
    offset =  -round( [ min( [ cp(1,:) 0 ] ) min( [ cp(2,:) 0 ] ) ] ); 
    for i = 1:h1
        for j = 1:w1
            if double(I1(i,j,:)) > 10
                Ip(i+offset(2),j+offset(1),:) = double(I1(i,j,:));
            end
        end
    end
%     Ip(1+offset(2):h1+offset(2),1+offset(1):w1+offset(1),:) = double(I1(1:h1,1:w1,:));  
    I1 = Ip;
    [h1 w1 d1] = size(I1); 
    figure; imshow(Ip,[])
    % 
end
