
clc
clear all
 
fprintf('Start at ')
fprintf(datestr(now))
fprintf('.\n\n\n')
 
 
%# of sequences
numimg = 10;
str1 = 'fig7.9.';
str2 = '.gif';
 
for seq = 1:numimg
    if seq < 10
        strimg1 = strcat(str1,int2str(seq),str2);
        strimg2 = strcat(str1,int2str(seq+1),str2);
    else
        if seq == 10
            strimg1 = strcat(str1,int2str(seq),str2);
        else
            strimg1 = strcat(str1,int2str(seq+1),str2);
        end
        strimg2 = strcat(str1,int2str(seq),str2);
    end
    
    %%%% Create Gaussian Pyramids %%%%
    % # of Gaussian Pyramids levels
    level = 2;
    I1_o = 0;
    I2_o = 0;
    I_1 = 0;
    I_2 = 0;
    I1 = 0;
    I2 = 0;
    Ii = 0;
    Ip = 0;
    % read images
    I1_o = im2double(imread(strimg1));
    I2_o = im2double(imread(strimg2));
 
    [R_o,C_o] = size(I1_o);
 
    % initialize parameters
    ac = [0;0;0;0;0;0;0;0];
 
    % track computation
    rec(1,1:7) = 0; %[# of loops, ac parameters]
 
    for z = level:-1:0
        %create Gaussian Pyramid. Starts from the lowest level
        for i = 1:R_o
            for j = 1:C_o
                if (i == 1) | ((i/2^z)==floor(i/2^z))
                    if (j == 1) | ((j/2^z)==floor(j/2^z))
                        I_1(ceil(i/2^z),ceil(j/2^z)) = I1_o(i,j);
                        I_2(ceil(i/2^z),ceil(j/2^z)) = I2_o(i,j);
                    end
                end
            end
        end
 
        [R,C] = size(I_1);
 
        %%% Image Smoothing by Neighborhood Everage %%%
        nR = 1;
        for i = 2:R-1
            nC = 1;
            for j = 2:C-1
                V = 0;
                V2 = 0;
                for m = 1:3
                    for n = 1:3
                        V = V + I_1(i-2+m,j-2+n)/9;
                        V2 = V2 + I_2(i-2+m,j-2+n)/9;
                    end
                end
                I1(nR,nC) = round(V*255);
                I2(nR,nC) = round(V2*255);
                nC = nC+1;
            end
            nR = nR+1;
        end
        [r,c] = size(I1);
 
        %%%% initial loop %%%%
        Ip = I1;
        ac(3,1) = ac(3,1)*2;
        ac(6,1) = ac(6,1)*2;
 
        %%%%%%% Transform image by ac %%%%%%%
        Ii(1:r,1:c) = 0;
 
        for m = 1:r
            for n = 1:c
                %%%% compute new X %%%%
                Xp = [n;m];
                X(1,1) = ac(3)+(ac(1)+1)*Xp(1)+ac(2)*Xp(2)+ac(7)*(Xp(1)^2)+ac(8)*Xp(1)*Xp(2);
                X(2,1) = ac(6)+ac(4)*Xp(1)+(ac(5)+1)*Xp(2)+ac(7)*Xp(1)*Xp(2)+ac(8)*(Xp(2)^2);
                %%%% transform image %%%%
                if X(1,1) > 1 && X(1,1) < c && X(2,1) > 1  && X(2,1)< r
                    %%%% interpolate intensity by area based %%%%
                    dx = X(1,1) - floor(X(1,1));
                    dy = X(2,1) - floor(X(2,1));
                    Ii(m,n) = (I1(floor(X(2,1)),floor(X(1,1)))*((1-dx)*(1-dy)))+ (I1(floor(X(2,1)),ceil(X(1,1)))*((dx)*(1-dy)))+ (I1(ceil(X(2,1)),floor(X(1,1)))*((1-dx)*(dy)))+ (I1(ceil(X(2,1)),ceil(X(1,1)))*((dx)*(dy)));
                else
                    Ii(m,n) = -1;
                end
            end
        end
        Ip = Ii;
        %figure, imshow(Ip,[])
 
        %%%% start loop %%%%%%%%%%
        %tolerance
        ep = 1;
        loop = 0;
        while ep > 10^(-6) & loop < 200
            loop = loop+1;
            B(1:8,1) = 0;
            A(1:8,1:8) = 0;
            for i = 2:r
                    for j = 2:c
                        if Ip(i,j) > -1
                            It = (Ip(i,j) - I2(i,j));
                            Iy = I2(i,j) - I2(i-1,j);
                            Ix = I2(i,j) - I2(i,j-1);
                            dI(1:2,1) = [Ix;Iy];
                            X(1:2,1:8) = [j, i, 1, 0, 0, 0, (j^2), j*i; 0, 0, 0, j, i, 1, j*i, (i^2)];
                            A = A + X'* dI* dI'* X;
                            B = B + (X'* dI)* It;
                        end
                    end
            end
            a0 = -inv(A)*B;
            
            x = 1;
            y = 1;
            % update
            t(3,1) = ac(3,1)+a0(3,1)+ac(1,1)*a0(3,1)+ac(2,1)*a0(6,1)+ac(7,1)*a0(3,1)^2+ac(8,1)*a0(3,1)*a0(6,1);
            t(1,1) = ac(8,1)*x*a0(6,1)+2*ac(7,1)*x*a0(3,1)+ac(1,1)*x*a0(1,1)+ac(2,1)*a0(4,1)*x+x*a0(1,1)+ac(1,1)*x+2*ac(7,1)*x*a0(1,1)*a0(3,1)+ac(8,1)*a0(3,1)*a0(4,1)*x;
            t(2,1) = ac(2,1)*y*a0(5,1)+ac(1,1)*a0(2,1)*y+a0(2,1)*y+2*ac(7,1)*a0(2,1)*y*a0(3,1)+ac(8,1)*a0(2,1)*y*a0(6,1)+ac(8,1)*a0(3,1)*y*a0(5,1)+ac(2,1)*y+ac(8,1)*a0(3,1)*y;
            tmp1 = 2*ac(7,1)*x^2*a0(1,1)+ac(8,1)*a0(4,1)*x^2+ac(1,1)*a0(7,1)*x^2+ac(7,1)*x^2*a0(1,1)^2+a0(7,1)*x^2+2*ac(7,1)*a0(3,1)*a0(7,1)*x^2+ac(8,1)*x^2*a0(1,1)*a0(4,1)+ac(8,1)*a0(7,1)*x^2*a0(6,1)+ac(7,1)*x^2;
            tmp2 = ac(7,1)*x*y+a0(7,1)*x*y+a0(4,1)*ac(8,1)*x*y+a0(5,1)*ac(7,1)*x*y+a0(7,1)*x*y*ac(5,1)+a0(7,1)*x*ac(1,1)*y+a0(7,1)*x*ac(1,1)*y*ac(5,1)+a0(7,1)*ac(2,1)*y*ac(4,1)*x+a0(7,1)*ac(3,1)*ac(7,1)*x*y+a0(7,1)*ac(8,1)*x*y*ac(6,1)+2*a0(8,1)*ac(4,1)*x*y+2*a0(8,1)*ac(4,1)*x*y*ac(5,1)+2*a0(8,1)*ac(6,1)*ac(7,1)*x*y;
            t(7,1) = (tmp1 + tmp2)/2;
            tmp1 = a0(8,1)*x*y+ac(8,1)*x*y+ac(8,1)*a0(2,1)*y*a0(4,1)*x+ac(1,1)*a0(8,1)*x*y+ac(2,1)*a0(7,1)*x*y+2*ac(7,1)*x*a0(2,1)*y+ac(8,1)*x*y*a0(5,1)+ac(8,1)*x*a0(1,1)*y+2*ac(7,1)*x*a0(1,1)*a0(2,1)*y+2*ac(7,1)*a0(3,1)*a0(8,1)*x*y+ac(8,1)*x*a0(1,1)*y*a0(5,1)+ ac(8,1)*a0(3,1)*a0(7,1)*x*y+ac(8,1)*a0(8,1)*x*y*a0(6,1);
            tmp2 = a0(7,1)*ac(2,1)*y^2+a0(8,1)*y^2+2*a0(8,1)*y^2*ac(5,1)+a0(7,1)*ac(2,1)*y^2*ac(5,1)+ac(8,1)*y^2+a0(7,1)*ac(3,1)*ac(8,1)*y^2+2*a0(8,1)*ac(6,1)*ac(8,1)*y^2+a0(5,1)*ac(8,1)*y^2+2*a0(8,1)*ac(6,1)*ac(8,1)*y^2+a0(5,1)*ac(8,1)*y^2+a0(8,1)*y^2*ac(5,1)^2;
            t(8,1) = (tmp1 + tmp2)/2;
            t(6,1) = a0(6,1)+ac(6,1)+ac(4,1)*a0(3,1)+ac(5,1)*a0(6,1)+ac(7,1)*a0(3,1)*a0(6,1)+ac(8,1)*a0(6,1)^2;
            t(4,1) = ac(7,1)*x*a0(6,1)+a0(4,1)*x+ac(4,1)*x*a0(1,1)+ac(5,1)*a0(4,1)*x+ac(7,1)*a0(3,1)*a0(4,1)*x+ac(7,1)*x*a0(1,1)*a0(6,1)+ac(4,1)*x+2*ac(8,1)*a0(4,1)*x*a0(6,1);
            t(5,1) = ac(4,1)*a0(2,1)*y+2*ac(8,1)*y*a0(6,1)+ac(5,1)*y*a0(5,1)+y*a0(5,1)+ac(7,1)*a0(2,1)*y*a0(6,1)+ac(7,1)*a0(3,1)*y*a0(5,1)+ac(5,1)*y+ac(7,1)*a0(3,1)*y+2*ac(8,1)*y*a0(5,1)*a0(6,1);
 
            ac = t;
%             ac = ac+a0;
 
            %%%%%%% Compute new Ip (intermediate image) %%%%%%%
            Ii(1:r,1:c) = 0;
 
            for m = 1:r
                for n = 1:c
                    %%%% compute new X %%%%
                    Xp = [n;m];
                    X(1,1) = ac(3)+(ac(1)+1)*Xp(1)+ac(2)*Xp(2)+ac(7)*(Xp(1)^2)+ac(8)*Xp(1)*Xp(2);
                    X(2,1) = ac(6)+ac(4)*Xp(1)+(ac(5)+1)*Xp(2)+ac(7)*Xp(1)*Xp(2)+ac(8)*(Xp(2)^2);
                    %%%% transform image %%%%
                    if X(1,1) > 1 && X(1,1) < c && X(2,1) > 1  && X(2,1)< r
                        %%%% interpolate intensity by area based %%%%
                        dx = X(1,1) - floor(X(1,1));
                        dy = X(2,1) - floor(X(2,1));
                        Ii(m,n) = (I1(floor(X(2,1)),floor(X(1,1)))*((1-dx)*(1-dy)))+ (I1(floor(X(2,1)),ceil(X(1,1)))*((dx)*(1-dy)))+ (I1(ceil(X(2,1)),floor(X(1,1)))*((1-dx)*(dy)))+ (I1(ceil(X(2,1)),ceil(X(1,1)))*((dx)*(dy)));
                    else
                        Ii(m,n) = -1;
                    end
                end
            end
            Ip = Ii;
            % compute tolerance        
            ep = a0'*a0;
    %         figure, imshow(Ip,[])
    %         figure, imshow(Ip2,[])
        end
        rec(level-z+1,1:9) = [loop,ac'];
    end
    %figure,imshow(Ip,[])
    transrec(seq,1:8) = ac';
    seq
    fprintf(datestr(now))
end
 
mosaic
 
fprintf('\n')
fprintf('Finish at ')
fprintf(datestr(now))
fprintf('\n')

%%



Im = 0;
 
maxR = r;
maxC = c;
shiftR = 1;
shiftC = 1;
 
for seq = 1:numimg
    if seq < 10
        strimg1 = strcat(str1,'0',int2str(seq-1),str2);
        strimg2 = strcat(str1,'0',int2str(seq),str2);
    else
        if seq == 10
            strimg1 = strcat(str1,'0',int2str(seq-1),str2);
        else
            strimg1 = strcat(str1,int2str(seq-1),str2);
        end
        strimg2 = strcat(str1,int2str(seq),str2);
    end
    
    I1 = round(im2double(imread(strimg1))*255);
    I2 = round(im2double(imread(strimg2))*255);
    
    if seq == 1
        Im = I1;
    end
    
    for m = 1:r
        for n = 1:c
            
            I = I2(m,n);
            
            %%%% compute new X %%%%
            Xp = [n;m];
            for loop = 1:seq
                ac = transrec(loop,1:8);
                Xp(1,1) = ac(3)+(ac(1)+1)*Xp(1)+ac(2)*Xp(2)+ac(7)*(Xp(1)^2)+ac(8)*Xp(1)*Xp(2);
                Xp(2,1) = ac(6)+ac(4)*Xp(1)+(ac(5)+1)*Xp(2)+ac(7)*Xp(1)*Xp(2)+ac(8)*(Xp(2)^2);
                
                [maxR,maxC] = size(Im);
                if floor(Xp(1))+shiftC < 1
                    shiftC = shiftC - floor(Xp(1));
                    maxC = maxC + floor(Xp(1))+shiftC;
                    Imp(1:maxR,1:shiftC) = 0;
                    Imp(1:maxR,shiftC+1:maxC) = Im;
                    Im = Imp;
                end
                if floor(Xp(2))+shiftR < 1
                    shiftR = shiftR - floor(Xp(2));
                    maxR = maxR + floor(Xp(2))+shiftR;
                    Imp(1:shiftR,1:maxC) = 0;
                    Imp(shiftR+1:maxR,1:maxC) = Im;
                    Im = Imp;
                end                
                if floor(Xp(1))+shiftC > maxC
                    maxC = floor(Xp(1))+shiftC;
                end 
                if floor(Xp(2))+shiftR > maxR
                    maxR = floor(Xp(2))+shiftR;
                end  
                
            end
            Im(floor(Xp(2))+shiftR,floor(Xp(1))+shiftC) = I;
 
        end
        
    end
    
end
figure, imshow(Im,[])

