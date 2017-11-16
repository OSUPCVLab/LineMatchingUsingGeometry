% Copyright (c) Mohammed Al-Shahri
% Photogrammetric Computer Vision Lab.
% The Ohio State University


% Main code for line matching algorithm

%%
addpath(genpath(pwd));

%% step 1: pre-processing

% (a) read images

image1=[pwd '\images\dataset_1\house000.png'];      % view 1
image2=[pwd '\images\dataset_1\house004.png'];      % view 2
line_img1=[pwd '\images\dataset_1\img1.txt'];       % lines in view 1
line_img2=[pwd '\images\dataset_1\img2.txt'];       % lines in view 2

images={image1;image2};
line_imgs={line_img1;line_img2};

I1=imread(image1);
I2=imread(image2);

w1=size(I1);
w2=size(I2);

Imgs={I1;I2};

%%  Extarct lines & compute angles between lines 

% Line extraction method: choose one method
% (LS:Line signature), or  (LP: Leverage Point)
LE='LP';        
[p1 p2 L1 L2 ang1 ang2]=PreProcess_LinesDetectionAndMerging1(images,line_imgs,w1,w2,LE);

[matches]=test_demo_ASIFT(images); % extract ASIFT matches

[F X ] = testfund(I1,I2,matches);  % estiamte F matrix

% check if dominante plane exist in the scene
[Hdom, inlier] = ransacfithomography(X(1:3,:),X(4:6,:), 0.001); 
DomP=length(inlier)>0.2*size(X,2);   
if DomP
    k = convhull(X(1,inlier),X(2,inlier));
    polygonA=[X(1,inlier(k));X(2,inlier(k))];
else
    polygonA=[];
end


%% step 2: form pair of lines
[comb1 comb2]=FormCombOfTwoLines(p1,p2,L1,L2,ang1,ang2,w1,w2);

%% step 3: compute local H /or coplanar points for each pair

[H inliers]=H_4_each_pair(comb1,X,p1,L1,DomP,Hdom,polygonA,I1);  


%% step 4: matching process
showFig=0;
[matchVa BestMatches GraphPairs]=TopDownApproach1(F,H,comb1,comb2,L1,L2,p1,p2,ang1,ang2,I1,I2,showFig);  


%% final results

% Affinity matrix: matchVa

gr=2;  % 1 or 2  (based on cost affinity marix)
per=0.10;  % greedy method, all weights >10% of affinity marix max weight

s=max(matchVa{gr}(:,1));
j=matchVa{gr}(:,1)>per*s;

% display correspondences
h=3;        % show best 1, or 2, or 3, .... matches. 
PlotCorrespLines2(I1,I2,p1,p2,L1,L2,BestMatches{gr}(j,1:h+1))




