% Copyright (c) Mohammed Al-Shahri
% Photogrammetric Computer Vision Lab.
% The Ohio State University

function [p1 p2 L1 L2 ang1 ang2]=PreProcess_LinesDetectionAndMerging1(inname,outname,w1,w2,LE)

% Detection lines : use 20 as min length of segment to be detected
% Delete or Merging : use 30 to combine parallel line within 30 pixels and
% keep the longer, which is more stronger

% [Pts] = DetectLines(Imgs,1,30);
% p1=Pts{1};
% p2=Pts{2};

if strcmp(LE, 'LP')
    % line extraction using "LSD method , used for LP method in my paper"
    [p1] = DetectLines_LSD(inname{1},outname{1},20);
    [p2] = DetectLines_LSD(inname{2},outname{2},20);
    
elseif strcmp(LE, 'LS')
    % line extraction using "Line signature paper method, used for LS method in my paper"
    [p1] = DetectLines_LSig(inname{1},outname{1},20);
    [p2] = DetectLines_LSig(inname{2},outname{2},20);
    
else
    disp('wrong line extraction method')
end


% removing lines extracted from the boundary of the image I1
outXp1s=all([p1(1,:)<10;p1(3,:)<10],1);
outXp1e=all([p1(1,:)>w1(2)-10;p1(3,:)>w1(2)-10],1);
outYp1s=all([p1(2,:)<10;p1(4,:)<10],1);
outYp1e=all([p1(2,:)>w1(1)-10;p1(4,:)>w1(1)-10],1);

outp1=any([outXp1s;outXp1e;outYp1s;outYp1e],1);
p1(:,outp1)=[];

% removing lines extracted from the boundary of the image I2
outXp2s=all([p2(1,:)<10;p2(3,:)<10],1);
outXp2e=all([p2(1,:)>w2(2)-10;p2(3,:)>w2(2)-10],1);
outYp2s=all([p2(2,:)<10;p2(4,:)<10],1);
outYp2e=all([p2(2,:)>w2(1)-10;p2(4,:)>w2(1)-10],1);

outp2=any([outXp2s;outXp2e;outYp2s;outYp2e],1);
p2(:,outp2)=[];


% angtol=0.75*pi/180;
% linkrad=12;
% linedevtol=2;
% sizep1=size(p1,2)+1;
% sizep2=size(p2,2)+1;
% while (sizep1~=size(p1,2) || sizep2~=size(p2,2))
%     
%     sizep1=size(p1,2);
%     sizep2=size(p2,2);
%     
%     p1 = (mergeseg(p1', angtol, linkrad, linedevtol))';
%     p2 = (mergeseg(p2', angtol, linkrad, linedevtol))';
%            
% end


if size(p1,1)==4
    p1=[p1(1:2,:);ones(1,size(p1,2));p1(3:4,:);ones(1,size(p1,2))];
    p2=[p2(1:2,:);ones(1,size(p2,2));p2(3:4,:);ones(1,size(p2,2))];
end

[L1,L2] = IntersectionPoint4lines_2(p1,p2);


% step2: angles between lines
[ang1]=ComputeAnglesBetweenLines(p1);
[ang2]=ComputeAnglesBetweenLines(p2);





       
       
end