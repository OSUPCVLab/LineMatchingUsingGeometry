function [ang int1 int2 x1n x2n]=CheckAngleAndIntersection2SegLines(x1,x2,I1,SegEnds)

% this function computes the clockwise angle between 2 lines (or line
% segment and line vecor) & check if thier intersection point is located
% within any segment or outside

% SegEnds:  0: if segments ends direction are unknown;
%          ~0: if direction is known

[r1 c1]=size(x1);
[r2 c2]=size(x2);

if r1==r2
    
    if r1==4
        x1=[x1(1:2);1;x1(3:4);1];
        x2=[x2(1:2);1;x2(3:4);1];
    end
    
    L1=cross(x1(1:3),x1(4:6));
    L2=cross(x2(1:3),x2(4:6));
    
    if SegEnds
        intP=cross(L1,L2);  intP=intP./intP(3);
    
        d1= sqrt(sum((x1(1:3)-x1(4:6)).^2));
        d1a=sqrt(sum((intP-x1(1:3)).^2));
        d1b=sqrt(sum((intP-x1(4:6)).^2));
        [m1 m2]=sort([d1a d1b],'descend');
%         [m1 m2]=max([d1a d1b]);

        x1n=[intP;x1(3*m2(1)-2:3*m2(1));x1(3*m2(2)-2:3*m2(2))];

        d2= sqrt(sum((x2(1:3)-x2(4:6)).^2));
        d2a=sqrt(sum((intP-x2(1:3)).^2));
        d2b=sqrt(sum((intP-x2(4:6)).^2));
        
        [m1 m2]=sort([d2a d2b],'descend');
%         [m1 m2]=max([d2a d2b]);

        x2n=[intP;x2(3*m2(1)-2:3*m2(1));x2(3*m2(2)-2:3*m2(2))];
    
        L1=cross(x1n(1:3),x1n(4:6));
        L2=cross(x2n(1:3),x2n(4:6));
        
        int1 = all(d1>[d1a d1b]);
        int2 = all(d2>[d2a d2b]);
        
    else
        int1=[];
        int2=[];
        x1n=[];
        x2n=[];
    end
        
    ang = mod( atan2( det([L1(1:2),L2(1:2)]) , dot(L1(1:2),L2(1:2)) ) , 2*pi )*180/pi;
    
    
    
%     figure(2),plot_edges_overlay_image(I1,[x1n x2n],1);hold on;plot(intP(1),intP(2),'Or');hold off
%     pause
else
    
    if r1==4
        x1=[x1(1:2);1;x1(3:4);1];

    end
    
    L1=cross(x1(1:3),x1(4:6));
    L2=x2;
    
    if ~SegEnds
        
        intP=cross(L1,L2);  intP=intP./intP(3);

        d1= norm(x1(1:3)-x1(4:6));
        d1a=norm(intP-x1(1:3));
        d1b=norm(intP-x1(4:6));

        [m1 m2]=sort([d1a d1b],'descend');
        x1n=[intP;x1(3*m2(1)-2:3*m2(1));x1(3*m2(2)-2:3*m2(2))];
   
        L1=cross(x1n(1:3),x1n(4:6));
        int1 = all(d1>[d1a d1b]);
    else
        int1=[];
    end
    
    ang = mod( atan2( det([L1(1:2),L2(1:2)]) , dot(L1(1:2),L2(1:2)) ) , 2*pi )*180/pi;
    int2=[];
    x2n=[];
    
%     figure(2),plot_edges_overlay_image(I1,x1([1 2 4 5]),1);hold on;
%             hline(L1,'r');
%             hline(x2,'g');
%             pause
end
end

    
    

