function [L1,L2] = IntersectionPoint4lines_2(p1,p2)

[r1 c1]=size(p1);
[r2 c2]=size(p2);

if r1==4
    p1=[p1(1:2,:);ones(1,c1);p1(3:4,:);ones(1,c1)];
    p2=[p2(1:2,:);ones(1,c2);p2(3:4,:);ones(1,c2)];
end

L1=cross(p1(1:3,:),p1(4:6,:));
L2=cross(p2(1:3,:),p2(4:6,:));


end
    