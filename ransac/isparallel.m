% ISPARALLEL - are 3 lines parellel



function r = isparallel(l1, l2, l3)


if size(l1,1)==4
    l1=[l1(1:2,:);1;l1(3:4,:);1];
    l2=[l2(1:2,:);1;l2(3:4,:);1];
    l3=[l3(1:2,:);1;l3(3:4,:);1];
end

if size(l1,1)==6
    l1=cross(l1(1:3,:),l1(4:6,:));
    l2=cross(l2(1:3,:),l2(4:6,:));
    l3=cross(l3(1:3,:),l3(4:6,:));
end


r=abs(det([l1 l2 l3]))<eps;

end