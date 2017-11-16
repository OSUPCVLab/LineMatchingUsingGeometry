function [comb1 comb2]=FormCombOfTwoLines(p1,p2,L1,L2,angles1,angles2,w1,w2)



[r1 c1]=size(L1);
[r2 c2]=size(L2);

comb1=combnk(1:c1,2);
comb2=combnk(1:c2,2);

[r1 c1]=size(comb1);
[r2 c2]=size(comb2);

idx1=false(r1,1);
idx2=false(r2,1);

inPoint1=zeros(r1,3);
inPoint2=zeros(r2,3);

for i=1:r1
   
    ang1 = angles1(comb1(i,1),comb1(i,2));
    
    isPar1=min(ang1,abs(180-ang1));
    
    if isPar1>30
        
        inPoint1(i,:)=cross(L1(:,comb1(i,1)),L1(:,comb1(i,2)));
        inPoint1(i,:)=inPoint1(i,:)./repmat(inPoint1(i,3),1,3);
        
        cond1a=inPoint1(i,:)>0;
        cond1b=inPoint1(i,:)<=([w1([2 1]) 1]);
    
        cond1=all([cond1a cond1b],2);
        
        if cond1
            idx1(i,1)=true;
        end
        
    end   
end

for i=1:r2
   
    ang2 = angles2(comb2(i,1),comb2(i,2));
    
    isPar2=min(ang2,abs(180-ang2));
    
    if isPar2>30
        
        inPoint2(i,:)=cross(L2(:,comb2(i,1)),L2(:,comb2(i,2)));
        inPoint2(i,:)=inPoint2(i,:)./repmat(inPoint2(i,3),1,3);
        
        cond2a=inPoint2(i,:)>0;
        cond2b=inPoint2(i,:)<=([w2([2 1]) 1]);
    
        cond2=all([cond2a cond2b],2);
        
        if cond2
            idx2(i,1)=true;
        end
        
    end   
end


index1=find(idx1);
index2=find(idx2);


comb1=comb1(index1,:);
comb2=comb2(index2,:);

end