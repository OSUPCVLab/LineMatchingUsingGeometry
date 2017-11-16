function [Angles]=ComputeAnglesBetweenLines(P)



% this function computes the smallest angle between every two lines

if size(P,1)==6
    L=cross(P(1:3,:),P(4:6,:));
    P=P([1 2 4 5],:);
elseif size(P,1)==4
    L=cross([P(1:2,:);ones(1,size(P,2))],[P(3:4,:);ones(1,size(P,2))]);
end

[r c]=size(L);

Angles=zeros(c,c);

for i=1:c
    
    for j=1:c
        
                
        intp1=cross(L(:,i),L(:,j)); intp1=intp1./intp1(3);
        
%         intp2=cross(L(:,Ksub(i,2)),L(:,Ksub(j,2))); intp2=intp2./intp2(3);
        
        da1=sqrt(sum(([P(1:2,i) P(3:4,i)]-[repmat(intp1(1:2),1,2)]).^2,1));
        db1=sqrt(sum(([P(1:2,j) P(3:4,j)]-[repmat(intp1(1:2),1,2)]).^2,1));
        
%         da2=sqrt(sum(([P(1:2,Ksub(i,2)) P(3:4,Ksub(i,2))]-[repmat(intp2(1:2),1,2)]).^2,1));
%         db2=sqrt(sum(([P(1:2,Ksub(j,2)) P(3:4,Ksub(j,2))]-[repmat(intp2(1:2),1,2)]).^2,1));
        
        if (da1(1)>da1(2))
            va1=(intp1(1:2)-P(1:2,i));
        else
            va1=(intp1(1:2)-P(3:4,i));
        end
        
        if (db1(1)>db1(2))
            vb1=(intp1(1:2)-P(1:2,j));
        else
            vb1=(intp1(1:2)-P(3:4,j));
        end
        
        
        
        theta=acosd((dot(va1,vb1))./(norm(va1)*norm(vb1)));
        
        Angles(i,j)=theta;
    end
end

end