function [H inliers]=H_4_each_pair(pairs,Xsift,p1,L1,domP,Hdom,polygonA,I1)

% omP : 1:if there is a dominante plane in the scene , 0: if not

alpha=2; beta=0.5;
[r c]=size(pairs);
[r2 c2]=size(Xsift);
samp=10;

H=cell(1,r);
inliers=cell(1,r);

for i=1:r
    [i r]
    % check if any part of pair i are located in the polygon of domP
    if domP
        
        IN = inpolygon([p1(1,pairs(i,:)) p1(4,pairs(i,:))],[p1(2,pairs(i,:)) p1(5,pairs(i,:))],polygonA(1,:),polygonA(2,:));
     
        if any(IN)
             H{i}=Hdom;
             continue;
        end
    end
     
     
    dp1sL=p1(1:3,pairs(i,1))'*L1(:,pairs(i,2))/norm(L1(1:2,pairs(i,2)));
    dp1eL=p1(4:6,pairs(i,1))'*L1(:,pairs(i,2))/norm(L1(1:2,pairs(i,2)));
    
    dp2sL=p1(1:3,pairs(i,2))'*L1(:,pairs(i,1))/norm(L1(1:2,pairs(i,1)));
    dp2eL=p1(4:6,pairs(i,2))'*L1(:,pairs(i,1))/norm(L1(1:2,pairs(i,1)));
    
    case1=(sign(dp1sL)==sign(dp1eL) && sign(dp2sL)==sign(dp2eL));
    case2=(sign(dp1sL)~=sign(dp1eL) && sign(dp2sL)==sign(dp2eL));
    case3=(sign(dp1sL)==sign(dp1eL) && sign(dp2sL)~=sign(dp2eL));
    case4=(sign(dp1sL)~=sign(dp1eL) && sign(dp2sL)~=sign(dp2eL));
    
    d1=sum(repmat(L1(:,pairs(i,1)),1,c2).*Xsift(1:3,:),1)./repmat(norm(L1(1:2,pairs(i,1))),1,c2);
    d2=sum(repmat(L1(:,pairs(i,2)),1,c2).*Xsift(1:3,:),1)./repmat(norm(L1(1:2,pairs(i,2))),1,c2);
       
    if case4
        
        maxX=max([p1(1,pairs(i,1:2)) p1(4,pairs(i,1:2))]);
        minX=min([p1(1,pairs(i,1:2)) p1(4,pairs(i,1:2))]);
        maxY=max([p1(2,pairs(i,1:2)) p1(5,pairs(i,1:2))]);
        minY=min([p1(2,pairs(i,1:2)) p1(5,pairs(i,1:2))]);
        
        s1=[(maxX-minX)/2+minX (maxY-minY)/2+minY 1]';
        
        L1s=[L1(1,pairs(i,1)) L1(2,pairs(i,1)) -(L1(1,pairs(i,1))*s1(1)+L1(2,pairs(i,1))*s1(2))]';
        L2s=[L1(1,pairs(i,2)) L1(2,pairs(i,2)) -(L1(1,pairs(i,2))*s1(1)+L1(2,pairs(i,2))*s1(2))]';
        
        d1=sum(repmat(L1s,1,c2).*Xsift(1:3,:),1)./repmat(norm(L1s(1:2)),1,c2);
        d2=sum(repmat(L2s,1,c2).*Xsift(1:3,:),1)./repmat(norm(L2s(1:2)),1,c2);
        
        t1=norm(L1s(1:2));
        t2=norm(L2s(1:2));
        
        j2=all([abs(d1)<(beta*t2) ; abs(d2)<(beta*t1)],1);
        
        ds1=sqrt(sum((repmat(s1,1,c2)-Xsift(1:3,:)).^2,1));
        [nd1 md1]=sort(ds1);
        
    elseif case1
                     
        j=all([sign(d1)==sign(dp2sL);sign(d2)==sign(dp1sL)],1);
        
        if sign(dp1sL)==sign(dp2sL)
            p1t=p1([4:6 1:3],pairs(i,2));
        else
            p1t=p1(:,pairs(i,2));
        end
        diag1=cross((p1(1:3,pairs(i,1))+p1(4:6,pairs(i,1)))/2,(p1t(1:3)+p1t(4:6))/2);
        diag2=cross((p1(1:3,pairs(i,1))+p1t(1:3))/2,(p1(4:6,pairs(i,1))+p1t(4:6))/2);
        
        s1=cross(diag1,diag2);
        s1=s1/s1(3);
        
        t1=norm(p1(1:3,pairs(i,1))-p1(4:6,pairs(i,1)));
        t2=norm(p1(1:3,pairs(i,2))-p1(4:6,pairs(i,2)));
        
        ds1=sqrt(sum((repmat(s1,1,c2)-Xsift(1:3,:)).^2,1));
        
        j2=all([j; abs(ds1)<(alpha*max(t1,t2))],1);
        
        ds1(~j2)=inf;
        [nd1 md1]=sort(ds1);
        
        
    elseif case2
        
        j1=(sign(d1)==sign(dp2sL));
        
        midp1=(p1(1:3,pairs(i,1))+p1(4:6,pairs(i,1)))/2;
        t1=norm(p1(1:3,pairs(i,1))-p1(4:6,pairs(i,1)));
        t2=norm(p1(1:3,pairs(i,2))-p1(4:6,pairs(i,2)));
        
        x3=midp1(1)+(p1(5,pairs(i,1))-p1(2,pairs(i,1)))./t1*20;
        y3=midp1(2)-(p1(4,pairs(i,1))-p1(1,pairs(i,1)))./t1*20;
        
        orL=cross(midp1,[x3 y3 1]');
        d2=sum(repmat(orL,1,c2).*Xsift(1:3,:),1)./repmat(norm(orL(1:2)),1,c2);
        
        j2=all([j1; abs(d1)<alpha*t2 ; abs(d2)<beta*t1],1);
        
        ds1=sqrt(sum((repmat(midp1,1,c2)-Xsift(1:3,:)).^2,1));
      
        ds1(~j2)=inf;
        [nd1 md1]=sort(ds1);
        
    elseif case3
        
        j1=(sign(d2)==sign(dp1sL));
        
        midp2=(p1(1:3,pairs(i,2))+p1(4:6,pairs(i,2)))/2;
        t1=norm(p1(1:3,pairs(i,1))-p1(4:6,pairs(i,1)));
        t2=norm(p1(1:3,pairs(i,2))-p1(4:6,pairs(i,2)));
        
        x3=midp2(1)+(p1(5,pairs(i,2))-p1(2,pairs(i,2)))./t2*20;
        y3=midp2(2)-(p1(4,pairs(i,2))-p1(1,pairs(i,2)))./t2*20;
        
        orL=cross(midp2,[x3 y3 1]');
        
        d1=sum(repmat(orL,1,c2).*Xsift(1:3,:),1)./repmat(norm(orL(1:2)),1,c2);
        
        j2=all([j1; abs(d1)<beta*t2 ; abs(d2)<alpha*t1],1);
        
        ds1=sqrt(sum((repmat(midp2,1,c2)-Xsift(1:3,:)).^2,1));
      
        ds1(~j2)=inf;
        [nd1 md1]=sort(ds1);
        
        
    end
%     [case1 case2 case3 case4]
    n=nnz(j2);
    
    if n>5
        [Hi,Inlier] = ransacfithomography(Xsift(1:3,j2),Xsift(4:6,j2), 0.01);           
    else
%         disp('not enough points to compute H ...')
        Hi=[];
        Inlier=[];
    end
    

    
    H{i}=Hi;
   
    
    j22=find(j2);
   
    inliers{i}=j22(Inlier);
   
   
    
%     figure,plot_edges_overlay_image(I1,p1(:,pairs(i,1:2)),1);hold on;hline(orL,'y');plot(Xsift(1,j2),Xsift(2,j2),'+c');
    
    
end


end
            
            
        
        
        
                    
                    
                
            
            
            
   