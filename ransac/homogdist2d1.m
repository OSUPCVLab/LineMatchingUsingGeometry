function [H,inliers,d2,dpts] = homogdist2d1(H, x, th1, th2, T12,type)
if size(x,1)==12
    % This is based on line end-points
    T1 = T12(1:3,:);      
    T2 = T12(4:6,:);
    H  = T2\H*T1;
    
    x1 = x(1:6,:);   % Extract x1 and x2 from x
    x2 = x(7:12,:);   
    
    x1 = [ T1\x1(1:3,:); T1\x1(4:6,:)];
    x2 = [ T2\x2(1:3,:); T2\x2(4:6,:)];
    
    % Calculate, in both directions, the transfered points    
    
    Hx1    = [H*x1(1:3,:); H*x1(4:6,:)];
    invHx2 = [(H)\x2(1:3,:); (H)\x2(4:6,:)] ;   
    
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
    
    
    % difference in x,y coordinates (not orthognal distance)
    dpt1s=((x1([1 2],:)- invHx2([1 2],:)).^2);
    dpt1e=((x1([4 5],:)- invHx2([4 5],:)).^2);
    
    dpt2s=((x2([1 2],:)- Hx1([1 2],:)).^2);
    dpt2e=((x2([4 5],:)- Hx1([4 5],:)).^2);
    
    dpts=(sqrt(sum((dpt1s + dpt1e + dpt2s + dpt2e),1)));

    % Point-to-line disatnce: another form of distance error estimation
    
    d11=abs(Distance2D_pt2line(invHx2(1:3,:),invHx2(4:6,:),x1(1:3,:),3));
    d12=abs(Distance2D_pt2line(invHx2(1:3,:),invHx2(4:6,:),x1(4:6,:),3));
    
    d21=abs(Distance2D_pt2line(Hx1(1:3,:),Hx1(4:6,:),x2(1:3,:),3));
    d22=abs(Distance2D_pt2line(Hx1(1:3,:),Hx1(4:6,:),x2(4:6,:),3));
    
    d2=max([d11;d12;d21;d22]);
    inliers = find((d2<th1) & ((dpts)<th2));
end

if size(x,1)==6
    % This is based on points || line equations
    T1 = T12(1:3,:);      
    T2 = T12(4:6,:);
    H  = T2\H*T1;
    
    x1 = x(1:3,:);   % Extract x1 and x2 from x
    x2 = x(4:6,:);   
    
    x1 = T1\x1(1:3,:);
    x2 = T2\x2(1:3,:);
   
    % Calculate, in both directions, the transfered points    

    Hx1    = H*x1(1:3,:);
    invHx2 = H\x2(1:3,:);   
    
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
   
    % difference in x,y coordinates (not orthognal distance)
    
    dpt1=((x1([1 2],:)- invHx2([1 2],:)).^2);
    dpt2=((x2([1 2],:)- Hx1([1 2],:)).^2);
    dpts=(sqrt(sum((dpt1 + dpt2),1)));
    



    d2=[];
    inliers = find(dpts<th2);
end
end