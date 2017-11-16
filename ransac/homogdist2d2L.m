function [H,inliers,d1,dpts] = homogdist2d2L(H, x, l , th1, th2, T12, type)

    % This is based on line end-points
    T1 = T12(1:3,:);      
    T2 = T12(4:6,:);
    
    Hm  = T2\H*T1;

    x1 = x(1:6,:);   % Extract x1 and x2 from x
    x2 = x(7:12,:);   
    
    l1 = l(1:3,:);   % Extract l1 and l2 from l
    l2 = l(4:6,:);   

    l1 = [ T1\l1(1:3,:)];
    l2 = [ T2\l2(1:3,:)];

    % Calculate, in both directions, the transfered points    
    % (l1, l2) lines used to compute H
    

    Hx1    = [Hm*x1(1:3,:); Hm*x1(4:6,:)];
    invHx2 = [(Hm)\x2(1:3,:); (Hm)\x2(4:6,:)] ;

    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
    
    
    % Transformed Line equations  using H: lines could be those used to
    % compute H (l1,l2), or lines computed from transformed points (L1,L2)
    

    % difference in x,y coordinates (not orthognal distance)
    dpt1=sum(sqrt((x1([1 2 4 5],:)- invHx2([1 2 4 5],:)).^2),1);
    dpt2=sum(sqrt((x2([1 2 4 5],:)-    Hx1([1 2 4 5],:)).^2),1);
%     [dpt1;dpt2]
    dpts=(dpt1 + dpt2);
%%
%%
    % Point-to-line disatnce: 1st-form of distance error estimation (using
    % lines computed from points [before they been transformed] then apply
    % tranformation to lines using computed H L1<--->L2)
    
    L1    = hnormalise(l1);
    L2    = hnormalise(l2);
    
    invHtL1    = [(Hm')\L1(1:3,:)];
    HtL2       = [(Hm')*L2(1:3,:)];
    invHtL1    = hnormalise(invHtL1);
    HtL2       = hnormalise(HtL2);

    d11=abs(Distance2D_pt2line(invHx2(1:3,:),invHx2(4:6,:),HtL2,x1(1:3,:),3));
    d12=abs(Distance2D_pt2line(invHx2(1:3,:),invHx2(4:6,:),HtL2,x1(4:6,:),3));
    d13=abs(Distance2D_pt2line(Hx1(1:3,:),Hx1(4:6,:),invHtL1,x2(1:3,:),3));
    d14=abs(Distance2D_pt2line(Hx1(1:3,:),Hx1(4:6,:),invHtL1,x2(4:6,:),3));
    [d11;d12;d13;d14];
    d1=max([d11;d12;d13;d14]);
    inliers = find((d1<th1) & ((dpts)<th2));
%%
%%
%     % Point-to-line disatnce: 2nd-form of distance error estimation (using
%     % lines computed from points [after they been transformed in both
%     % directions x1<--->x2] )
%     Lx1    = cross(invHx2(1:3,:), invHx2(4:6,:));
%     Lx2    = cross(Hx1(1:3,:), Hx1(4:6,:));
%     
%     Lx1    = hnormalise(Lx1);
%     Lx2    = hnormalise(Lx2);
%     
%     d21=abs(Distance2D_pt2line(invHx2(1:3,:),invHx2(4:6,:),Lx1,x1(1:3,:),3));
%     d22=abs(Distance2D_pt2line(invHx2(1:3,:),invHx2(4:6,:),Lx1,x1(4:6,:),3));
%     d23=abs(Distance2D_pt2line(Hx1(1:3,:),Hx1(4:6,:),Lx2,x2(1:3,:),3));
%     d24=abs(Distance2D_pt2line(Hx1(1:3,:),Hx1(4:6,:),Lx2,x2(4:6,:),3));
%     [d21;d22;d23;d24]
%     d2=max([d21;d22;d23;d24]);



end