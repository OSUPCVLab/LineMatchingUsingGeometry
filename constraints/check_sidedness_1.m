function [over]=check_sidedness_1(x1,x2,type1)

% This function checks the overlap constraint

if strcmp('P2P', type1)
    [r1 c1]=size(x1);
    [r2 c2]=size(x2);
    
    

    if r1==4
        x1=[x1(1:2,:);ones(1,c1);x1(3:4,:);ones(1,c1)];
    elseif r1==2
        x1=[x1(1:2,:);ones(1,c1)];        
    end
    
    if r2==4
        x2=[x2(1:2,:);ones(1,c2);x2(3:4,:);ones(1,c2)];
    end

    L2=cross(x2(1:3,:),x2(4:6,:));
    
    cond1=sign(x1(1:3,:)'*L2);
    cond2=sign(x1(4:6,:)'*L2);
    
    cond12=[cond1;cond2];
    
    i1=find(cond12==0);
    
    if ~isempty(i1)
        i2=repmat(sum(cond12),2,1);
        cond12(i1)=i2(i1);
    end
    
    over=~(all([cond12(1,:)==cond12(2,:) cond12(:,1)'==cond12(:,2)']));

    
    
elseif strcmp('P2L', type1)
    
    [r1 c1]=size(x1);

    
    

    if r1==4
        x1=[x1(1:2,:);ones(1,c1);x1(3:4,:);ones(1,c1)];
    elseif r1==2
        x1=[x1(1:2,:);ones(1,c1)];        
    end
    
    L2=x2;
    
    cond1=sign(x1(1:3,:)'*L2);
    cond2=sign(x1(4:6,:)'*L2);
    
    cond12=[cond1;cond2];
    
    i1=find(cond12==0);
    
    if ~isempty(i1)
        i2=repmat(sum(cond12),2,1);
        cond12(i1)=i2(i1);
    end
    
    over=~(all([cond12(1,:)==cond12(2,:) cond12(:,1)'==cond12(:,2)']));

end

end