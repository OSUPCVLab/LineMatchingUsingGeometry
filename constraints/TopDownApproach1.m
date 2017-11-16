% Copyright (c) Mohammed Al-Shahri
% Photogrammetric Computer Vision Lab.
% The Ohio State University

function [matchVa BestMatches MatchError]=TopDownApproach1(F,H,comb1,comb2,L1,L2,p1,p2,angles1,angles2,I1,I2,showFig)

% This function measure the similarity of two corresponding pairs in two
% images based on:
% 1: epipolar geometry
% 2: overlap
% 3: the local H gatherd from set of corresponding points 

% just for the figure
curFig=gcf+1;
%------------

[S1 S2 S3]=svd(F);
e1=S3(:,end); %e1=e1./e1(end);
e2=S1(:,end); %e2=e2./e2(end);
e1e2=cross(e1,e2);


[r1 c1]=size(L1);
[r2 c2]=size(L2);

MatchesError1 = zeros(c1,c2);
MatchesError2 = zeros(c1,c2);

GraphPairs=[];

c1a=(L1(2,comb1(:,1)).*L1(3,comb1(:,2))-L1(3,comb1(:,1)).*L1(2,comb1(:,2)));
c2a=(L1(3,comb1(:,1)).*L1(1,comb1(:,2))-L1(1,comb1(:,1)).*L1(3,comb1(:,2)));
c3a=(L1(1,comb1(:,1)).*L1(2,comb1(:,2))-L1(2,comb1(:,1)).*L1(1,comb1(:,2)));
c1a=c1a./c3a; c2a=c2a./c3a; c3a=c3a./c3a;

inPoint1=[c1a' c2a' c3a'];

c1b=(L2(2,comb2(:,1)).*L2(3,comb2(:,2))-L2(3,comb2(:,1)).*L2(2,comb2(:,2)));
c2b=(L2(3,comb2(:,1)).*L2(1,comb2(:,2))-L2(1,comb2(:,1)).*L2(3,comb2(:,2)));
c3b=(L2(1,comb2(:,1)).*L2(2,comb2(:,2))-L2(2,comb2(:,1)).*L2(1,comb2(:,2)));
c1b=c1b./c3b; c2b=c2b./c3b; c3b=c3b./c3b;

inPoint2=[c1b' c2b' c3b'];


for i=1:size(comb1,1)
    
    if isempty(H{i})
        continue;
    end
    
    
    x1=inPoint1(i,:)';
    Fx1  = F*x1;
    
    D_inPoint2_Fx1 = abs(inPoint2*Fx1)/norm(Fx1(1:2));
        
%     [Y2,Y22] = sort(D_inPoint2_Fx1);                 % sort by distance from current epipolar line (Fx1)
    
    id2=(D_inPoint2_Fx1<10);
    
    D2=D_inPoint2_Fx1(id2);
    comb22=comb2(id2,:);
    inPoint22=inPoint2(id2,:);

    [i size(comb1,1)]
%     size(comb22,1), pause
    
    corr=[];
    
    for j=1:size(comb22,1)
               
        x2=inPoint22(j,:)';

        Ftx2 = F'*x2;
         
        D1 = abs(x1'*Ftx2)/norm(Ftx2(1:2));
%         [D1 D2(j)]                  
%         x1a=[cross(Ftx2,L1(:,comb1(index1(i),1))) cross(Ftx2,L1(:,comb1(index1(i),2)))];
%         x1a=x1a./[repmat(x1a(3,:),3,1)];
% 
%         x2a=[cross(Fx1, L2(:,comb2(index2(j),1))) cross(Fx1, L2(:,comb2(index2(j),2)))];  
%         x2a=x2a./[repmat(x2a(3,:),3,1)];
%                         
%         % point to point distance
%         D1a=sqrt(sum((x1a-repmat(x1,1,2)).^2));
%         D1b=sqrt(sum((x2a-repmat(x2,1,2))).^2);
%                
%         D1=sum(D1a) + sum(D1b);
        
        
        % point to line distance
%         D2=abs(sum(Fx1.*x2,1))./sqrt(sum(Fx1(1:2).^2)) + abs(sum(Ftx2.*x1))./sqrt(sum(Ftx2(1:2).^2));

%                 [D1 D2(j)]
        if (all([D1 D2(j)]<10))
    
%             [D1 D2(j)]

            angL1EL1 = acosd(dot(L1(1:2,comb1(i,1)),Ftx2(1:2))/(norm(L1(1:2,comb1(i,1)))*norm(Ftx2(1:2))));
            angL1EL2 = acosd(dot(L1(1:2,comb1(i,2)),Ftx2(1:2))/(norm(L1(1:2,comb1(i,2)))*norm(Ftx2(1:2))));

            angL1EL1=min(angL1EL1,abs(180-angL1EL1));
            angL1EL2=min(angL1EL2,abs(180-angL1EL2));

            angL2EL1 = acosd(dot(L2(1:2,comb22(j,1)),Fx1(1:2))/(norm(L2(1:2,comb22(j,1)))*norm(Fx1(1:2))));
            angL2EL2 = acosd(dot(L2(1:2,comb22(j,2)),Fx1(1:2))/(norm(L2(1:2,comb22(j,2)))*norm(Fx1(1:2))));

            angL2EL1=min(angL2EL1,abs(180-angL2EL1));
            angL2EL2=min(angL2EL2,abs(180-angL2EL2));

            angLEL=[angL1EL1 angL1EL2; angL2EL1 angL2EL2];
            AngOK=any(all([angLEL(:,1)>angLEL(:,2) angLEL(:,1)<angLEL(:,2)],1));

            if ~AngOK
                comb22(j,:)=comb22(j,[2 1]);
            end



            endPoint1=reshape(p1(:,comb1(i,:)),3,4);
            endPoint2=reshape(p2(:,comb22(j,:)),3,4);

            epL1=F'*endPoint2;
            epL2=F*endPoint1;
            

            
            [Ang1 , ~ , ~ , cpp1 , cqq1]=CheckAngleAndIntersection2SegLines(p1(:,comb1(i,1)),p1(:,comb1(i,2)),I1,1);
            [Ang2 , ~ , ~ , cpp2 , cqq2]=CheckAngleAndIntersection2SegLines(p2(:,comb22(j,1)),p2(:,comb22(j,2)),I2,1);
            
            sr1I1=sqrt(sum((cpp1(1:3)-cpp1(7:9)).^2))*sqrt(sum((cpp1(4:6)-cpp1(7:9)).^2))/sum((cpp1(4:6)-cpp1(7:9)).^2);
            sr2I1=sqrt(sum((cqq1(1:3)-cqq1(7:9)).^2))*sqrt(sum((cqq1(4:6)-cqq1(7:9)).^2))/sum((cqq1(4:6)-cqq1(7:9)).^2);
            
%             sr1I1b=sqrt(sum((cpp1(1:3)-cpp1(4:6)).^2))*sqrt(sum((cpp1(4:6)-cpp1(7:9)).^2))/sum((cpp1(4:6)-cpp1(7:9)).^2);
%             sr2I1b=sqrt(sum((cqq1(1:3)-cqq1(4:6)).^2))*sqrt(sum((cqq1(4:6)-cqq1(7:9)).^2))/sum((cqq1(4:6)-cqq1(7:9)).^2);
            
            sr1I2=sqrt(sum((cpp2(1:3)-cpp2(7:9)).^2))*sqrt(sum((cpp2(4:6)-cpp2(7:9)).^2))/sum((cpp2(4:6)-cpp2(7:9)).^2);
            sr2I2=sqrt(sum((cqq2(1:3)-cqq2(7:9)).^2))*sqrt(sum((cqq2(4:6)-cqq2(7:9)).^2))/sum((cqq2(4:6)-cqq2(7:9)).^2);
            

      
            [int1ab]=check_sidedness_1(p1(:,comb1(i,1)),epL1(:,1:2),'P2L');
            [int1cd]=check_sidedness_1(p1(:,comb1(i,2)),epL1(:,3:4),'P2L');
            

             
%             over1=[int1ab int1cd]
            
            if (all([int1ab int1cd]))
                
                [int1ab]=check_sidedness_1(p2(:,comb22(j,1)),epL2(:,1:2),'P2L');
                [int1cd]=check_sidedness_1(p2(:,comb22(j,2)),epL2(:,3:4),'P2L');
 
                
%                 over2=[int1ab int1cd]
                        if ~(all([int1ab int1cd]))
                            
                            continue;
                        end
                        

                [angI1a]=CheckAngleAndIntersection2SegLines(cpp1(1:6),Ftx2,I1,1);
                [angI1b]=CheckAngleAndIntersection2SegLines(cqq1(1:6),Ftx2,I1,1);
                
                [angI2a]=CheckAngleAndIntersection2SegLines(cpp2(1:6),Fx1*-1,I1,1);
                [angI2b]=CheckAngleAndIntersection2SegLines(cqq2(1:6),Fx1*-1,I1,1);
                
%                 [angI1a angI1b D1;  angI2a angI2b D2(j)]
                
                angI1a=mod(angI1a,360); 
                angI1b=mod(angI1b,360);
                
                angI2a=mod(angI2a,360);
                angI2b=mod(angI2b,360);
                
                anglesI1I2=[angI1a angI1b;  angI2a angI2b];
%                 COSANGratio=[cosd(angI1a/angI1b)    cosd(angI2a/angI2b)   angI1a/angI1b    angI2a/angI2b]
%                 SINANGratio=[sind(angI1a/angI1b)    sind(angI2a/angI2b)   angI1a/angI1b    angI2a/angI2b]
%                 [(anglesI1I2(1)-180)*(anglesI1I2(2)-180)  (anglesI1I2(3)-180)*(anglesI1I2(4)-180)]

                angratio1=[anglesI1I2(:,2) - anglesI1I2(:,1)]';
                fullrot=angratio1<0;
                
                angratio1=angratio1+360.*fullrot;
                
                angratio1=min(abs(angratio1))/max(abs(angratio1))*100;
                      
%                 angratio2=abs(anglesI1I2(:,1)-anglesI1I2(:,2))
%                 angratio2=min(angratio2)/max(angratio2)*100
                
                sr=[sr1I1 sr2I1 ; sr1I2 sr2I2];
%                 diffSR=[sr(1,:)-sr(2,:)];
%                 sumSR=[sum(abs(sr(1,:)-sr(2,:)))]
%                 [sr1I1/sr1I2 sr2I1/sr2I2; sr1I1b/sr1I2b sr2I1b/sr2I2b]
%                 [min(dc1,dc2)/max(dc1,dc2)*100]
                
                if ~all([abs(sr(1,:)-sr(2,:))<10 ])    %if ~all([abs(sr(1,:)-sr(2,:))<0.4 angratio1>80])
                    continue;
                end
                    
               
            else
                continue;
            end
            
            Hx1=H{i}*x1;    Hx1=Hx1/Hx1(end);
            Htx2=H{i}\x2;   Htx2=Htx2/Htx2(end);
            
            SymErr=max([norm(x1-Htx2) norm(x2-Hx1)]);
            
            if SymErr>15
                continue;
            end
            
            % estimate lines error
%             H_line=inv(H{i}');
            L1_test=(H{i}')*L2(:,comb22(j,:));
            L2_test=(H{i}')\L1(:,comb1(i,:));

            d1s=abs(sum(L1_test.*p1(1:3,comb1(i,:)),1))./sqrt(sum(L1_test(1:2,:).^2,1));
            d1e=abs(sum(L1_test.*p1(4:6,comb1(i,:)),1))./sqrt(sum(L1_test(1:2,:).^2,1));


            d2s=abs(sum(L2_test.*p2(1:3,comb22(j,:)),1))./sqrt(sum(L2_test(1:2,:).^2,1)); 
            d2e=abs(sum(L2_test.*p2(4:6,comb22(j,:)),1))./sqrt(sum(L2_test(1:2,:).^2,1)); 


            LineErr=max([d1s d1e d2s d2e]);
            
            if (LineErr>15)
                continue;
            end
            % end-of lines error

                
%             ddd3=732887328
            corr=cat(1,corr,[comb1(i,:) comb22(j,:) D1  D2(j) SymErr LineErr]);
        
% 
            if showFig~=0
                figure(curFig),  subplot(1,2,1), imshow(I1), hold on,
                            plot(x1(1),x1(2),'ro',x1(1),x1(2),'r+','MarkerSize',10);
%                             plot(x1a(1),x1a(2),'co',x1a(1),x1a(2),'c+','MarkerSize',10);
%                             plot(x1a(4),x1a(5),'bo',x1a(4),x1a(5),'b+','MarkerSize',10);
                            hline(L1(:,comb1(i,1)),'r');
                            hline(L1(:,comb1(i,2)),'y');

                            hline(Ftx2(:,1),'w');

                            hline(epL1(:,1),'g');
                            hline(epL1(:,2),'g');
                            hline(epL1(:,3),'c');
                            hline(epL1(:,4),'c');
% 
%                             hline(e1e2,'b');

                            plot([p1(1,comb1(i,1)); p1(4,comb1(i,1))],[p1(2,comb1(i,1)); p1(5,comb1(i,1))],'LineWidth',3);
                            plot([p1(1,comb1(i,2)); p1(4,comb1(i,2))],[p1(2,comb1(i,2)); p1(5,comb1(i,2))],'LineWidth',3);
                            hold off;

                            subplot(1,2,2), imshow(I2), hold on,
                            plot(x2(1),x2(2),'ro',x2(1),x2(2),'r+','MarkerSize',10);
%                             plot(x2a(1),x2a(2),'co',x2a(1),x2a(2),'c+','MarkerSize',10);
%                             plot(x2a(4),x2a(5),'bo',x2a(4),x2a(5),'b+','MarkerSize',10);
                            hline(L2(:,comb22(j,1)),'r');
                            hline(L2(:,comb22(j,2)),'y');

                            hline(Fx1(:,1),'w');

                            hline(epL2(:,1),'g');
                            hline(epL2(:,2),'g');
                            hline(epL2(:,3),'c');
                            hline(epL2(:,4),'c');
% 
%                             hline(e1e2,'b');

                            plot([p2(1,comb22(j,1)); p2(4,comb22(j,1))],[p2(2,comb22(j,1)); p2(5,comb22(j,1))],'LineWidth',3);
                            plot([p2(1,comb22(j,2)); p2(4,comb22(j,2))],[p2(2,comb22(j,2)); p2(5,comb22(j,2))],'LineWidth',3);

                            hold off;
                            pause;
            end

        end
    end
    
    if isempty(corr)
        continue;
    end
    
    [m1 m2]=sort((corr(:,5)+corr(:,6)));
    
    if length(m2)>1
        corr=corr(m2(1:end),:);
    end
    
    for z=1:1
        Indx = sub2ind(size(MatchesError1), corr(z,1:2), corr(z,3:4));
        MatchesError1(Indx)= MatchesError1(Indx)+1;
%         MatchesError2(Indx)= MatchesError2(Indx)+exp(-D1/10)+exp(-sum(SymErr/2)/40);
        MatchesError2(Indx)= MatchesError2(Indx)+ exp(-(corr(z,5)))+exp(-(corr(z,6)))+exp(-(corr(z,7)))+exp(-(corr(z,8)));
    end
    
    GraphPairs = cat(1,GraphPairs,corr);

%     GraphPairs1(sub2ind(size(GraphPairs1),corr(:,1),corr(:,2)))=GraphPairs1(sub2ind(size(GraphPairs1),corr(:,1),corr(:,2)))+1;
    
    
    
end
                         
MatchError={MatchesError1;MatchesError2};

[matchVa1,BestMatches1]=sort(MatchesError1,2,'Descend');
[matchVa2,BestMatches2]=sort(MatchesError2,2,'Descend');

BestMatches{1,1}=[(1:size(L1,2))' BestMatches1];
BestMatches{2,1}=[(1:size(L1,2))' BestMatches2];

matchVa{1,1}=matchVa1;
matchVa{2,1}=matchVa2;


end



        

