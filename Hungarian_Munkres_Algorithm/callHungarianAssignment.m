%% call optimization using Hungarian code (assignmentoptimal)
%[assignment, cost] = assignmentoptimal(distMatrix)
gr=1;  % 1 or 2  (based on cost affinity marix addembled , see main1.m)
per=0.09;
s=max(matchVa{gr}(:,1));

j=matchVa{gr}(:,1)>per*s;
% j=all([matchVa{gr}(:,1)>per*s (matchVa{gr}(:,1)-matchVa{gr}(:,2))>per*s/2],2);
nnz(j)

B=matchVa{gr};
for i=1:size(B,1)
    B(i,BestMatches{gr}(i,2:end))=B(i,:);
end;

 BB=1./B;   % 0 values assigned inf
 
 BB2=BB(j,:);
 
[assignment, cost] = assignmentoptimal(BB2);
% [assignment, cost] = assignmentsuboptimal1(BB2);


matchHung=[BestMatches{gr}(j,1) assignment];

ok=all([find(j) assignment>0],2);

[nnz(assignment) nnz(ok)]

PlotCorrespLines2(I1,I2,p1,p2,L1,L2,matchHung(ok,:))