function [c,p]=allcosts(C)
%ALLCOSTS Calculate all costs for an assignment problem.
%
%[c,p]=allcosts(C)
%c returns the costs, p the corresponding permutations.

% v1.0  95-07-18. Niclas Borlin, niclas@cs.umu.se.

p=allperm(size(C,1));

c=zeros(size(p,1),1);

I=eye(size(C,1));

for i=1:size(p,1)
	c(i)=sum(C(logical(sparse(p(i,:),1:size(C,1),1))));
end
