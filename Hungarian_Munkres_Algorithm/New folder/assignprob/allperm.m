function p=allperm(n)
%ALLPERM All permutation matrix.
%
%p=allperm(n)
%Returns a matrix with all permutations of 1:n stored row-wise.

% v1.0  95-07-18. Niclas Borlin, niclas@cs.umu.se.

if (n<=1)
	p=1;
else
	q=allperm(n-1);
	p=[];
	for i=1:n
		p=[p;i*ones(size(q,1),1) q+(q>=i)];
	end
end
