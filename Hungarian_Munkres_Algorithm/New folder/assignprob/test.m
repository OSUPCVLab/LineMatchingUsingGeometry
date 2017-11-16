disp('Testing hungarian...');
A=magic(10);
B=A(4:7,4:7);
[c,t]=hungarian(B);
if (any(c~=[2 1 3 4]))
	disp('Wrong coupling!');
elseif (t~=77)
	disp('Wrong solution!');
else
	disp('Hungarian appears OK.');
end

disp('Testing condass...');
[k,c1,t1,c2,t2]=condass(B);
if (any(c1~=[2 1 3 4]))
	disp('Wrong best coupling!');
elseif (t1~=77)
	disp('Wrong lowest cost!');
elseif (any(c2~=[1 2 3 4]))
	disp('Wrong second best coupling!');
elseif (t2~=102)
	disp('Wrong second lowest cost!');
elseif (k~=t1/(t2-t1))
	disp('Wrong condition number!');
else
	disp('condass appears OK.');
end

disp('Testing allcosts...');
[c,p]=allcosts(B);
cTrue=[	102
		127
		202
		227
		222
		222
		77
		102
		182
		202
		202
		197
		177
		202
		182
		202
		302
		297
		202
		202
		207
		202
		307
		302
];
pTrue=[ 1     2     3     4
		1     2     4     3
		1     3     2     4
		1     3     4     2
		1     4     2     3
		1     4     3     2
		2     1     3     4
		2     1     4     3
		2     3     1     4
		2     3     4     1
		2     4     1     3
		2     4     3     1
		3     1     2     4
		3     1     4     2
		3     2     1     4
		3     2     4     1
		3     4     1     2
		3     4     2     1
		4     1     2     3
		4     1     3     2
		4     2     1     3
		4     2     3     1
		4     3     1     2
		4     3     2     1
	  ];
if (any(c~=cTrue))
	disp('Wrong costs!');
elseif (any(p~=pTrue))
	disp('Wrong permutations!');
else	
	disp('allcosts appears OK.');
end
