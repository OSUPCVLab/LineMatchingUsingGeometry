A=magic(10);
B=A(4:7,4:7);
disp('Cost matrix:')
disp(B)
disp('Calculating best assignment...');
[c,t]=hungarian(B);
disp('Best assignment (as row indices):')
disp(c)
disp('Best assignment (as logical matrix):')
disp(logical(full(sparse(c,1:4,1))))
disp('Lowest cost:')
disp(t)

disp(sprintf('\nCalculating condition number for solution...'));
[k,c1,t1,c2,t2]=condass(B);
disp('Lowest cost (should be same as above): ')
disp(t1)
disp('corresponding assignment (should be same as above):')
disp(c1)
disp('Second lowest cost: ')
disp(t2)
disp('corresponding assignment:')
disp(c2)
disp('Condition number for solution:')
disp(k)

disp(sprintf('\nCalculating all possible costs...'));
[c,p]=allcosts(B);
% Sort by cost.
[y,i]=sort(c);
disp('The three lowest costs:')
disp(c(i(1:3)))
disp('Corresponding assignments:')
disp(p(i(1:3),:))
disp('The three highest costs:')
disp(c(i(end+[-2:0])))
disp('Corresponding assignments:')
disp(p(i(end+[-2:0]),:))
