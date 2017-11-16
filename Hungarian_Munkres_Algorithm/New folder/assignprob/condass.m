function [k,C1,T1,C2,T2]=condass(A)
%CONDASS Calculate condition number of the assigment problem.
%
%[k,C1,T1,C2,T2]=condass(A)
%A  - A square cost matrix.
%k  - The condition number of the assigment problem.
%C1 - The best assigment.
%T1 - The lowest cost.
%C2 - The second best assignment.
%T2 - The second lowest cost.
%
%The condition number is calculated as the relative difference between
%the best and second best solutions, as described in Nystrom, Soderkvist,
%and Wedin, "A Note on some Identification Problems Arising in Roentgen
%Stereo Photogrammetric Analysis", J of Biomechanics, 27(10):1291-1294,
%1994.

% v1.0  96-09-14. Niclas Borlin, niclas@cs.umu.se.

% A substantial effort was put into this code. If you use it for a
% publication or otherwise, please include an acknowledgement and notify
% me by email. /Niclas

% Create a large number used to block selected assignments.
big=sum(sum(A))+1;

% Get best assigment.
[C1,T1]=hungarian(A);

% Initialize second best solution.
T2=inf;
C2=zeros(size(C1));

% Create a work matrix.
B=A;
for i=1:length(C1)
    % Block assigment in column i.
    B(C1(i),i)=big;
    % Get best assigment with this one blocked.
    [C,T]=hungarian(B);
    if (T<T2)
		% Remember it if it's the best so far.
		T2=T;
		C2=C;
    end
    % Remove blocking in column i.
    B(C1(i),i)=A(C1(i),i);
end

% Calculate difference...
mu=T2-T1;

% ...and condition number.
k=T1/mu;
