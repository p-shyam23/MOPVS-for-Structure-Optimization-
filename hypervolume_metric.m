function HV= hypervolume_metric( A)
% This function to calculate the volume ( in the objective space) covered
% by members of a non-dominated set of solution A
%   HV= Pi( vi)
% vi is constructed with a reference point W and the soulution Xi as the
% digonal corners of the hypercube. The reference point can be found by
% constructing a vector of worst objective function values
% Supposed we are considering the minimization problem. It means that all
% objective functions are minimized
[m,n]=size(A); % m is number of solutions in set A
% n is number of objective function;
B=sortrows(A,n); % sort
p=B(m,:); % the worst point
d=zeros(1,m); % distance for vi to worst point p
v=zeros(1,m); % volume values of vi
for i=1:m
    d(i)=sqrt(sum((B(i,:)-p).^2));
    if d(i)==0
        v(i)=1;
    else
        v(i)=d(i);
    end    
end % end for
HV= sum(v); % hypervolume value
end