function Hv=fhypervolume(fReal,W)
% output file name
f=fReal;% Extract the results from the last iteration
[f1min,ns]=sort(f(1,:));% Sort the first objective value in ascending order
f=f(:,ns);
Idl=[];% Collect Pareto points that exceed the RF point
for i=1:size(f,2)
    if f(1,i)>W(1)
        Idl=[Idl i];
    elseif f(2,i)>W(2)
        Idl=[Idl i];
    end
end
f(:,Idl)=[];% delete all the points in Idl
for i=1:size(f,2)
    if i==1
        V(i)=abs(W(1)-f(1,i))*abs(W(2)-f(2,i));
    else
        V(i)=abs(W(1)-f(1,i))*abs(f(2,i-1)-f(2,i));
    end
    
end
if size(f,2)==0
    Hv=0;
else
    Hv=sum(V);
end

