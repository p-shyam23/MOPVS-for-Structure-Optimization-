%% Cited from NSGA-II All rights reserved.
function f = evaluate_objective(Fobj,x)

%% function f = evaluate_objective(x)
% Function to evaluate the objective functions for the given input vector
% x. x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% D is the number of decision variables. 
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and D matches your initial user
% input.
% A set of testing function is stored in folder TEST
%% Retrieve function from folder TEST
% % for i=1:size(x,1)
% %     [obj,const]=feval(Fobj,x(i,:));
% %     f(i,:) = fpenal01(obj,const,0); % change the name of function we can use different functions
% % end

[obj,const]=feval(Fobj,x);
const_vio=sum(const);
c0=25; %stress limit

f = fpenal01(obj,const_vio,25);

function fp=fpenal01(f,g,c0)
e1=3;
e2=3;
if g > 0
    fp=f+f.*(1+e1*g/c0)^e2;
else
    fp=f;
end