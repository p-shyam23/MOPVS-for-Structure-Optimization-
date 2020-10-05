function f = MOSOS02(Fobj,D,M,LB,UB,Pop,ecosize,GEN,ishow,narchive)
%% Multiple Objective Symbiotic Organisms Search (MOSOS)
% MOSOS is developed by Doddy Prayogo and Duc-Hoc Tran, MOSOS is based on
% the SOS algorithm which is one of the most recently defined algorithms
% f - optimal fitness
% X - optimal solution
% D  Dimensional of problem at hand
% M Number of objective function
% eco is ecosystem which is a matrix consists of all individuals
% ecosize is number of individual in ecosystem
% LB lower boundary constraint
% UB upper boundary constraint
%% Algorithm Variables
K = D+M;
eco = Pop(:,1:K+1);
eco_ad = zeros(ecosize,K);
%% Optimization Circle
g = 1;
fout=eco(:,1:K);
while g<=GEN % for each generation    
    for i = 1:ecosize   % Organisms' Looping (for each individual)
        %% Mutualism Phase
        % Choose organism j randomly other than organism i
        j = floor(rand()* ecosize) + 1;
        while j==i
            j = floor(rand()* ecosize) + 1;
        end
        % Determine Mutual Vector & Beneficial Factor
        mutualVector=mean([eco(i,1:D);eco(j,1:D)]);
        BF1=round(1+rand); BF2=round(1+rand);        
        % randomly select the best organism from first non-dominated front of eco
        Ffront = eco((find(eco(:,K+1)==1)),:); % first front
        ri =  floor(size(Ffront,1)*rand())+1; % ri is random index
        bestOrganism = eco(ri,1:D);
%         fitness_i = eco(i,D+1);
%         fitness_j = eco(j,D+2);
%         bestFitness_i= eco(ri,D+1);
%         bestFitness_j= eco(ri,D+2);
%         BF1=fitness_i/bestFitness_i;  % Adaptive BF1 
%         if BF1>2
%             BF1=2;
%         elseif BF1<1
%             BF1=1;
%         end
%         BF2=fitness_j/bestFitness_j; % Adaptive BF2
%         if BF2>2
%             BF2=2; 
%         elseif BF2<1
%             BF2=1;
%         end

        % Calculate new solution after Mutualism Phase
        ecoNew1 = eco(i,1:D)+rand(1,D).*(bestOrganism-BF1.*mutualVector);
        ecoNew2 = eco(j,1:D)+rand(1,D).*(bestOrganism-BF2.*mutualVector);
        % Handling constraints
        ecoNew1 = bound(ecoNew1(:,1:D),UB,LB); 
        ecoNew2 = bound(ecoNew2(:,1:D),UB,LB);
        % Evaluate function at ecoNew1 and ecoNew2
        ecoNew1(:,D + 1: K) = evaluate_objective(Fobj,ecoNew1(:,1:D));
        ecoNew2(:,D + 1: K) = evaluate_objective(Fobj,ecoNew2(:,1:D));
        % Nondomination checking of 2 trial individuals
        % For the first trail ecoNew1
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (ecoNew1(:,D+k)<eco(i,D+k))
                dom_less = dom_less + 1;
            elseif (ecoNew1(:,D+k)== eco(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector ecoNew1 dominates
            % target vector Xi. Replace Xi by ecoNew1 in current ecosystem and
            % add Xi to advanced population 2
            eco_ad(i,1:K) = eco(i,1:K); % add Xi to advanced eco
            eco(i,1:K) = ecoNew1(:,1:K); % replace Xi by ecoNew1            
        else % else Add Xi (trial vector) to advanced eco
            eco_ad(i,1:K)= ecoNew1;
        end % end if
        % For the second trail ecoNew2
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (ecoNew2(:,D+k)<eco(j,D+k))
                dom_less = dom_less + 1;
            elseif (ecoNew2(:,D+k)== eco(j,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector ecoNew1 dominates
            % target vector Xi. Replace Xj by ecoNew2 in current ecosystem and
            % add Xj to advanced population
            eco_ad(j,1:K) = eco(j,1:K); % add Xi to advanced eco
            eco(j,1:K) = ecoNew2(:,1:K); % replace Xi by ecoNew1            
        else % else Add Xi (trial vector) to advanced eco
            eco_ad(j,1:K)= ecoNew2;
        end % end if
        %% Commensialism Phase
        % Choose organism j randomly other than organism i
        j = floor(rand()* ecosize) + 1;
        while j==i
            j = floor(rand()* ecosize) + 1;
        end
        % Calculate new solution after Commensalism Phase
        ecoNew1 = eco(i,1:D)+(rand(1,D)*2-1).*(bestOrganism-eco(j,1:D));
        ecoNew1 = bound(ecoNew1(:,1:D),UB,LB);
        % Evaluate function at ecoNew1 and ecoNew2
        ecoNew1(:,D + 1: K) = evaluate_objective(Fobj,ecoNew1(:,1:D));
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (ecoNew1(:,D+k)<eco(i,D+k))
                dom_less = dom_less + 1;
            elseif (ecoNew1(:,D+k)== eco(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector ecoNew1 dominates
            % target vector Xi. Replace Xi by ecoNew1 in current ecosystem and
            % add Xi to advanced population
            eco_ad(i,1:K) = eco(i,1:K); % add Xi to advanced eco
            eco(i,1:K) = ecoNew1(:,1:K); % replace Xi by ecoNew1            
        else % else Add Xi (trial vector) to advanced eco
            eco_ad(i,1:K)= ecoNew1;
        end % end if
        %% Parasitism Phase
        % Choose organism j randomly other than organism i
        j = floor(rand()* ecosize) + 1;
        while j==i
            j = floor(rand()* ecosize) + 1;
        end
        % Determine Parasite Vector
        parasiteVector=eco(i,1:D);
        seed=randperm(D);
        pick=seed(1:ceil(rand*D));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(UB(pick)-LB(pick))+LB(pick);        
        % Evaluate the Parasite Vector
        parasiteVector(:,D + 1: K) = evaluate_objective(Fobj,parasiteVector(:,1:D));
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (parasiteVector(:,D+k)<eco(j,D+k))
                dom_less = dom_less + 1;
            elseif (parasiteVector(:,D+k)== eco(j,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector ecoNew1 dominates
            % target vector Xi. Replace Xi by ecoNew1 in current ecosystem and
            % add Xi to advanced population
            eco_ad(j,1:K) = eco(j,1:K); % add Xi to advanced eco
            eco(j,1:K) = parasiteVector(:,1:K); % replace Xi by ecoNew1            
        else % else Add Xi (trial vector) to advanced eco
            eco_ad(j,1:K)= parasiteVector;
        end % end if 
    end % end for i
    if rem(g, ishow) == 0
        fprintf('Generation: %d\n', g);        
    end
    eco_com = [eco(:,1:K) ; eco_ad];
    fout=pareto_selection(eco_com,fout,narchive,D,K);
    intermediate_eco = non_domination_sort_mod(eco_com, M, D);
    Pop  = replace_chromosome(intermediate_eco, M,D,ecosize);
    eco=Pop(:,1:K+1); %
    g = g+1;
end % end g
f= fout;
%% Subfunctions
% Check the boundary limit
function a=bound(a,ub,lb)
a(a>ub)=ub(a>ub); a(a<lb)=lb(a<lb);

function fout=pareto_selection(eco_com,fout0,narchive,D,K);

fpareto=fout0(:,D+1:K);
xpareto=fout0(:,1:D);
x1=eco_com(:,1:D);
f1=eco_com(:,D+1:K);

x=[xpareto' x1'];
f=[fpareto' f1'];
[m1,n1]=size(x);

for i=1:n1
    xi=x(:,i);
    fi=f(:,i);
    gi=0;
    A(i,i)=0;
    for j=i+1:n1
        xj=x(:,j);
        fj=f(:,j);
        gj=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
        A(i,j)=p_count1;
        A(j,i)=p_count2;
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=sum(A,1);
Indm=[];
for i=1:n1
    if B(i)==0
        Indm=[Indm i];
    end
end
nndm=length(Indm);
pareto1=x(:,Indm);
fpareto1=f(:,Indm);

if nndm > narchive
    nsl=farchive22(fpareto1,narchive);
    pareto1=pareto1(:,nsl);
    fpareto1=fpareto1(:,nsl);
end
fout=[pareto1' fpareto1'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p1,p2]=fdominated(f1,g1,f2,g2)
n=length(f1);
mg1=max(g1);
mg2=max(g2);

icount11=0;
icount12=0;
icount21=0;
icount22=0;

if mg1<=0&mg2<=0
    for i=1:n
        if f1(i) <= f2(i)
            icount11=icount11+1;
        end
        if f1(i) < f2(i)
            icount12=icount12+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if f2(i) <= f1(i)
            icount21=icount21+1;
        end
        if f2(i) < f1(i)
            icount22=icount22+1;
        end
    end
    if icount11 == n & icount12 > 0
        p1=1;
    else
        p1=0;
    end
    if icount21 == n & icount22 > 0
        p2=1;
    else
        p2=0;
    end
elseif mg1 <=0 & mg2 > 0
    p1=1;p2=0;
elseif mg2 <=0 & mg1 > 0
    p1=0;p2=1;
else
    if mg1 <= mg2
        p1=1;p2=0;
    else
        p1=0;p2=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iselect=farchive22(f,narchive)

[m,n]=size(f);

fmax=max(f,[],2);fmin=min(f,[],2);
fdel=max(fmax-fmin,1e-5);
for i=1:n
    f(:,i)=(f(:,i)-fmin)./fdel;
end
nn=narchive;
x1=linspace(0,1,nn);
x2=1-x1;
iselect0=1:n;
for i=1:narchive
    df=[f(1,:)-x1(i);f(2,:)-x2(i)];
    d2=std(df);
    [d2min,imin]=min(d2);
    iselect(i)=iselect0(imin);
    iselect0(imin)=[];
    f(:,imin)=[];
end

%%%%%%%%%%%%%%%%%%%%%%% End of file %%%%%%%%%%%%%%%%%%%%%%%%%%%
