function NSHTS1()
clc;
clear all;
warning off;
pop=100; 
gen=500;  
nshts(pop,gen)

function [number_of_objectives, number_of_decision_variables, min_range_of_decesion_variable, max_range_of_decesion_variable] = objective_description_function()
number_of_objectives = 2;
min_range_of_decesion_variable =ones(1, 10).*1;
max_range_of_decesion_variable =ones(1, 10).*41;
number_of_decision_variables = length(min_range_of_decesion_variable);

function f = evaluate_objective(x, M, V)
    f=p10bar(x(1:10));
% f = [];
% f(1) = 1 - exp(-4*x(1))*(sin(6*pi*x(1)))^6;
% sum = 0;
% for i = 2 : 6
%     sum = sum + x(i)/4;
% end
% g_x = 1 + 9*(sum)^(0.25);
% f(2) = g_x*(1 - ((f(1))/(g_x))^2);

function f = initialize_variables(N, M, V, min_range, max_range)
min = min_range;
max = max_range;
K = M + V;
for i = 1 : N
       for j = 1 : V
        f(i,j) = min(j) + (max(j) - min(j))*rand(1);
    end
      f(i,V + 1: K) = evaluate_objective(f(i,:), M, V);
end
function f = non_domination_sort_mod(x, M, V)
[N, m] = size(x);
clear m
front = 1;
F(front).f = [];
individual = [];
for i = 1 : N
    individual(i).n = 0;
    individual(i).p = [];
    for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M
            if (x(i,V + k) < x(j,V + k))
                dom_less = dom_less + 1;
            elseif (x(i,V + k) == x(j,V + k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= M
            individual(i).p = [individual(i).p j];
        end
    end   
    if individual(i).n == 0
        x(i,M + V + 1) = 1;
        F(front).f = [F(front).f i];
    end
end
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f)
       if ~isempty(individual(F(front).f(i)).p)
        	for j = 1 : length(individual(F(front).f(i)).p)
            	individual(individual(F(front).f(i)).p(j)).n = ...
                	individual(individual(F(front).f(i)).p(j)).n - 1;
        	   	if individual(individual(F(front).f(i)).p(j)).n == 0
               		x(individual(F(front).f(i)).p(j),M + V + 1) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + V + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end
current_index = 0;
for front = 1 : (length(F) - 1)
    distance = 0;
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);
    end
    current_index = current_index + i;
    sorted_based_on_objective = [];
    for i = 1 : M
        [sorted_based_on_objective, index_of_objectives] = ...
            sort(y(:,V + i));
        sorted_based_on_objective = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
        end
        f_max = ...
            sorted_based_on_objective(length(index_of_objectives), V + i);
        f_min = sorted_based_on_objective(1, V + i);
        y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)...
            = Inf;
        y(index_of_objectives(1),M + V + 1 + i) = Inf;
         for j = 2 : length(index_of_objectives) - 1
            next_obj  = sorted_based_on_objective(j + 1,V + i);
            previous_obj  = sorted_based_on_objective(j - 1,V + i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + V + 1 + i) = Inf;
            else
                y(index_of_objectives(j),M + V + 1 + i) = ...
                     (next_obj - previous_obj)/(f_max - f_min);
            end
         end
    end
    distance = [];
    distance(:,1) = zeros(length(F(front).f),1);
    for i = 1 : M
        distance(:,1) = distance(:,1) + y(:,M + V + 1 + i);
    end
    y(:,M + V + 2) = distance;
    y = y(:,1 : M + V + 2);
    z(previous_index:current_index,:) = y;
end
f = z();
function f  = hts_algorithm(chromosome, M, V, mu, mum, l_limit, u_limit)
[pop, variables] = size(chromosome);
rank = chromosome(:,variables - 1);
distance = variables;
[N,m] = size(chromosome);
clear m
p = 1;
for i = 1 : N
    child = [];
    child_1=[];
    Mean_pop=mean(chromosome);
    parent_1 = i;
    p2 = round(N*rand(1));
    if p2 < 1
        p2 = 1;
    end
    while isequal(chromosome(parent_1,:),chromosome(p2,:))
        p2 = round(N*rand(1));
        if p2 < 1
            p2 = 1;
        end
    end
    parent_1 = chromosome(parent_1,:);
    parent_2 = chromosome(p2,:);
    R=rand;
    for j = 1 : V
       if R<=0.3333  
            if rank(i)<rank(p2)
                if rand < 0.5
                    child(j)=(1-R^2)*parent_1(j);
                else
                    child(j)=(1-rand)*parent_1(j);
                end
            else
                if rand < 0.5
                    child(j)=(1-R^2)*parent_2(j);
                else
                    child(j)=(1-rand)*parent_2(j);
                end
            end
        elseif 0.3333< R && R <=0.6666 
            if rank(i)<rank(p2)
                if rand < 0.5              
                    child(j)=parent_1(j)+ (R)*(parent_1(j)-parent_2(j));
                else
                    child(j)=parent_1(j)+ (rand)*(parent_1(j)-parent_2(j));
                end
            else
                if rand < 0.5
                    child(j)=parent_1(j)+ (R)*(parent_2(j)-parent_1(j));
                else
                    child(j)=parent_1(j)+ (rand)*(parent_2(j)-parent_1(j));
                end
            end
        else 
            if rand < 0.5
               child(j)=parent_1(j)+ (rand)*(parent_1(j)-Mean_pop(j)*abs(R-rand));
            else
                child(j)=parent_1(j)+ (rand)*(parent_1(j)-Mean_pop(j)*round(1+rand));
            end
            
        end
        if child(j) > u_limit(j)
            child(j) = u_limit(j);
        elseif child(j) < l_limit(j)
            child(j) = l_limit(j);
        end  
    end
        child(:,V + 1: M + V) = evaluate_objective(child, M, V);
        child1(i,:) = child;
end
        f = child1;
 function nshts(pop,gen)
     no_of_runs=100;
     for ii=1:no_of_runs
format long;
[M, V, min_range, max_range] = objective_description_function();
chromosome = initialize_variables(pop, M, V, min_range, max_range);
chromosome = non_domination_sort_mod(chromosome, M, V);
for i = 1 : gen
    mu = 20;
    mum = 20;
    offspring_chromosome = hts_algorithm(chromosome,M, V, mu, mum, min_range, max_range);
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    clear temp
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = offspring_chromosome;
    intermediate_chromosome = non_domination_sort_mod(intermediate_chromosome, M, V);
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    if ~mod(i,20)
     fprintf('%d generations completed\n',i);
    end
end
Pareto_MOHTS{ii}=chromosome;
save('Pareto_MOHTS10bar.mat','Pareto_MOHTS')
     end
% if M == 2
%     plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
% elseif M ==3
%     plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
% end
 function f = tournament_selection(chromosome, pool_size, tour_size)
[pop, variables] = size(chromosome);
rank = variables - 1;
distance = variables;
for i = 1 : pool_size
     for j = 1 : tour_size
       candidate(j) = round(pop*rand(1));
       if candidate(j) == 0
            candidate(j) = 1;
        end
        if j > 1
            while ~isempty(find(candidate(1 : j - 1) == candidate(j)))
                candidate(j) = round(pop*rand(1));
                if candidate(j) == 0
                    candidate(j) = 1;
                end
            end
        end
    end
    for j = 1 : tour_size
        c_obj_rank(j) = chromosome(candidate(j),rank);
        c_obj_distance(j) = chromosome(candidate(j),distance);
    end
    min_candidate = ...
        find(c_obj_rank == min(c_obj_rank));
      if length(min_candidate) ~= 1
        max_candidate = ...
        find(c_obj_distance(min_candidate) == max(c_obj_distance(min_candidate)));
        if length(max_candidate) ~= 1
            max_candidate = max_candidate(1);
        end
        f(i,:) = chromosome(candidate(min_candidate(max_candidate)),:);
    else
           f(i,:) = chromosome(candidate(min_candidate(1)),:);
    end
end
   function f  = replace_chromosome(intermediate_chromosome, M, V,pop)
[N, m] = size(intermediate_chromosome);
[temp,index] = sort(intermediate_chromosome(:,M + V + 1));
clear temp m
for i = 1 : N
    sorted_chromosome(i,:) = intermediate_chromosome(index(i),:);
end
max_rank = max(intermediate_chromosome(:,M + V + 1));
previous_index = 0;
for i = 1 : max_rank
    current_index = max(find(sorted_chromosome(:,M + V + 1) == i));
    if current_index > pop
        remaining = pop - previous_index;
        temp_pop = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        [temp_sort,temp_sort_index] = ...
            sort(temp_pop(:, M + V + 2),'descend');
        for j = 1 : remaining
            f(previous_index + j,:) = temp_pop(temp_sort_index(j),:);
        end
        return;
    elseif current_index < pop
          f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
    else
         f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        return;
    end
     previous_index = current_index;
end