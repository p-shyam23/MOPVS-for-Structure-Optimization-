function [f] = p10bar(y)
% global panelty ;
% x;
% clc

% Counter of number of run
% load data.mat;
% kkk = kkk+1;
% save data.mat *kkk;

format shortg
% TLBO result

% y=[];
y=round(y);

Sec=[1.62 1.80 1.99 2.13 2.38 2.62 2.88 2.93 3.09 3.13 3.38 3.47 3.55 3.63...
    3.84 3.87 3.88 4.18 4.22 4.49 4.59 4.80 4.97 5.12 5.74 7.22 7.97 11.5...
    13.5 13.9 14.2 15.5 16.0 16.9 18.8 19.9 22.0 22.9 26.5 30.0 33.5];
x=Sec(y);
ro=0.1;
    [A,L,maxdisplacement,g1] = constrain_toporep1(x,ro); 
    g1(g1<0)=0;
%     f(1,1)= (sum(((A).*L))*ro); % obj fun 1 
%     f(2,1)=maxdisplacement; % obj fun 2 and 
    f(1,1)= (sum(((A).*L))*ro); % obj fun 1
    f(1,2)= maxdisplacement; % obj fun 2 and
    g=g1; % g1(=stress-25) is constraint violation allowable is 25 ksi
e1=3;
e2=3;
if sum(g) > 0
    f(1,1)=f(1,1)+f(1,1).*(1+e1*sum(g/25))^e2;
    f(1,2)=f(1,2)+f(1,2).*(1+e1*sum(g/25))^e2;
end



function [A,L,maxdisplacement,g1] = constrain_toporep1(x,ro)
A=x;

% E; modulus of elasticity
% A: area of cross section
% L: length of bar

% E=6.98e10; % N.m2
E=1e4; % Ksi

%EA=E*A;

% generation of coordinates and connectivities

nodeCoordinates=[360*2 360;360*2 0;360 360;360 0 ;0 360;0 0];% (inch)
elementNodes=[3 5;1 3;4 6;2 4;3 4;1 2;3 6;4 5;1 4;2 3];
numberElements=size(elementNodes,1);
numberNodes=size(nodeCoordinates,1);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

% for structure:
% displacements: displacement vector
% force : force vector
% stiffness: stiffness matrix

GDof=2*numberNodes; % GDof: total number of degrees of freedom
% displacements=zeros(GDof,1);
force=zeros(GDof,1);

force(4)=-100; % in Klb
force(8)=-100;

% computation of the system stiffness matrix
[stiffness,mass,L]=formStiffness2Dtruss(GDof,numberElements,...
    elementNodes,xx,yy,E,A,ro);
prescribedDof=(9:12);
[displacements,activeDof]=solution(GDof,prescribedDof,stiffness,force);
    
% stresses at elements
[sigma]=stresses2Dtruss(numberElements,elementNodes,xx,yy,displacements,E);


    for e=1:numberElements
        if sigma(e) < 0
            sigma(e)=-sigma(e);
        end
    end
    for e=1:GDof
        if displacements(e) < 0
            displacements(e)=-displacements(e);
        end
    end
    
%     c = [sigma'.*B'-25;displacements(3:8)-2];
% fit1= (sum(((A).*L))*ro)
ws=[25]; % Constraint Stress linite to 25 ksi

    for e=1:numberElements
        g1(e)=sigma(e)-ws(1);
    end
%     for e=1:length(displacements)
%         g2(e)=displacements(e)-ws(2);
%     end

    g1=[g1]; % constraint violation 
% Penalty function need to add

maxdisplacement=max(displacements); % Maximum displacement obj 2
% end

function [stiffness,mass,L]=formStiffness2Dtruss(GDof,numberElements,...
    elementNodes,xx,yy,E,A,ro)
stiffness=zeros(GDof);
mass=zeros(GDof);
% computation of the system stiffness matrix
for e=1:numberElements;
    
% elementDof: element degrees of freedom (Dof)
indice=elementNodes(e,:) ;
elementDof=[ indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
xa=xx(indice(2))-xx(indice(1));
ya=yy(indice(2))-yy(indice(1));
L(e)=sqrt(xa*xa+ya*ya);
C=xa/L(e);
S=ya/L(e);
k1=E*(A(e))/L(e)*[C*C C*S -C*C -C*S; C*S S*S -C*S -S*S;...
    -C*C -C*S C*C C*S;-C*S -S*S C*S S*S];
stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+k1;
end

function [displacements,activeDof] =solution(GDof,prescribedDof,stiffness,force)
% function to find solution in terms of global displacements
activeDof=setdiff([1:GDof]',[prescribedDof]);
U=stiffness(activeDof,activeDof)\force(activeDof);
displacements=zeros(GDof,1);
displacements(activeDof)=U;

function [sigma]=stresses2Dtruss(numberElements,elementNodes,xx,yy,displacements,E)
% stresses at elements
for e=1:numberElements
indice=elementNodes(e,:);
elementDof=[indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
xa=xx(indice(2))-xx(indice(1));
ya=yy(indice(2))-yy(indice(1));
L=sqrt(xa*xa+ya*ya);
C=xa/L;
S=ya/L;
sigma(e)=E/L*[-C -S C S]*displacements(elementDof);
end
