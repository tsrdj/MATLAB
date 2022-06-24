clc;clear all;close all;
%-------------------------------------------------
% Canonical correlations for iris data
%-------------------------------------------------
load data_iris  % x is 4 x 150 with 3 groups  
s=cov(x');      % VCM(x)
xstd=diag(diag(s))^-0.5*x;  % x is standardized
[p n]=size(x);    q=2;  
%-------------manipulation--------------------
s11=s(1:q,1:q);         s12=s(1:q,q+1:p);
s21=s12';               s22=s(q+1:p,q+1:p);
A=s22^-0.5*s21*s11^-0.5;
R1=A'*A;                R2=A*A';
[E L1]=eig(R1);        [F L2]=eig(R2);
U=s11^-0.5*E;          V=s22^-0.5*F;
ro=sqrt(diag(L1));
%----------o\p---------------------------------
fprintf('\n\t  First canonical correlation coeff=%f',ro(2));
fprintf('\n\t  The coefficients of canonical variables are (in column)\n');
disp([U(:,2) V(:,2)]);

fprintf('\n\t  Second canonical correlation coeff=%f',ro(1));
fprintf('\n\t  The coefficients of canonical variables are (in column)\n');
disp([U(:,1) V(:,1)]);
%----------Verification---------------------------------
fprintf('\n\t  The V-C-M of canonical variables is\n');
ccmat=U'*s12*V;;        disp(ccmat)
%----------Scatter plot of first pair of Canonical variables---------------
l=U(:,2);           m=V(:,2);
x1=x(1:2,:);        x2=x(3:4,:);
u=l'*x1;            v=m'*x2;
fprintf('\n\t  The V-C-M of first pair of CV.\n');
disp(corrcoef(u,v));
c=[ones(1,50) 2*ones(1,50) 3*ones(1,50)];
set(gcf,'color',[1 1 1]);
scatter(u,v,5,c); 
xlabel('u1->');     ylabel('v1');
title('Scatter plot of first pair of Canonical variables')


