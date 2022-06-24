clc;clear all;close all;
%-------------------------------------------------
% Canonical correlations based on V-C-M
%-------------------------------------------------
s=[8    2   3   1
    2   5   -1  3
    3   -1  6   -2
    1   3   -2  7];
p=4;    q=2;    
%-------------------------------------------------
% Canonical correlations based on V-C-M
%-------------------------------------------------
% p=6;    q=3; 
% s=[12 4 1 2 5 1
% 4 10 1 1 2 1
% 1 1 5 2 2 2
% 2 1 2 7 2 3
% 5 2 2 2 8 1
% 1 1 2 3 1 5];
%-------------manipulation--------------------
s11=s(1:q,1:q);         s12=s(1:q,q+1:p);
s21=s12';               s22=s(q+1:p,q+1:p);
R1=s11^-0.5*s12*inv(s22)*s21*s11^-0.5;
R2=s22^-0.5*s21*inv(s11)*s12*s22^-0.5;
[P1 L1]=eig(R1);        [P2 L2]=eig(R2);
U=s11^-0.5*P1;          V=s22^-0.5*P2;
ro=sqrt(diag(L1));
%----------o\p---------------------------------
fprintf('\n\t  First canonical correlation coeff=%f',ro(1));
fprintf('\n\t  The coefficients of canonical variables are (in column)\n');
disp([U(:,1) V(:,1)]);

fprintf('\n\t  Second canonical correlation coeff=%f',ro(2));
fprintf('\n\t  The coefficients of canonical variables are (in column)\n');
disp([U(:,2) V(:,2)]);
%----------Verification---------------------------------
ccmat=U'*s12*V;
fprintf('\n\t  The V-C-M of canonical variables is\n');
disp(ccmat)
