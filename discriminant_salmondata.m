clc; clear all; close all;
%-------------------Discriminant analysis for Salnom data--------------------------
load data_salmon;
 x1=x(:,1:50);
 x2=x(:,51:100);
 n1=50; n2=50;
%---------Sample data 1:----------------------------
% x1=[191	185	200	173	171	160	188	186	174	163
%     131	134	137	127	128	118	134	129	131	115
%     53	50	52	50	49	47	54	51	52	47];
% x2=[186	211	201	242	184	211	217	223	208	199
%     107	122	144	131	108	118	122	127	125	124
%     49	49	47	54	43	51	49	51	50	46];
% x=[x1 x2];
% n1=10;n2=10;
% x0=[194  124  49]';
%---------Sample data 2:----------------------------
% x1=[90	115.5	94.8	91.5	117	140.1	138	112.8	99	123	81	111
%     18.4	16.8	21.6	20.8	23.6	19.2	17.6	22.4	20	20.8	22	20];
% x2=[105	82.8	94.8	73.2	114	79.2	89.4	96	77.4	63	81	93
%     19.6	20.8	17.2	20.4	17.6	17.6	16	18.4	16.4	18.8	14	14.8];
x=[x1 x2];
%n1=12; n2=12;
x0=[85,18.8]';
%------------------Manipulation----------------------------------------------------
m1=mean(x1,2);
m2=mean(x2,2);
s1=cov(x1');
s2=cov(x2');
A=(n1-1)*s1+(n2-1)*s2;
s=A/(n1+n2-2);
%-----------------FDF--------------------------------------------------------------
d=inv(s)*(m1-m2);                              %since for sample,so d instead of del.
c=(m1-m2)'*inv(s)*(m1+m2)/2;
for i=1:n1+n2
   if( d'*x(:,i)>c)
       pop(i)=1;
   else
       pop(i)=2;
       
   end
end
n2_1=sum(pop(1:n1)==2);
n1_2=sum(pop(n1+1:n1+n2)==1);
n=n2_1+n1_2;
p2_1=n2_1/n1;
p1_2=n1_2/n2;
fprintf('\n\t coefficint of FDF is:');
disp(d');
fprintf('\n\t constantc=%f\n',c);
fprintf('\n\t n(2/1)=%f\n',n2_1);
fprintf('\n\t n(1/2)=%f\n',n1_2);
fprintf('\n\t total missclassified =%f\n',n);
fprintf('\n\t p(2/1)=%f\n',p2_1);
fprintf('\n\t p(1/2)=%f\n',p1_2);
fprintf('\n\t missclassification prob=%f\n',n/(n1+n2));
%--Classification of x0---------------------------------------------------
if( d'*x0>c)
    fprintf('\n\t Given individual belogs to Population 1');
else
    fprintf('\n\t Given individual belogs to Population 2');
end





    
    