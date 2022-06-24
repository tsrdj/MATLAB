%---------------------------------
%---------------------------------
clc;clear all;
x=[5.935	14.2	2.265	2.27	2.91
1.523	13.1	0.597	0.75	2.62
2.599	12.7	1.237	1.11	1.72
4.009	15.2	1.649	0.81	3.02
4.687	14.7	2.312	2.50	2.22
8.044	15.6	3.641	4.51	2.36
2.766	13.3	1.244	1.03	1.97
6.538	17.0	2.618	2.39	1.85
6.451	12.9	3.147	5.52	2.01
3.314	12.2	1.606	2.18	1.82
3.777	13.0	2.119	2.83	1.80
1.530	13.8	0.798	0.84	4.25
2.768	13.6	1.336	1.75	2.64
6.585	14.9	2.763	1.91	3.17
]';
[p n]=size(x);
%-----------Manupulation-------------------------
s=cov(x');      R=(diag(diag(s))^-0.5)*s*(diag(diag(s))^-0.5);
[P L]=eig(R);   l=diag(L);
P=P(:,[1 2 5 4 3]);
l=l([1 2 5 4 3]);
TSV=sum(l);     pv=l/TSV;       cpv=cumsum(pv*100);
%-----------Output-----------------------------------
fprintf('\n\t Coefficients of %d  Principle Components (in column) are :\n',p);
disp(P);
fprintf('\n\t Variance of PC                :\n'); disp(l');
fprintf('\n\t proportion of Variance of PC  :\n'); disp(pv');
fprintf('\n\t Cumulative percent of Var. of PC:\n'); disp(cpv');
set(gcf,'color',[1 1 1]);   
plot(cpv,'-bx')