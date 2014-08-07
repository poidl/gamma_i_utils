close all

% number of grid points X*Y*Z;
x=[6*3,12*6,18*9,24*12,30*15,36*18]*1e4; 


% execution time for gamma^i
tg=[0.2,1.,2.5,4.3,6.7,9.8];

% execution time for omega
to=[4.5,20.2,47,90,154,226];

gg=polyfit(x,tg,1);
oo=polyfit(x,to,1);

plot(x,tg)
hold on
plot(x,tg,'*')



figure
plot(x,to)
hold on
plot(x,to,'*')


