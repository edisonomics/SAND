% test the even for perfect sin curve fft doesn't produce perfect peak
% because of finite sampling size even though there is no decay
t=1:0.1:10;
ntime=length(t);
y=2*sin(2*pi*t);%.*exp(-1*t);
plot(t,y);
fx=(0:(ntime/2))/ntime;
fy=fft(y);
P2=abs(fy/ntime);
P1=P2(1:ntime/2+1);
P1(2:end-1)=2*P1(2:end-1);
plot(fx,P1);
