%   Raul Diaz
%   CECS 463 SOC II Spr 2015
%   Lecture Assignment 10 Due: 5/6/2015

%
clc; clear all; format compact; 

%% Problem 1 Sonar Detection 

%Load Files/ Set sampling freq & support vectors
load sonar1; Ts=0.01; Txtime=0.25; vel = 5000; 
n=0:length(x)-1; m=0:length(y)-1;
% Generate Impulse Response
h=fliplr(x);
s=conv(y,h); ns=0:length(s)-1;
%Determine Distance
[smax,nmax]=max(s);

d=vel*(nmax*Ts-Txtime)/2;
fprintf('\nIndex of Max: %.f \nDistance: %.f\n',nmax,d); 

%Plot Transmitted &  Recieved Signal 
figure (1); clf(1);
subplot(2,1,1); plot(n*Ts,x); title('Xmit Signal X'); grid on;
xlabel('t (sec)'); ylabel('x(t)');

subplot(2,1,2); plot(m*Ts,y); title('Received Signal Y'); grid on;
xlabel('t (sec)'); ylabel('y(t)');

%Plot Match Filter 
figure (2); clf(2);
plot(n*Ts,h); title('Matched Filter H'); grid on;
xlabel('t (sec)'); ylabel('h(t)');

%Plot Correlation Results
figure (3); clf(3);
plot(ns*Ts,s); title('Correlator Output'); grid on;
xlabel('t (sec)'); ylabel('s(t)');


%% Problem 2 Digital Message Reception 


format compact; clear all; clc;
load sonar2; Ts=0.01; xmitTime=0.5;
n=0:length(x)-1; m=0:length(y)-1;

figure (1); 
subplot(2,1,1); plot(n*Ts,x); title('Xmit Binary 1 Signal X'); grid on;
xlabel('t (sec)'); ylabel('x(t)');
subplot(2,1,2); plot(m*Ts,y); title('Received Signal Y'); grid on;
xlabel('t (sec)'); ylabel('y(t)');

h=fliplr(x);
s=conv(y,h); ns=0:length(s)-1; threshold=0.70*abs(max(s));
figure (2); 
plot(n*Ts,h); title('Matched Filter H'); grid on;
xlabel('t (sec)'); ylabel('h(t)');

figure(3); 
plot(ns*Ts,s); title('Correlator Output'); grid on;
xlabel('t (sec)'); ylabel('s(t)');

N=length(y)/length(x); ch=0; msg=[];
for i=1:N
    if (s(i*length(x))>=threshold)
        fprintf('1'); ch=ch+2^(7-mod(i-1,8));
    else
        fprintf('0');
    
    end
    if(mod(i,8)==0)
        fprintf(' ');    
        display(char(ch));
        msg=[msg,char(ch)];
        ch=0;
    end
end
fprintf('Message: %s\n\n',msg)


