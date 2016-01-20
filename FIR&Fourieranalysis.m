
% Raul Diaz
% CECS 463 SOC II SPRING 2015
% Lab #6  due Due: 5/6/2015
clc; clear all; clf; format compact;

%% Problem 1
clear all; figure(1); clf(1); hold on;
fprintf('\n\n\nSTART ASSIGN #10');
fprintf('\nProblem 1\n');
f1=250; f2=600; fs=3000;    %Frequencies
T1=1/f1; T2=1/f2; Ts=1/fs; t=[0:Ts:4*max([T1,T2])];
x=cos(2*pi*f1*t)+ cos(2*pi*f2*t);
plot(t,x,'b--');grid;title('Plot of sinusoids at 250 and 600 Hz');
xlabel('Time (sec)'); ylabel('Signal');
pause(1.0);

r=exp(1j*2*pi*f2/fs);       %Zero location (conjugate pairs)
b=poly([r,conj(r)]);        %Denominator polynomial with real coefficients
rb=roots(b);                %Show roots

Hw=abs(sum(b.*exp(-1j*2*pi*f1/fs).^[0:2])); %Gain at f2
N=length(x);                %%Implement FIR filter
    for i=1:N
        if(i>2) y(i)= b(1)*x(i) + b(2)*x(i-1) + b(3)*x(i-2);
        else y(i)=0;
    end
    end

plot(t,y/Hw,'r');grid;title('Problem 1. Plot of filtered signal');
xlabel('Time (sec)'); ylabel('Signal');grid;
legend('Input','Output','Location','NorthWest');

%Find period
N=length(y);
for n=10:N
    if(y(n+1)-y(n)>=0)
        if(y(n+2)-y(n+1)<=0)
            n1=n;
            break
        end
    end
end

for n=n1+1:N
    if(y(n+1)-y(n)>=0)
        if(y(n+2)-y(n+1)<=0)
            n2=n;
            break
        end
    end
end
f=1/((n2-n1)*Ts); fprintf('Ouput Signal Frequency (Hz) = %6.2f\n',f);
plot(n1*Ts,y(n1+1)/Hw,'ko',n2*Ts,y(n2+1)/Hw,'ko');
pause(1.0);


%% Problem 2
clear all; figure(2); clf(2); hold on;
fprintf('\nProblem 2\n');
f1=500; f2=1250; f3=1800; fs=25000;     %Frequencies
T1=1/f1; T2=1/f2; T3=1/f3; Ts=1/fs; t=[0:Ts:4*max([T1,T2,T3])];
x=cos(2*pi*f1*t)+ cos(2*pi*f2*t)+cos(2*pi*f3*t);
plot(t,x,'b--');grid;title('Plot of sinusoids at 500, 1250, and 1800 Hz');
xlabel('Time (sec)'); ylabel('Signal');
pause(1.0);

r=[exp(1j*2*pi*f1/fs), exp(1j*2*pi*f3/fs)]; %Zero locations (conjugate pairs) of 500,1800
b=poly([r,conj(r)]);        %Denominator polynomial with real coefficients
rb=roots(b);                %Show roots
M=length(b);                %Number of filter coefficients
Hw=abs(sum(b.*exp(-1j*2*pi*f2/fs).^[0:M-1])); %Gain at f2
N=length(x); k=0;           %Implement FIR filter
for n=1:N
     y(n)=0;
        
    for k=1:M
        if(n-k > M-1) y(n)= y(n)+b(k)*x(n-k+1);
        else y(n)=0;
        end
    end
end

hold on;
plot(t,y/Hw,'r');grid;title('Problem 2. Plot of filtered signal');
xlabel('Time (sec)'); ylabel('Signal');grid;
legend('Input','Output','Location','NorthWest');

%Find period
N=length(y);
for n=50:N
    if(y(n+1)-y(n)>=0)
        if(y(n+2)-y(n+1)<=0)
            n1=n;
            break
        end
    end
end

for n=n1+1:N
    if(y(n+1)-y(n)>=0)
        if(y(n+2)-y(n+1)<=0)
            n2=n;
            break
        end
    end
end
f=1/((n2-n1)*Ts); fprintf('Output Signal Frequency (Hz) = %6.2f\n',f);
plot(n1*Ts,y(n1+1)/Hw,'ko',n2*Ts,y(n2+1)/Hw,'ko');
pause(1.0);

%% Problem 3
figure(3); clear all; clf(3); hold on;
fprintf('\nProblem 3\n');
f1=500; f2=1250; f3=1800; fs=25000;     %Frequencies
T1=1/f1; T2=1/f2; T3=1/f3; Ts=1/fs; t=[0:Ts:4*max([T1,T2,T3])];
x=cos(2*pi*f1*t)+ cos(2*pi*f2*t)+cos(2*pi*f3*t);
plot(t,x,'b--');grid;title('Plot of sinusoids at 500, 1250, and 1800 Hz');
xlabel('Time (sec)'); ylabel('Signal');
pause(1.0);
r=exp(1j*2*pi*f2/fs);       %Zero locations (conjugate pairs) of 1250 Hz
b=poly([r,conj(r)]);        %Denominator polynomial with real coefficients
rb=roots(b);                %Show roots
M=length(b);                %Number of filter coefficients
Hw=abs(sum(b.*exp(-1j*2*pi*f1/fs).^[0:M-1])); %Add gains |H(f1)| + |H(f3)|
Hw=(Hw+abs(sum(b.*exp(-1j*2*pi*f3/fs).^[0:M-1]))); %Approximate(?) adjustment
N=length(x); k=0; %Implement FIR filter
for n=1:N
    y(n)=0;
    for k=1:M
        if(n-k > M-1) y(n)= y(n)+b(k)*x(n-k+1);
        else y(n)=0;
        end
    end
end

hold on;
plot(t,y/Hw,'r');grid;title('Problem 3. Plot of filtered signal');
xlabel('Time (sec)'); ylabel('Signal');grid;
legend('Input','Output','Location','NorthWest');

%Find period
N=length(y);
for n=50:N
    if(y(n+1)-y(n)>=0)
        if(y(n+2)-y(n+1)<=0)
            n1=n;
            break
        end
    end
end

n2=0;n3=0;
for n=n1+1:N
    if(y(n+1)-y(n)>=0)
        if(y(n+2)-y(n+1)<=0)
            if(n2==0) n2=n; 
            end
            if(abs(y(n+1)-y(n1))<0.01)n3=n; 
                break; 
            end
        end
    end
end
fprintf(' Measurements only approximate since two tones present:\n');
fx=1/((n2-n1)*Ts); fprintf(' Output Signal Frequency (Hz) = %6.2f\n',fx);
fy=1/((n3-n1)*Ts); fprintf(' Output Signal Frequency (Hz) = %6.2f\n',fy);
plot(n1*Ts,y(n1+1)/Hw,'ko',n2*Ts,y(n2+1)/Hw,'ko',n3*Ts,y(n3+1)/Hw,'ko');
pause(1.0);

%% Plot DTFT of input x(n) and output y(n) 
figure(4); clf(4); hold on;
w=pi*[0:1:500]/500; f=(fs/2)*[0:500]/500;
X=dtft(x,[0:length(x)-1],w);
Y=dtft(y,[0:length(y)-1],w);

subplot(2,1,1); plot(f,abs(X)/max(abs(X)),'b'); %Normalize
axis([0 fs/10 0 1]); grid on;
title('Plot of |X(f)| showing 500, 1200 and 1800 Hz frequencies');
xlabel('frequency (Hz)'); ylabel('|X(f)|');

subplot(2,1,2); plot(f,abs(Y)/max(abs(Y)),'r'); %Normalize
axis([0 fs/10 0 1]); grid on;
title('Plot of |Y(f)| showing suppression of 1250Hz signal in output');
xlabel('frequency (Hz)'); ylabel('|Y(f)|');
hold off;