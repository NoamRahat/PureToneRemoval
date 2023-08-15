
close all
clc
clear
%% Generation of noisy signal
% Enter your ID 
ID = 205918360;
[inputSignal,fs,SNR_in] = inputSignalBuilder(ID);
soundsc(inputSignal,fs)
figure();plot(0:length(inputSignal)-1,inputSignal);
xlabel('n','fontsize',16);
ylabel('signal','fontsize',16);

audiowrite(['Input_' num2str(ID) '.wav'],inputSignal,fs)
[x, fs]= audioread('about_time.wav');
SNR_in = 10*log10(mean(x.^2)/mean((inputSignal-x).^2));

%% Noise frequency detection
%   In this part the input signal is examined. We assume a pure tone
%   disturbance of cos(w_0*n), and have to locate w_0=(2*pi/N)*k0
%   the cosine wave is periodic with N=512
%% Fourier Discrete Transform - Detect w_0 from the last frame of the signal
Nframe=512;
x_last_frame=inputSignal((end-Nframe+1):end);
% plot x_last_frame DTFT to detect w_0
n=1:512;
[X,omega]=my_DTFT(x_last_frame,n,512);
plot(omega/pi,abs(X));
[~,w0_ind] = max(abs(X));
w_0 = abs(omega(w0_ind));
% dont forget to define w axis to plot the DTFT in right scale

%% Band stop with III implementations
%%%%% Implemetation I : perfect filtering (FIR)
N = 1000;
n = -N:N;
B = pi/65;

tic;
h_1 = (2*cos(w_0*n).*sin(B*n))./(pi*n);
h_1(N+1)=B/pi;

% Note- you can use conv() function to filter the signal. 
% use the option 'same' to get the same output length.
% for example if you have: input-x filter-h:
% y= conv(x,h,'same') 

y_1 = inputSignal-conv(inputSignal,h_1,'same');
toc;

audiowrite(['Output_I_' num2str(ID) '.wav'],y_1,fs)
SNR_out = 10*log10(mean(x.^2)/mean((y_1-x).^2));

%%DTFT h1
n=1:2001;
[H1,omega1]=my_DTFT(h_1,n,512);
plot(omega1/pi,abs(H1));
xlabel('frequency in units of pi');
ylabel('magnitude');
title('H_1 DTFT');

%%% Freaquency response- H_1(e^jw)

%TODO
%%
%%%%% Implemetation II : ZOH design (FIR)
N = 100;
n = -N:N;


h_2 = 2*cos(w_0*n)/(2*N+1);

y_2 = inputSignal-conv(inputSignal,h_2,'same'); 



audiowrite(['Output_II_' num2str(ID) '.wav'],y_2,fs)
SNR_out = 10*log10(mean(x.^2)/mean((y_2-x).^2))

%%DTFT h2
n=1:201;
[H2,omega2]=my_DTFT(h_2,n,512);
plot(omega2/pi,abs(H2));
xlabel('frequency in units of pi');
ylabel('magnitude');
title('H_2 DTFT');

%%% Freaquency response- H_2(e^jw)

%TODO

%%
%%%%% Implemetation III : recursive design (IIR)

alpha = 0.999;
z_1=0;z_2=0; % initial rest

for n=1:length(inputSignal)
    z_1 =alpha*exp(1j*w_0)*z_1+(1-alpha)*inputSignal(n) ;
    z_2 =alpha*exp(-1j*w_0)*z_2+(1-alpha)*inputSignal(n);

    y_3(n,1) =inputSignal(n)-z_1-z_2;
end


audiowrite(['Output_III_' num2str(ID) '.wav'],y_3,fs)
SNR_out = 10*log10(mean(x.^2)/mean((y_3-x).^2))

%%DTFT of h3:
n = -length(inputSignal)/2+1:1:length(inputSignal)/2;
h_3 = -2.*(1-alpha).*(alpha.^n).*heaviside(n).*cos(n*w_0);
h_3(14400) = 1;
[H_3,omega_3] = my_DTFT(h_3,n,5000);
figure; plot((omega_3/pi),H_3)
        xlabel('frequency in units of pi');
        ylabel('magnitude');
        title('H_3 DTFT');


%% Performace evaluation:

[Grade, SNR_out_ref]= GradeMyOutput(ID,y_1,1);
[Grade, SNR_out_ref]= GradeMyOutput(ID,y_2,2);
[Grade, SNR_out_ref]= GradeMyOutput(ID,y_3,3);
