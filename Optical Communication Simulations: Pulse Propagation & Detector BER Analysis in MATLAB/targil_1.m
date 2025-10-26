%defining parameters
b2=-20;%[ps^2/km]
L=50;%[km]
t=linspace(-500,500,1001);
w=linspace(-pi,pi,1001);
H_w_fiber=exp(-1j*b2*L*w.^2/2);

%% Q1

%defin pulse
T_0=25;%[ps]
A_0=1;
A_T = A_0*exp(-(t.^2/(2*T_0^2)));
A_W = fftshift(fft(A_T));

%finding A_out in the frequency domain-convelution in time=multiple in
%freq...
A_Wout = H_w_fiber.*A_W;
A_Tout = ifft(ifftshift(A_Wout));
P_OUT = (A_Tout).^2;
P= (A_T).^2;

%finding the pulse width beyond exp(-1) limit:
t_in = t(abs(P) >= exp(-1));
calc_T_0_in = round(t_in(end));
fprintf('The T_0 for the input of A0𝑒𝑥𝑝(−𝑡2/2𝑇0 2) is: %d [ps]\n',calc_T_0_in )

t_out = t(abs(P_OUT) >= max(abs(P_OUT))*exp(-1));
calc_T_0_out = round(t_out(end));
fprintf('The T_0 for the Output of A0𝑒𝑥𝑝(−𝑡2/2𝑇0 2) is: %d [ps]\n',calc_T_0_out )

%plots
figure(1)
hold on
plot(t,abs(P),'-','color','black');
plot(t,abs(P_OUT),'-','color','g');
xlabel('Time [ps]');
ylabel('Magnitude');
legend('Input signal','Output signal')
title('question 1:pulse in time domain')
figure(2)
hold on
plot(w,abs(A_W))
plot(w,abs(A_Wout))
xlabel('w')
ylabel('Magnitude')
legend('Input signal','Output signal')
title('question 1:pulse in frequency domain')
axis([-pi pi 0 70])
%% Q2

%define pulse
T_0=25;%[ps]
A_0=1;
A_T2 = A_0*sech(t/T_0);
A_W2 = fftshift(fft(A_T2));

%A_out
A_W2out = H_w_fiber.*A_W2;
A_T2out = ifft(ifftshift(A_W2out));
P_OUT = (A_T2out).^2;
P= (A_T2).^2;

%T_pulse
t_in = t(abs(P) >= interp1(t, (P), 25));
calc_T_0_in = round(t_in(end));
fprintf('The T_0 for the input of 𝐴0sech(t/T0) is: %d [ps]\n',calc_T_0_in )

t_out = t(abs(P_OUT) >= max(abs(P_OUT))*interp1(t, (P), 25));
calc_T_0_out = round(t_out(end));
fprintf('The T_0 for the Output of 𝐴0sech(t/T0) is: %d [ps]\n',calc_T_0_out )

%plots
figure(3)
hold on
plot(t,abs(P),'-','color','black');
plot(t,abs(P_OUT),'-','color','g');
xlabel('Time [ps]');
ylabel('Magnitude');
legend('Input signal','Output signal')
title('question 2:pulse in time domain,sech function')
figure(4)
hold on
plot(w,abs(A_W2))
plot(w,abs(A_W2out))
xlabel('w')
ylabel('Magnitude')
legend('Input signal','Output signal')
title('question 2:pulse in frequency domain,sech function')
axis([-pi pi 0 90])

%% Q3
% Parameters
alpha = 0.2; %[dB/km]
L_loss = 20; %[dB]
T_0 = 25; %[ps]
beta_2 = -20; %[ps^2/km]
L = L_loss / alpha; %[km]
fprintf('The fiber length L is: %.2f km\n', L);
z = L;
T_z = T_0 * sqrt(1 + (z * beta_2 / T_0^2)^2)%[ps];
fprintf('The broadened pulse width T(z) is: %.2f ps\n', T_z);
B = 1/(T_z*sqrt(8))*10^3; %[Gbit/s](convert 1/pico(-12) to G(9)
fprintf('The maximum bit rate B is: %.2f Gbit/s\n', B);
%% Q4

% parameters
t = -30*(10^6):1:30*(10^6);%[sec] 
L = 5;%[km] 
b2=-20;%[ps^2/km]
b1 = 4.83*(10^6);%[ps/km]
A_0 = 1; 
a = 30*(10^6);
w = 2*pi*(-a/(2*a+1):1/(2*a+1):a/(2*a+1));

%T_0 = 5
T_0 = 5;
A_T = A_0*exp(-(t.^2/(2*T_0^2))); 
A_w = fftshift(fft(A_T));
% transfer with beta1:
H_w_fiber = exp(-1j*b2/2*L*w.^2).*exp(-1j*b1*w*L);
A_wout = H_w_fiber .* A_w; 
A_Tout = ifft(ifftshift(A_wout)); 
P_out_5 = abs(A_Tout).^2; 
P_5 = abs(A_T).^2;
%finding width between e^-1 Tin:
threshold_in = exp(-1);
t_in = t(P_5 >= threshold_in);
T0_input = round(t_in(end));
fprintf('The pulse width T_0 for the input (T_0 = 5) is: %d [ps]\n', T0_input);
%finding width between e^-1 Tout:
threshold_out = max(abs(P_out_5)) * exp(-1);
t_out = t(abs(P_out_5) >= threshold_out);
mid_index = round(length(t_out) / 2);
T0_output = round(abs(t_out(end) - t_out(mid_index)));
fprintf('The pulse width T_0 for the output (T_0 = 5) is: %d [ps]\n', T0_output);

%T_0 = 25
T_0 = 25;
A_T = A_0*exp(-(t.^2/(2*T_0^2))); 
A_w = fftshift(fft(A_T)); 
% transfer with beta1:
H_w_fiber = exp(-1j*b2/2*L*w.^2).*exp(-1j*b1*w*L);
A_wout = H_w_fiber .* A_w;
A_Tout = ifft(ifftshift(A_wout)); 
P_out_25 = abs(A_Tout).^2; 
P_25 = abs(A_T).^2;
%finding width between e^-1 Tin:
threshold_in = exp(-1);
t_in = t(P_25 >= threshold_in);
T0_input = round(t_in(end));
fprintf('The pulse width T_0 for the input (T_0 = 25) is: %d [ps]\n', T0_input);
%finding width between e^-1 Tout:
threshold_out = max(abs(P_out_25)) * exp(-1);
t_out = t(abs(P_out_25) >= threshold_out);
mid_index = round(length(t_out) / 2);
T0_output = round(abs(t_out(end) - t_out(mid_index)));
fprintf('The pulse width T_0 for the output (T_0 = 25) is: %d [ps]\n', T0_output);

%plots
figure(5)
hold on
plot(t, P_5, 'b', 'LineWidth', 1.5)
plot(t, abs(P_out_5), 'r', 'LineWidth', 1.5)
xlabel('Time [ps]')
ylabel('Power')
title('T_0 = 5')
legend('Input Power', 'Output Power')


figure(6)
hold on
plot(t, P_25, 'b', 'LineWidth', 1.5)
plot(t, abs(P_out_25), 'r', 'LineWidth', 1.5)
xlabel('Time [ps]')
ylabel('Power')
title('T_0 = 25')
legend('Input Power', 'Output Power')
grid on

%% Q5

%defining parameters
a = 2000; 
T_0 = 25; 
L = 50;
t = -a:1:a;
w = 2*pi*(-a/(2*a+1):1/(2*a+1):a/(2*a+1));
B = 1/100;%[Gbit/ps](10G(10^10)/ps(10^12))
my_word=[1 0 1 0 0 0 1 1 0 1 0 0 1 0 1 0 0 0 1 0 0 0 0 1 0 1 0 0 1 0 1 0 ];

%define pulse
A_T = zeros(1,2*a+1);
A_T2 = zeros(1,2*a+1);
for i=1:32
A_T = A_T + my_word(i)*A_0*exp(-((t-(i-16)/B).^2)/(2*T_0^2));
A_T2 = A_T2 + my_word(i)*sech((t-(i-16)/B)/T_0);
end
A_W = fftshift(fft(A_T));
A_W2 = fftshift(fft(A_T2));
H_w_fiber=exp(-1j*b2/2*L*w.^2);

%A_out
A_Wout = H_w_fiber.*A_W; 
A_Tout = ifft(ifftshift(A_Wout));
A_Wout2 = H_w_fiber.*A_W2; 
A_Tout2 = ifft(ifftshift(A_Wout2)); 
P_OUT = (A_Tout).^2;
P= (A_T).^2;
P_OUT2 = (A_Tout2).^2;
P2= (A_T2).^2;


%plots
figure(7)
hold on
plot(t,abs(P),'b')
plot(t,abs(P_OUT),'r')
xlabel('Time [ps]')
ylabel('Magnitude')
legend('input','output')
title('gaussian')

figure(8)
hold on
plot(t,abs(P2),'b')
plot(t,abs(P_OUT2),'r')
xlabel('Time [ps]')
ylabel('Magnitude')
legend('input','output')
title('sech')
