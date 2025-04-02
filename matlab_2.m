clear all
close all
%question 1:finding p using parameters:
delta_f=10^10;%10G[bps]
T=300;%[k]
R_L=50;%[ohm]
BER=10^-9;
R=1;%[w/A]
Fn=3;%[lin]
q = sqrt(2) * erfcinv(2 * BER);%BER=1/2erfc(q/sqrt(2))
K_B=1.38*10^-23;%constant boltzman
p_rec=(q/R)*sqrt(4*K_B*T*Fn*delta_f/R_L)%lin
fprintf('Prec_at_BER=10^-9 = %.4g\n',p_rec);
%% question 2:ploting p as function of q:
q=linspace(2,8,10^3);
BER=0.5*erfc(q./sqrt(2));
p_rec=q.*sqrt(4*K_B*T*Fn*delta_f/R_L)/R;%lin
figure(2)
plot(p_rec, BER, 'LineWidth', 1.5);
grid on;
xlabel('p_rec [W]');
ylabel('BER');
title('p_rec as a function of BER');
%% question 3
BER_3 = 10^-3;
[min_diff, target_idx] = min(abs(BER - BER_3));
P_rec3 = p_rec(target_idx);
P_in3 = P_rec3*2;
fprintf('Prec_at_BER=10^-3 = %.4g\n',P_in3);
symbols = randi([0, 1], 1, 100000);
I_3 = symbols.*R*2*P_rec3; %[A] current for each symbol without noise
I_d3 = P_in3/2; %[A] decision current
sig_3 = sqrt((4*K_B*T/R_L)*Fn*delta_f); % only thermal noice
I_3_noise = I_3 + normrnd(0,sig_3,1,100000); %[A] current with guassian noise for each symbol

errors = 0;
for i=1:100000
 if symbols(i) == 0
     if I_3_noise(i)>I_d3
         errors = errors + 1;
     end
 end
 if symbols(i) == 1
     if I_3_noise(i)<I_d3
         errors = errors +1;
     end
 end
end
fprintf('Number of errors question 3: %d\n',errors )

%plot
figure(2)
hold on
histogram(I_3_noise,100)
xlabel('I_{noise}[A]')
ylabel('Nuber of Events')
title('Histogram For P0, P1 : PIN Detector')
xline(I_d3,'--','I_d')
hold off
%% question 4 + 5
%parmeters
M = 10;
k_a = 0.7;
F_a = k_a*M + (1-k_a)*(2-1/M);
e = 1.6*10^-19; %[j]

%תוחלת
I_0 = 0;
I_1 = 2*M*R*P_rec3;
%שונות
sig_0 = sqrt((4*K_B*T/R_L)*Fn*delta_f);
sig_1 = sqrt((4*K_B*T/R_L)*Fn*delta_f + 4*e*M^2*R*P_rec3*F_a*delta_f);
%זרם החלטה
I_d4 = (sig_0*I_1 + sig_1*I_0)/(sig_0 + sig_1);
%הסתברות שגיאה
Q = (I_1 - I_0)/(sig_0 + sig_1);
BER_4 = 0.5*erfc(Q/sqrt(2));

symbols = randi([0, 1], 1, 100000);
I_4 = symbols.*I_1; %[A] current for each symbol without noise
I_4_noise = zeros(1,100000); 
errors2 = 0;
for i=1:100000
 if symbols(i) == 0
     I_4_noise(i) = I_4(i) + normrnd(0,sig_0);
     if I_4_noise(i)>I_d4
         errors2 = errors2 + 1;
     end
 end
 if symbols(i) == 1
     I_4_noise(i) = I_4(i) + normrnd(0,sig_0+sig_1);
     if I_4_noise(i)<I_d4
         errors2 = errors2 + 1;
     end
 end
end
fprintf('Number of errors question 5: %d\n',errors2 )

%plot
figure(3)
hold on
histogram(I_4_noise,100)
xlabel('I_{noise}[A]')
ylabel('Nuber of Events')
title('Histogram For P0, P1 : PAD Detector')
xline(I_d4,'--','I_d')
hold off
