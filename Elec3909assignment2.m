%%
% Constants
c = 3*10^8;
e_o = 8.854*10^-12;
u_o = 4*pi*10^-7;
%%
j = sqrt(-1);
ZL = 15+j*35;
Z0 = 50;
betaL = linspace(0, 2*pi, 100);
Zin = Z0*(ZL*cos(betaL)+j*Z0*sin(betaL))./(Z0*cos(betaL)+j*ZL*sin(betaL));

figure(1);
plot(betaL/pi, real(Zin),'r', 'Linewidth',3);
hold on;
plot(betaL/pi, imag(Zin),'b--', 'Linewidth',3);
plot(betaL/pi, 0*real(Zin), 'g');
hold off;
xlabel('\beta l /\pi')
ylabel('Re\{Z_{in}\}, Im\{Z_{in}\}')

a = abs(((50+j*100)/(150+j*100)))^2;

hyperc =@(x) cosh(x);
hypers =@(x) sinh(x);
%%
%Calculations for B. Antenna Radome Design


lambda = (3*10^8) / (sqrt(4.9)*10*10^9); % lambda = 2pi/B --> B=w/Vp where vp = c/root(ur*er)
zero = lambda/2; % sin(x) = 0 @ intervals of m*lambda/2
zerom1 = zero;
zerom2 = 2*zero;
zerom3 = 3*zero; % third m value gets a thickness of at least 2cm

figure(2)
fplot(hyperc, [-10 10], 'r');
hold on;
fplot(hypers, [-10 10], 'b--');
hold off;

beta =@(f) 2*pi*f*sqrt(4.9*e_o*u_o/2)*sqrt(sqrt(1)+1); % simplified beta equation
n_fiberglass = sqrt(u_o/(4.9*e_o)); % impedance of fiberglass
% sinh(x) = -jsin(jx) allows j term infront of sine
nfiber =@(f) (n_fiberglass*(n_fiberglass*j*sin(beta(f)*zerom3)+377*cos(beta(f)*zerom3))/(n_fiberglass*cos(beta(f)*zerom3)+377*j*sin(beta(f)*zerom3))); % using equation derived in question 1 to solve for gamma
gamma_reflected =@(f) (nfiber(f) - 377)/(nfiber(f) + 377); % reflection calculation
reflected =@(f) abs(gamma_reflected(f))^2; % reflected ratio of pref/pinc
figure(3);

% Pref = |gamma|^2 * Pinc
% where Pinc is the incident wave (|Vo|^2)/2Z
fplot(reflected, [8*10^9 12*10^9]);
xlabel('Frequency (Hz)')
ylabel('Percentage Reflected')
%%
% Question 3
ZL = 100+j*100;
%z_nl = 100/50 + 100/50*j;

%Calculate a value of multiple of lambda to find a value of length that
%makes the imaginary portion zero.
lambda_3 = (3*10^8)/((5*10^9)*sqrt(10.2));
lc_test =@(test) test*lambda_3;
beta_3 = 2*pi*(5*10^9)*sqrt(10.2*e_o*u_o);
z_lc_test =@(test) (100*(100*j*sin(beta_3*lc_test(test))+ZL*cos(beta_3*lc_test(test)))/(100*cos(beta_3*lc_test(test))+ZL*j*sin(beta_3*lc_test(test))));
real_z =@(test) real(z_lc_test(test));
imag_z =@(test) imag(z_lc_test(test));
%the plot shows first value occurs before 0.1, so used 0.08 and found 0.088
%was the value
imag_0 = fzero(imag_z, 0.1);

%Using that value to calculate Z_lc
lc = imag_0*lambda_3;
z_lc = (100*(100*j*sin(beta_3*lc)+ZL*cos(beta_3*lc))/(100*cos(beta_3*lc)+ZL*j*sin(beta_3*lc)));

%Smith chart to find load => z = 4.5
zt = sqrt(50*real(z_lc));

% using RT/duroid 6010 @ er = 10.2 and h = 1.905mm

%
beta_var =@(freq) 2*pi*freq*sqrt(10.2*e_o*u_o);
z_lc_1 =@(freq) (100*(100*j*sin(beta_var(freq)*lc)+ZL*cos(beta_var(freq)*lc))/(100*cos(beta_var(freq)*lc)+ZL*j*sin(beta_var(freq)*lc)));

%
%Calculated values of width from online calculator, and lengths calculated
%using epsilon effective
length_QW = c/(4*(5*10^9)*sqrt(6.114));
width_QW = 0.137*10^-3;
length_lc = imag_0*c/((5*10^9)*sqrt(6.209));
width_lc = 0.2395*10^-3;
%

%width = 0.137218129513*10^-3;
%length = 6.14665051232*10^-3;

%
f_l = 3*10^9;
f_u = 7*10^9;
steps = linspace(3*10^9, 7*10^9, 400);
f = f_l;
step_size = 0.01*10^9;
i=1;
while (f < f_u)
    beta_2(i) = 2*pi*f*sqrt(10.2*e_o*u_o); % simplified beta equation
    zot(i) = (100*(100*j*sin(beta_2(i)*lc)+ZL*cos(beta_2(i)*lc))/(100*cos(beta_2(i)*lc)+ZL*j*sin(beta_2(i)*lc)));
    ztrans(i) = (zt*(zt*j*sin(beta_2(i)*lambda_3/4)+zot(i)*cos(beta_2(i)*lambda_3/4))/(zt*cos(beta_2(i)*lambda_3/4)+zot(i)*j*sin(beta_2(i)*lambda_3/4))); % using equation derived in question 1 to solve for gamma
    reflected_trans(i) = (ztrans(i) - 50)/(50 + ztrans(i));
    transmission(i) = 2*ztrans(i)/(50 + ztrans(i));
    f_data(i) = f;
    log_data(i) = 20*log(reflected_trans(i));
    i = i+1;
    f = f + step_size;
end


figure(4);
plot(f_data, ztrans);


figure(5);
plot(f_data, real(reflected_trans));
hold on;
plot(f_data, imag(reflected_trans));
hold off;
xlabel('Frequency (Hz)')
ylabel('Real and Imaginary Reflection Coefficient')

figure(6);
plot(f_data, 20*log(transmission));
xlabel('Frequency (Hz)')
ylabel('Transmission Coefficient (dB)')

figure(7);
fplot(real_z, [0 2]);
hold on;
fplot(imag_z, [0 2]);
hold off;
xlabel('Frequency (Hz)')
ylabel('Real and Imaginary Impedance (Ohms)')

figure(8);
plot(f_data, real(reflected_trans));
xlabel('Frequency (Hz)')
ylabel('Reflection Coefficient (dB)')

figure(9);
plot(f_data, real(ztrans));
hold on;
plot(f_data, imag(ztrans));
hold off;
