clear all;
close all;
clc

%% Param?tres de simulation
M = 2;
nb = log2(M);               % nb de bits/symboles
Ts = 0.05;
Fe = 1/Ts;
N = 128;
Nu = 128;
Nb_OFDM_Symb = 500;         % nb symbole OFDM par trame
Nb_symb = Nb_OFDM_Symb*Nu;  % nb de symb au total
Nb_bits = Nb_symb*nb;       % nb de bits au total

snr = 0:15;
sigma = 10.^(-snr/10);
teb = zeros(1,length(snr));

CP = 16;
L = 16;

itr_max = 100000;
nb_err_max = 10000;

%% Calcul et tracé du TEB
i=1;
egal = 0;
cyc = 0;
for i=1:length(snr)
    itr = 0;
    nb_err = 0;
    
    while itr < itr_max && nb_err < nb_err_max
        itr = itr + 1;
        
        % Sequence binaire aléatoire
        bits = randi([0,1],Nb_bits/nb,nb);
        
        %Canal
        h = [1];
        bits_est = ofdm(bits, h, sigma(i), N, Nb_OFDM_Symb, M, CP, egal, cyc);

        nb_err = nb_err + sum(abs(bits-bits_est));
    end
    teb(i) = nb_err / itr / Nb_bits;
end

teb

figure
semilogy(snr, teb, 'r');
title('TEB');
hold on
semilogy(snr, 1/2*erfc(sqrt(10.^(snr/10))));
hold on
legend('TEB theorique', 'TEB experimental');

%% Test du préfixe circulaire

CP = 16;
L = 16;

i=1;
egal = 0;
cyc = 1;
nb_err = 0;

% Sequence binaire aléatoire
bits = randi([0,1],Nb_bits/nb,nb);

% Canal de Rayleigh
h = [1];

bits_est = ofdm(bits, h, 0, N, Nb_OFDM_Symb, M, CP, egal, cyc);

nb_err = sum(abs(bits-bits_est));

teb = nb_err / Nb_bits

%% Test de l'égaliseur pour L = 16 et sigma = 0

CP = 16;
L = 16;

i=1;
egal = 1;
cyc = 0;
nb_err = 0;

% Sequence binaire aléatoire
bits = randi([0,1],Nb_bits/nb,nb);

% Canal de Rayleigh
h = (randn(1,L) + 1i * randn(1,L)) / sqrt(2*L);

bits_est = ofdm(bits, h, 0, N, Nb_OFDM_Symb, M, CP, egal, cyc);

nb_err = sum(abs(bits-bits_est));

teb = nb_err / Nb_bits


%% Egaliseur de canal - CP = L = 16

CP = 16;
L = 16;

i=1;
egal = 1;
cyc = 1;
for i=1:length(snr)
    itr = 0;
    nb_err = 0;
    
    while itr < itr_max && nb_err < nb_err_max
        itr = itr + 1;
        
        % Sequence binaire aléatoire
        bits = randi([0,1],Nb_bits/nb,nb);
        
        % Canal de Rayleigh
        h = (randn(1,L) + 1i * randn(1,L)) / sqrt(2*L);
        
        bits_est = ofdm(bits, h, sigma(i), N, Nb_OFDM_Symb, M, CP, egal, cyc);

        nb_err = nb_err + sum(abs(bits-bits_est));
    end
    teb(i) = nb_err / itr / Nb_bits;
end

teb

fig = openfig('bertool_bpsk.fig');
title('TEB pour L = CP = 16');
hold on
semilogy(snr, teb, 'r');
hold on
semilogy(snr, 1/2*erfc(sqrt(10.^(snr/10))));
hold on
legend('TEB theorique', 'TEB experimental');

%% Egaliseur de canal - CP =8 et L = 16

CP = 8;
L = 16; 

i=1;
egal = 1;
cyc = 1;
for i=1:length(snr)
    itr = 0;
    nb_err = 0;
    
    while itr < itr_max && nb_err < nb_err_max
        itr = itr + 1;
        
        % Sequence binaire aléatoire
        bits = randi([0,1],Nb_bits/nb,nb);
        
        % Canal de Rayleigh
        h = (randn(1,L) + 1i * randn(1,L)) / sqrt(2*L);
        
        bits_est = ofdm(bits, h, sigma(i), N, Nb_OFDM_Symb, M, CP, egal, cyc);

        nb_err = nb_err + sum(abs(bits-bits_est));
    end
    teb(i) = nb_err / itr / Nb_bits;
end

teb

fig = openfig('bertool_bpsk.fig');
title('TEB pour L = 16 et CP = 8');
hold on
semilogy(snr, teb, 'r');
hold on
semilogy(snr, 1/2*erfc(sqrt(10.^(snr/10))));
hold on
legend('TEB theorique', 'TEB experimental');

%% Figures
%figure(1)
%plot([0:length(sk)-1]*Ts,real(sk))
%figure(2)
%hist(real(sk(:)),100);
%figure(3)
%hist(imag(sk(:)),100);
figure(4)
auto_cor_real = xcorr(real(sk(:)));
l = length(auto_cor_real);
plot((-l/2)+0.5:(l/2) ,auto_cor_real);
figure(5)
auto_cor_imag = xcorr(imag(sk(:)));
plot((-l/2)+0.5:(l/2), auto_cor_real);
figure(6)
int_cor_real_imag = xcorr(real(sk(:)), imag(sk(:)));
plot((-l/2)+0.5:(l/2), int_cor_real_imag);















