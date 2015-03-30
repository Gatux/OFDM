function [ output ] = ofdm( input, h, sigma, N, Nb_OFDM_Symb, M, CP, egal, cyc)

% Emetteur
dcm = bi2de(input,'left-msb');

% M-PSK
Sk = zeros(N, Nb_OFDM_Symb);
Sk = pskmod(dcm,M);
Sk = reshape(Sk, [N, Nb_OFDM_Symb]);
sk = sqrt(N)*ifft(Sk,N);

if cyc == 0
    sk_cyc = sk;
else
    % Ajout du préfixe cyclique
    sk_cyc = zeros(CP+N, Nb_OFDM_Symb);
    sk_cyc(1:CP, :) = sk(end-CP+1:end, :);
    sk_cyc(CP+1:end, :) = sk;
end

% Canal 
noise = sqrt(sigma/2)*(randn(size(sk_cyc))+1i*randn(size(sk_cyc)));

if cyc == 0
    rk = filter(h,1,sk_cyc) + noise;
else
    Ncp = N + CP;
    H1 = zeros(Ncp,Ncp);
    H2 = zeros(Ncp,Ncp);
    for i = 1:Ncp
        for j = 1:Ncp
            if (i >= j && i-j+1 <= length(h))
                H1(i,j) = h(i-j+1);
            end
            if (i-j+Ncp+1 <= length(h))
                H2(i,j) = h(i-j+Ncp+1);
            end
        end
    end
    
    m = H1 + H2;
    rk = m*sk_cyc + noise;
end

% Recepteur
if cyc ~= 0
    % Suppression du préfixe
    rk(1:CP,:) = [];
end

Rk = fft(rk)/sqrt(N);

if egal == 0
    Rk_zf = Rk;
else
    % Zero-forcing
    H = fft(h,N);
    Rk_zf = zeros(size(Rk));
    for j = 1:N
        Rk_zf(j,:) = Rk(j,:)/H(j);
    end
end
    
Sk_est = pskdemod(Rk_zf,M);
output = de2bi(Sk_est,'left-msb');

end

