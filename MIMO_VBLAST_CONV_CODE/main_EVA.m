clear;
close all;
clc;

sampleRate = 1e6 * 2;                                % Sample rate (Hz)
sampleTime = 1/sampleRate;
delayVector = [0 2 3 4 7 11 17 25]*1e-7;                          % Discrete delays of four-path channel (s)
gainVector  = [0 -1.4 -3.6 -0.6 -9.1 -7 -12 -17];
max_doppler = 0;

nt = 2;

CHAN = rayleighchan( sampleTime , max_doppler, delayVector , gainVector );

NP = 1;    %%%  noise power  db

cplength = 16;
ber1 = zeros(1,5,1);
ber2 = zeros(1,5,1);
%ber3 = zeros(1,7,1);

Nseeds = 1000;
dB = 0;
tic;
ind3 = 1;
for M_order = [64]
    for seed = (1:0+Nseeds)
        ind2 = 1;
        
        for SP = (16:5:36)
            ind1 = 1;
            for nr = [2]
                ber1(ind1,ind2,ind3)= ofdm_mimo_zf( nt , nr , M_order , 256 , seed , CHAN , SP , NP ,dB, cplength) + ber1(ind1,ind2,ind3);
                
                ber2(ind1,ind2,ind3) =fbmc_mimo_zf( nt , nr , M_order , 256 , seed , CHAN , SP , NP ,dB)+ ber2(ind1,ind2,ind3);
                
                ind1 = ind1 + 1;
            end
            ind2 = ind2 + 1;
        end
        seed
    end
    M_order
    ind3 = ind3 + 1;
end
toc;
n= (15:5:35);
ber1 = ber1 ./ Nseeds;

ber2 = ber2 ./ Nseeds;

%ber3 = ber3 ./ Nseeds;

ber = [ber1;ber2];
%save('EVA_64_mimo_4_4_SNR_15_35.mat','ber');
%exit;

plot(n,log10(ber1.'),'--',n,log10(ber2.'),'-*');
%legend('ofdm-1','ofdm-2','ofdm-4','ofdm-8','fbmc94-1','fbmc94-2','fbmc94-4','fbmc94-8','fbmc92-1','fbmc92-2','fbmc92-4','fbmc92-8');
%xlabel('SNR(dB)')
%ylabel('log(BER)')