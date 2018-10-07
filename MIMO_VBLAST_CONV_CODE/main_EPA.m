clear;
close all;
clc;
sampleRate = 1e6 * 20;                                % Sample rate (Hz)
sampleTime = 1/sampleRate;
delayVector = [0 3 7 9 11 19 41]*1e-8;                          % Discrete delays of four-path channel (s)
gainVector  = [0 -1 -2 -3 -8 -17 -21];
max_doppler = 0;
%M_order = 2;
nt = 2;
%seed = 4;
CHAN = rayleighchan( sampleTime , max_doppler, delayVector , gainVector );
%CHAN = channel( sampleTime ,  delayVector , gainVector,seed);
%SP = 31;   %%% signal power  db
NP = 1;    %%%  noise power  db
%ind2 = 1;
%maxDelay=delayVector(end);
% cplength = ceil( 2 * maxDelay / sampleTime );
cplength = 12;
ber1 = zeros(1,7,1);
ber2 = zeros(1,7,1);
%ber3 = zeros(1,7,1);

Nseeds = 2000;
dB = 0;
tic;
ind3 = 1;
for M_order = [64]
    for seed = (1:0+Nseeds)
        ind2 = 1;
        
        for SP = (1 : 5 : 31)
            ind1 = 1;
            for nr = [2]
                ber1(ind1,ind2,ind3)= ofdm_mimo_zf( nt , nr , M_order , 64 , seed , CHAN , SP , NP ,dB, cplength) + ber1(ind1,ind2,ind3);
                
                ber2(ind1,ind2,ind3) =fbmc_mimo_zf( nt , nr , M_order , 64 , seed , CHAN , SP , NP ,dB)+ ber2(ind1,ind2,ind3);
                
                %ber3(ind1,ind2,ind3) =fbmc_mimo_zf( 3, nt , nr , M_order , 64 , seed , CHAN , SP , NP ,dB)+ ber3(ind1,ind2,ind3);
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
n= (0:5:30);
ber1 = ber1 ./ Nseeds;

ber2 = ber2 ./ Nseeds;

%ber3 = ber3 ./ Nseeds;

ber = [ber1;ber2];
% save('ber_16_64_256_div_2_4.mat','ber');
% exit;

plot(n,log10(ber1.'),'--',n,log10(ber2.'),'-*');
%legend('ofdm-1','ofdm-2','ofdm-4','ofdm-8','fbmc94-1','fbmc94-2','fbmc94-4','fbmc94-8','fbmc92-1','fbmc92-2','fbmc92-4','fbmc92-8');
%xlabel('SNR(dB)')
%ylabel('log(BER)')
