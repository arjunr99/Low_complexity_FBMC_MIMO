clear;
close all;
clc;
sampleRate = 1e6*2;                                % Sample rate (Hz)
sampleTime = 1/sampleRate;
delayVector = [0 1 2 5 16 23 50]*1e-7;                          % Discrete delays of four-path channel (s)
gainVector  = [-1 -1 0 0 -3 -5 -7];
max_doppler = 0;
%M_order = 2;

%seed = 4;
CHAN = rayleighchan( sampleTime , max_doppler, delayVector , gainVector );
%CHAN = channel( sampleTime ,  delayVector , gainVector,seed);
%SP = 31;   %%% signal power  db
NP = 1;    %%%  noise power  db

% maxDelay=delayVector(end);
% cplength = ceil( 2 * maxDelay / sampleTime );
cplength = 16;
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
        
        for SP = (0:5:31)
            ind1 = 1;
            for nr = [1]
                ber1(ind1,ind2,ind3) = ofdm_mimo_ebf( nr , M_order , 128 , seed , CHAN , SP , NP ,dB, cplength) + ber1(ind1,ind2,ind3);
                
                ber2(ind1,ind2,ind3) = fbmc_mimo_ebf(nr , M_order , 128 , seed , CHAN , SP , NP, dB)+ ber2(ind1,ind2,ind3);
                
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
n= (0:5:31);
ber1 = ber1 ./ Nseeds;

ber2 = ber2 ./ Nseeds;

%ber3 = ber3 ./ Nseeds;

ber = [ber1;ber2];
%save('ETU_64_siso_128_carr_sync.mat','ber');
%exit;

plot(n,log10(ber1.'),'--',n,log10(ber2.'),'-*');
% legend('ofdm-1','ofdm-2','ofdm-4','ofdm-8','fbmc97-1','fbmc97-2','fbmc97-4','fbmc97-8','fbmc97new-1','fbmc97new-2','fbmc97new-4','fbmc97new-8');
% xlabel('SNR(dB)')
% ylabel('log(BER)')
