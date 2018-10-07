function [avg_ber,rec_msg_p] = fbmc_mimo_ebf(Nr,M_order,Ncarriers,seed,CHAN,SP,NP,NP_offset)
%QAM order
Nt = 1;
%N_max=64; 
div=min(Nt,Nr);

%channel estimation
pilot=ones(Ncarriers,1);

[fb1] = fbmc_fb_k32_33(Ncarriers);
fb2 = fb1;
sent1 = fbmc_trans(pilot,Ncarriers,fb1);
sca = 10^(SP/20);
sent1 = sca*sent1;
sent = [sent1; zeros(Nt-1,length(sent1))];                   %transmit from 1 antenna and zero from others

y = mimo_chan(Nt , Nr , CHAN, sent.', seed);                          %transmit from another antenna

n_sd_est = 10^((NP+NP_offset)/20);
y = y + n_sd_est*randn(size(y));

    for n_rec = 1:Nr
        [temp] = fbmc_rec(y(:,n_rec).',Ncarriers,fb2);
        %shift(n_rec,:)=sh;                  %inverse fbmc on signal received from n_rec when transmitted from n_tr
        H(n_rec,:) = temp/sca;                              % gives channel between n_tr tranmit and n_rec antenna   
    end


BC = 50*log2(M_order);
msg_raw = logical(randi([0 1],div*BC*Ncarriers,1));
size(msg_raw);
sca = 10^(SP/20);

hEnc = comm.BCHEncoder(15,5);
enc = step(hEnc,msg_raw);

hMod = comm.RectangularQAMModulator('ModulationOrder',M_order,'BitInput',true,'NormalizationMethod','Average Power');
msg = sca*step(hMod,enc);
msg = reshape(msg,[div,3*BC*Ncarriers/log2(M_order)]);
%msg = sca * qammod(msg_raw,M_order,'UnitAveragePower',true);
%msg_p = zeros(Nt,BC*Ncarriers);

sent_p = fbmc_trans(msg , Ncarriers , fb1);

y_rec = mimo_chan(Nt , Nr , CHAN, sent_p.', seed);

%% noise addition
n_sd = 10^(NP/20);
s = size(y_rec);
y_rec = y_rec + n_sd/sqrt(2) * (randn(s(1),s(2)) + 1j*randn(s(1),s(2)));

for n_rec = 1:Nr
    rec_msg(n_rec,:)=reshape(fbmc_rec(y_rec(:,n_rec).',Ncarriers,fb2).',[1,3*BC*Ncarriers/log2(M_order)])/sca;
end

%H1 = repmat(H,1,3*BC);
rec_msg_p = zeros(1,3*BC*Ncarriers/log2(M_order));
for nc=1:Ncarriers
   rec_msg_p(nc:Ncarriers:end) = H(:,nc)'*rec_msg(:,nc:Ncarriers:end);
    
end

%rec_msg_p=(diag(H1'*rec_msg)./diag(H1'*H1));

% rx = reshape(rec_msg_p,[1 BC*Ncarriers*div]);
% scatter(real(rx),imag(rx));
% axis equal;
%rec_msg_raw = qamdemod(rec_msg_p,M_order,'UnitAveragePower',true);
hDeMod = comm.RectangularQAMDemodulator('ModulationOrder',M_order,'BitOutput',true,'NormalizationMethod','Average Power');
rec_msg_p = reshape(rec_msg_p,[3*div*BC*Ncarriers/log2(M_order),1]);
rec_msg_raw = step(hDeMod,rec_msg_p);
hDec = comm.BCHDecoder(15,5);
rec = step(hDec,rec_msg_raw);
%rec_msg_raw = reshape(rec_msg_raw,[div,BC*Ncarriers]);

avg_ber = nnz(msg_raw-rec)/numel(msg);