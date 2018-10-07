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

BC = 50 * log2(M_order);
msg_raw = logical(randi([0 1],BC*Ncarriers,1));
size(msg_raw);
sca = 10^(SP/20);

%trellis = poly2trellis([5 4],[23 35 0; 0 5 13]);
%traceBack = 16;
trellis = poly2trellis(7,[171 133]);
traceBack = 32;

convEncoder = comm.ConvolutionalEncoder('TrellisStructure',trellis);
vitDecoder = comm.ViterbiDecoder('TrellisStructure',trellis, ...
    'InputFormat','Hard','TracebackDepth',traceBack);

enc = step(convEncoder,msg_raw);

code_r = 2;

hMod = comm.RectangularQAMModulator('ModulationOrder',M_order,'BitInput',true,'NormalizationMethod','Average Power');
msg = sca * step(hMod,enc);
%msg = ones(1,512);

msg = reshape(msg,[Ncarriers,code_r*BC/log2(M_order)]);

xx = randperm(Ncarriers);
yy = eye(Ncarriers);
yy = yy(:,xx);

msg = yy*msg;

%msg = sca * qammod(msg_raw,M_order,'UnitAveragePower',true);
%msg_p = zeros(Nt,BC*Ncarriers);

sent_p = fbmc_trans(msg , Ncarriers , fb1);

y_rec = mimo_chan(Nt , Nr , CHAN, sent_p.', seed);

%% noise addition
n_sd = 10^(NP/20);
s = size(y_rec);
y_rec = y_rec + n_sd/sqrt(2) * (randn(s(1),s(2)) + 1j*randn(s(1),s(2)));

rec_msg = zeros(1,code_r*BC*Ncarriers/log2(M_order));

for n_rec = 1:Nr
    rec_msg(n_rec,:)=reshape(fbmc_rec(y_rec(:,n_rec).',Ncarriers,fb2).',[1,code_r*BC*Ncarriers/log2(M_order)])/sca;
end

rec_msg_p = zeros(1,code_r*BC*Ncarriers/log2(M_order));

for nc=1:Ncarriers
   rec_msg_p(nc:Ncarriers:end) = H(:,nc)'*rec_msg(:,nc:Ncarriers:end)/norm(H(:,nc)).^2;
    
end

rec_msg_p = reshape(rec_msg_p,[Ncarriers , code_r * BC/log2(M_order)]);

rec_msg_p = yy'*rec_msg_p;

hDeMod = comm.RectangularQAMDemodulator('ModulationOrder',M_order,'BitOutput',true,'NormalizationMethod','Average Power');
rec_msg_p = reshape(rec_msg_p,[code_r*div*BC*Ncarriers/log2(M_order),1]);
rec_msg_raw = step(hDeMod,rec_msg_p);
rec = step(vitDecoder,rec_msg_raw);
avg_ber = biterr(msg_raw(1:end-traceBack),rec(traceBack+1:end))/(length(msg_raw)-traceBack);

end