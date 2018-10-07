function [avg_ber,rec_msg_raw,msg_raw] = ofdm_mimo_ebf(Nr,M_order,Ncarriers,seed,CHAN,SP,NP,NP_offset,cp_length)
%QAM order
Nt = 1;
%N_max=64; 
div=min(Nt,Nr);

%channel estimation
pilot=ones(Ncarriers,1);

%sent1 = ofdm_trans(pilot,Ncarriers,cp_length);
ofdmMod = comm.OFDMModulator('FFTLength',Ncarriers,'CyclicPrefixLength',cp_length,'NumGuardBandCarriers',[0 ; 0]);
sca = 10^(SP/20);
sent = sca*sqrt(Ncarriers)*step(ofdmMod,pilot).';                   %transmit from 1 antenna and zero from others

y = mimo_chan(Nt , Nr , CHAN, sent.', seed);

n_sd_est = 10^((NP+NP_offset)/20);
y = y + n_sd_est*randn(size(y));

for n_rec = 1:Nr
    
    ofdmDemod = comm.OFDMDemodulator('FFTLength',Ncarriers,'CyclicPrefixLength',cp_length,'NumGuardBandCarriers',[0 ; 0]);
    temp = step(ofdmDemod,y(:,n_rec))/sqrt(Ncarriers)/sca;
    H(n_rec,:) = temp;
    % gives channel between n_tr tranmit and n_rec antenna
    
end

BC = 50*log2(M_order);
msg_raw = (randi([0 1],div*BC*Ncarriers,1));
%msg_raw = ones(Ncarriers,1);
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
%code_r = 1;

hMod = comm.RectangularQAMModulator('ModulationOrder',M_order,'BitInput',true,'NormalizationMethod','Average Power');
msg = sca*step(hMod,enc);

msg = reshape(msg,[Ncarriers,code_r*div*BC/log2(M_order),1]);

xx = randperm(Ncarriers);
yy = eye(Ncarriers);
yy = yy(:,xx);

msg = yy*msg;
%msg = sca * qammod(msg_raw,M_order,'UnitAveragePower',true);
%msg_p = zeros(Nt,BC*Ncarriers);

%sent_p(n_tr,:) = ofdm_trans(msg_p(n_tr,:) , Ncarriers , cp_length);
ofdmMod = comm.OFDMModulator('FFTLength',Ncarriers,'NumSymbols',code_r*BC/log2(M_order),'CyclicPrefixLength',cp_length,'NumGuardBandCarriers',[0 ; 0]);
sent_p=sqrt(Ncarriers)*step(ofdmMod,msg);

y_rec = mimo_chan(Nt , Nr , CHAN, sent_p, seed);
%y_rec = mimo_chan_c(Nt , Nr , CHAN, sent_p.', seed);

%y_rec = sent_p * ones(1,Nr);
%% noise addition
n_sd = 10^(NP/20);
s = size(y_rec);
y_rec = y_rec + n_sd/sqrt(2) * (randn(s(1),s(2)) + 1j*randn(s(1),s(2)));
ofdmDemod = comm.OFDMDemodulator('FFTLength',Ncarriers,'NumSymbols',code_r*BC/log2(M_order),'CyclicPrefixLength',cp_length,'NumGuardBandCarriers',[0 ; 0]);

rec_msg = zeros( Nr , code_r*div*BC*Ncarriers/log2(M_order) );

for n_rec = 1:Nr
    %rec_msg(n_rec,:)=reshape   (ofdm_rec(y_rec(:,n_rec).',Ncarriers,cp_length).',[1,BC*Ncarriers]);
    temp = step(ofdmDemod,y_rec(:,n_rec));
    %temp = yy'*temp;
    rec_msg(n_rec,:)= reshape(temp,[1,code_r*div*BC*Ncarriers/log2(M_order)])/sqrt(Ncarriers)/sca;    
end

%POSTCODING
rec_msg_p = zeros(1,code_r*BC*Ncarriers/log2(M_order));

for nc=1:Ncarriers
   rec_msg_p(nc:Ncarriers:end) = H(:,nc)'*rec_msg(:,nc:Ncarriers:end)/norm(H(:,nc)).^2;
    
end

%rec_msg_p=(diag(H'*rec_msg)./diag(H'*H));


% rx = reshape(rec_msg_p,[1 BC*Ncarriers*div]);
% scatter(real(rx),imag(rx));
% axis equal;
hDeMod = comm.RectangularQAMDemodulator('ModulationOrder',M_order,'BitOutput',true','NormalizationMethod','Average Power');
rec_msg_p = reshape(rec_msg_p,[Ncarriers, code_r*div*BC/log2(M_order)]);

rec_msg_p = yy'*rec_msg_p;

rec_msg_p = reshape(rec_msg_p,[code_r*div*BC*Ncarriers/log2(M_order),1]);
rec_msg_raw = step(hDeMod,rec_msg_p);
rec = step(vitDecoder,rec_msg_raw);
%errorRate = comm.ErrorRate('ReceiveDelay',2*traceBack);

%avg_ber = step(errorRate,msg_raw,rec);
%avg_ber = avg_ber(1);
avg_ber = biterr(msg_raw(1:end-traceBack),rec(traceBack+1:end))/(length(msg_raw)-traceBack);
end