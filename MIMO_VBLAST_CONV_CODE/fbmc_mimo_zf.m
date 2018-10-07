function [avg_ber,rec_msg_raw,msg_raw] = fbmc_mimo_zf(Nt,Nr,M_order,Ncarriers,seed,CHAN,SP,NP,NP_offset)
%QAM order

%channel estimation
pilot=ones(Ncarriers,1);

[fb1] = fbmc_fb_k32_33(Ncarriers);

fb2=fb1;
sent1 = fbmc_trans(pilot,Ncarriers,fb1);
sca = 10^(SP/20);
sent1 = sca*sent1;
sent = [sent1; zeros(Nt-1,length(sent1))];                   %transmit from 1 antenna and zero from others

for i=1:Nt
    y(:,:,i) = mimo_chan(Nt , Nr , CHAN, sent.', seed);
    %y(:,:,i) = mimo_chan_c(Nt , Nr , CHAN, sent.', seed);
    sent = circshift(sent,[1,0]);                            %transmit from another antenna
end

n_sd_est = 10^((NP+NP_offset)/20);
y = y + n_sd_est*randn(size(y));


for n_rec=1:Nr
    for n_tr = 1:Nt
        [temp] = fbmc_rec(y(:,n_rec,n_tr).',Ncarriers,fb2);      %inverse fbmc on signal received from n_rec when transmitted from n_tr
        H(n_rec,n_tr,:) = temp/sca;
    end
end

%% Eigen Beamforming obtaining U S V

BC = 50*log2(M_order);

msg_raw = randi([0 1],Nt*BC*Ncarriers,1);
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
msg = reshape(msg,[Ncarriers,code_r*BC/log2(M_order),Nt]);

xx = randperm(Ncarriers);
yy = eye(Ncarriers);
yy = yy(:,xx);

for l=1:Nt
    
   msg(:,:,l) = yy*msg(:,:,l);
    
end


for n_tr = 1:Nt
     sent_p(n_tr,:) = fbmc_trans(msg(:,:,n_tr) , Ncarriers , fb1);
end
y_rec = mimo_chan(Nt , Nr , CHAN, sent_p.', seed);
%y_rec = mimo_chan_c(Nt , Nr , CHAN, sent_p.', seed);
%% noise addition
n_sd = 10^(NP/20);
s = size(y_rec);
y_rec = y_rec + n_sd/sqrt(2) * (randn(s(1),s(2)) + 1j*randn(s(1),s(2)));


%%%-- average shift
%%%


for n_rec = 1:Nr
    rec_msg(n_rec,:)=reshape(fbmc_rec(y_rec(:,n_rec),Ncarriers,fb2),[1,BC*Ncarriers*code_r/log2(M_order)])/sca;
end

%POSTCODING
rec_msg_p = zeros(code_r*BC*Ncarriers/log2(M_order),Nt);

nn = BC*code_r/log2(M_order);

for sc = 1:Ncarriers
        rec_msg_p(sc:Ncarriers:end,:)=vblast_arr(rec_msg(:,sc:Ncarriers:end),H(:,:,sc)).';
end

hDeMod = comm.RectangularQAMDemodulator('ModulationOrder',M_order,'BitOutput',true,'NormalizationMethod','Average Power');

rec_msg_p = reshape(rec_msg_p,[Ncarriers , BC*code_r/log2(M_order) , Nt]);

for l=1:Nt
    
   rec_msg_p(:,:,l) = yy'*rec_msg_p(:,:,l);
    
end

rec_msg_p = reshape(rec_msg_p,[Nt*BC*Ncarriers*code_r/log2(M_order),1]);
rec_msg_raw = step(hDeMod,rec_msg_p);
rec = step(vitDecoder,rec_msg_raw);

avg_ber = biterr(msg_raw(1:end-traceBack),rec(traceBack+1:end))/(length(msg_raw)-traceBack);


end