clear;
close all;
clc;

M = 16;
k = log2(M);

trellis = poly2trellis([5 4],[23 35 0; 0 5 13]);
traceBack = 16;
codeRate = 2/3;

convEncoder = comm.ConvolutionalEncoder('TrellisStructure',trellis);
vitDecoder = comm.ViterbiDecoder('TrellisStructure',trellis, ...
    'InputFormat','Hard','TracebackDepth',traceBack);

qamModulator = comm.RectangularQAMModulator('BitInput',true,'NormalizationMethod','Average power');
qamDemodulator = comm.RectangularQAMDemodulator('BitOutput',true,'NormalizationMethod','Average power');

errorRate = comm.ErrorRate('ReceiveDelay',2*traceBack);

ebnoVec = 0:2:10;
errorStats = zeros(length(ebnoVec),3);

for m = 1:length(ebnoVec)
    snr = ebnoVec(m) + 10*log10(k*codeRate);
    %while errorStats(m,2) <= 100 && errorStats(m,3) <= 1e7
        dataIn = randi([0 1],1000000,1);
        dataEnc = step(convEncoder,dataIn);
        txSig = step(qamModulator,dataEnc);
        rxSig = awgn(txSig,snr);
        demodSig = step(qamDemodulator,rxSig);
        dataOut = step(vitDecoder,demodSig);
        errorStats(m,:) = step(errorRate,dataIn,dataOut);
    %end
    reset(errorRate)
end

berUncoded = berawgn(ebnoVec','qam',M);

semilogy(ebnoVec,[errorStats(:,1) berUncoded])
grid
legend('Coded','Uncoded')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
