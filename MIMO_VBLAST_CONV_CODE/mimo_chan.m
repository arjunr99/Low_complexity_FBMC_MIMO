function [rxSig,seed] = mimo_chan(numTrans , numRec , chan, input,seed)

%CHAN(1:numRec , 1:numTrans) = chan;
rxSig = zeros(size(input,1),numRec);
rng(seed);

for l=1:numRec
    for m=1:numTrans
        rxSig(:,l)= rxSig(:,l) + filter(chan , input(:,m));
    end
end

seed = rng();
end