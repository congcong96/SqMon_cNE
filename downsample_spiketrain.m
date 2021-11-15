function dsspktrain = downsample_spiketrain(spktrain, dsfactor)

[M, N] = size(spktrain);
dsN = ceil(N/dsfactor);
dsspktrain = zeros(M, dsN);

%append zeros to end for division
rem = mod(N, dsfactor);
if rem ~= 0
    spktrain = [spktrain zeros(M, dsfactor - rem)];
end

for i = 1:M
    summat = reshape(spktrain(i,:),dsfactor, dsN);
    dsspktrain(i,:) = sum(summat);  
end