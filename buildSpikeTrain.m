% buildSpikeTrain.m
%
% Inputs    spikeTimes: a vector or a cell array of vectors containing 
%               spike times. 
%           zeroTime: the time of the event to which you wish to lock.
%           startTime: how much time before the zeroTime do you wish to
%               include?
%           endTime: how much time after the zeroTime do you wish to
%               include?
%           dt: bin size for the spike train
%
% Output    Sm: the spike train (vector for single trial or matrix for
%               multiple trials).
%
% DKW, 1.4.16

function Sm = buildSpikeTrain(spikeTimes,zeroTime,startTime,endTime,dt)

if iscell(spikeTimes)
    numtrials=length(spikeTimes);
elseif isvector(spikeTimes)
    spikeTimes = {spikeTimes};
    numtrials=1;
elseif ismatrix(spikeTimes)
    error('spikeTimes input needs to be either a vector or a cell array of vectors, not a matrix')
end

for tr=1:numtrials
    
    spTs = spikeTimes{tr};
    
    spTs=spTs-double(zeroTime);
    spTs=spTs(spTs>startTime & spTs<endTime);
    spTs=spTs+abs(startTime);
    
    T = endTime-startTime;
    if startTime<=0, T=T+1; end
    
    S = sparse(zeros(1,ceil(T/dt)));
    for j=1:length(spTs)
        S(ceil((spTs(j))/dt)) = S(ceil(spTs(j)/dt))+1;
    end
    
    Sm(tr,:)=S;
end

Sm=full(Sm);
