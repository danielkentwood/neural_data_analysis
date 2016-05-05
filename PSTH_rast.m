
function handles = PSTH_rast(spikeTimes,time_params,varargin)
% TIME_PARAMS is a struct with the following fields.
%     time_params(1).zero_time=0;
%     time_params(1).start_time=-450;
%     time_params(1).end_time=250;
%     time_params(1).dt=5;
%     
%     NOTE: there should be as many instances of TIME PARAMS as there are
%     spike trains you want to plot.


% make sure spikeTimes is a cell array of cell arrays
if ~iscell(spikeTimes{1})
    numConds=1;
    spikeTimes = {spikeTimes};
else
    numConds=length(spikeTimes);
end

%% varargin default values (varargin is a struct with the following possible fields)
sep         = .075; % proportion of vertical separation between raster and psth
rasterSpace = .4; % vertical proportion of total plot that the raster takes up
axesDims    = [.15 .15 .7 .7]; % [left bottom width height]
MColor      = [1 0 0; 0 0 1; 0 1 .2; 1 .65 0]; % default colormap
errBars     = 1; % to turn psth error bars on and off
useSEs      = 0; % use SEM instead of 95% CI for the psth
smoothflag  = 1; % use smoothing
smoothtype  = 'gauss'; % can be either 'gauss' or 'spline'
gauss_sigma = 15;
splineOrder = 35; % for smoothing the psth with a spline
names       = {};


Pfields = {'sep', 'rasterSpace', 'axesDims', 'MColor', 'useSEs','splineOrder', 'errBars','smoothflag','smoothtype',...
    'gauss_sigma','names'};
for i = 1:length(Pfields) % if a params structure was provided as an input, change the requested fields
    if ~isempty(varargin)&&isfield(varargin{1}, Pfields{i}), eval(sprintf('%s = varargin{1}.(Pfields{%d});', Pfields{i}, i)); end
end
if ~isempty(varargin)  % if there is a params input
    fnames = fieldnames(varargin{1}); % cycle through field names and make sure all are recognized
    for i = 1:length(fnames)
        recognized = max(strcmp(fnames{i},Pfields));
        if recognized == 0, fprintf('fieldname %s not recognized\n',fnames{i}); end
    end
end

%% plot figure
psthSpace = 1-sep-rasterSpace;
hf = figure; hold on;
sp1 = subplot('position',[axesDims(1) (axesDims(2)+psthSpace*axesDims(4)+...
    sep*axesDims(4)) axesDims(3) rasterSpace*axesDims(4)]);
sp2 = subplot('position',[axesDims(1:2) axesDims(3) psthSpace*axesDims(4)]);
set(hf, 'color', [1 1 1]);

% start looping through conditions
j=0;
for c = 1:numConds
    % get current condition
    zt = time_params(c).zero_time;
    st = time_params(c).start_time;
    et = time_params(c).end_time;
    dt = time_params(c).dt;
    
    % get time vector
    tv = (st:dt:et)';
    
    % convert spike timestamps to spike train histogram
    sTms = spikeTimes{c};
    spikeTrains = buildSpikeTrain(sTms,zt,st,et,dt);
    sm_raw = spikeTrains;
    
    % compute the mean spike
    sm_mu = mean(sm_raw).*1000/dt;
    
    if smoothflag
        switch smoothtype
            case 'spline'
                sm_mu = lsSplineSmooth(tv,sm_mu,splineOrder)';
            case 'gauss'
                % create a function that convolves a spike train with a
                % gaussian kernel
                sm_gauss = gauss_spTrConvolve( spikeTrains, dt, gauss_sigma );
                sm_mu = mean(sm_gauss).*1000/dt;
        end
    end
    
    % compute the error bars
    if errBars
        if useSEs % use standard error
            sm_err = repmat(std(sm_raw)./sqrt(size(sm_raw,1)),2,1)' .*1000/dt;
            if smoothflag
                switch smoothtype
                    case 'spline'
                        sm_err = repmat(lsSplineSmooth(tv,std(sm_raw)./sqrt(size(sm_raw,1)),splineOrder),2,1)' .*1000/dt;
                    case 'gauss'
                        sm_err = repmat(std(sm_gauss)./sqrt(size(sm_gauss,1)),2,1)' .*1000/dt;
                end
            end
        else % use 95% CI 
            sm_err = bootci(1000,@mean,sm_raw)'.*1000/dt;
            if smoothflag
                switch smoothtype
                    case 'spline'
                        sm_err = lsSplineSmooth(tv,bootci(1000,@mean,sm_raw)',splineOrder)' .*1000/dt;
                    case 'gauss'
                        sm_err = bootci(1000,@mean,sm_gauss)'.*1000/dt;
                end
            end
            sm_err = sm_err-repmat(sm_mu',1,2);
            sm_err(:,2) = sm_err(:,2).*-1;
        end
    end
    
    % add the raster
    subplot(sp1)
    for t = 1:length(sTms)
        spTimes = sTms{t}';
        spTimes = spTimes(spTimes>st & spTimes<et);
        nsp = length(spTimes);
        ypl = repmat([j; j+.9], 1,  nsp);
        j=j+1;
        plot([spTimes'; spTimes'], ypl, 'color',MColor(c,:), 'linewidth', 1); hold on;
    end
    xlim([st et])
    % keep tabs of how many total rows there are in the raster
    rastRows(c) = size(sm_raw,1);
    
    % add the psth
    subplot(sp2)
    if errBars
        hold all
        [h(c).l,h(c).p] = boundedline(tv,sm_mu,...
            sm_err*-1,'cmap',MColor(c,:),'alpha');
    else
        hold all
        h(c).l = plot(tv,sm_mu,'color',MColor(c,:));
    end
    set(h(c).l,'linewidth',2);
end

% finalize the psth plot
subplot(sp2); axis tight
cxlm = get(sp2,'xlim');
cylm = get(sp2,'ylim');
axis manual
set(sp2,'ylim',[cylm(1)-(diff(cylm)*.025) cylm(2)+(diff(cylm)*.025)]);
set(sp2,'tickdir','out','YTick',[0 10.*max(unique(round(get(sp2,'YTick')./10)))],...
    'xlim',cxlm,'XTick',100.*unique(ceil(get(sp2,'XTick')./100)));
ylabel('Spikes/S')
xlabel('Time (ms)')
if ~isempty(names)
    lhPSTH = legend([h.l],names,'location','best');
    set(lhPSTH,'box','off')
end

% finalize the raster plot
totalRastRows=sum(rastRows);
subplot(sp1);
axis manual
set(sp1,'visible', 'off','ylim',[0 totalRastRows],'xlim',cxlm);

% create the output struct
handles.raster=sp1;
handles.psth=sp2;
handles.lines=h;
handles.figure=hf;







