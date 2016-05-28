% [out, h] = inferTuning(x,y,fr,varargin)
%
% Inputs    x: vector of x-values.
%           y: vector of y-values.
%           fr: vector of instantaneous or average firing rates
%           varargin: a struct with optional parameters
%                 - xwidth: the width of the plot window (default is to
%                 include all data points)
%                 - ywidth: the height of the plot window (default is to
%                 include all data points)
%                 - filtsize: size of the gaussian filter for smoothing
%                 (default: [3 3])
%                 - filtsigma: standard deviation of the gaussian filter
%                 (default: 0.5)
%                 - plotflag: plot or not (default is yes)  
%
% Output    out: a struct with the image matrix, x, and y
%           h: handle to the axes and figure (if plotting is turned on).
%
% DKW, Jan 2016

% to do: 
% * add thresholding in order to estimate a peak and return info about that
% peak (e.g., centroid, ellipse object fit, etc.)

function [out, h] = inferTuning(x,y,fr,varargin)

% varargin default values (varargin is a struct with the following possible fields)
xwidth      = max(abs(x))+1;
ywidth      = max(abs(y))+1;
filtsize    = [2 2];
filtsigma   = 0.5;
plotflag    = 1;
fig_Handle   = [];
axes_Handle  = [];

Pfields = {'xwidth','ywidth','filtsize','filtsigma','plotflag','fig_Handle','axes_Handle'};
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


%% build the image matrix
% make sure the data fit within the bounds
bound = y>-ywidth & y<ywidth & x>-xwidth & x<xwidth;
y=y(bound);
x=x(bound);
fr=fr(bound);

% create x and y vectors
xvec = -xwidth:xwidth;
yvec = -ywidth:ywidth;
act_map=NaN(length(yvec),length(xvec));



%% smoothing

if 1
% use k-nearest neighbors smoothing algorithm from Josh
% create matrix of saccade/firing rates
rallx=ceil(x);
rally=ceil(y);
% get x,y position of each pixel on the grid
[xg,yg]=meshgrid(xvec,yvec);
knn=round(length(fr)*.1); % number of neighbors that will influence each other
decay=-.7; % influence of neighbors decays over distance
peaks = rf_smooth([rallx' rally'],[xg(:) yg(:)],fr,decay,knn);
Ig=reshape(peaks,size(act_map));
end


if 0
% Convolve with gaussian filter, ignore NaNs
rallx=ceil(x)+round(xwidth);
rally=ceil(y)+round(ywidth);
act_map(sub2ind(size(act_map),rally,rallx))=ones(1,length(rallx)).*fr;
G = fspecial('gaussian',filtsize,filtsigma);
% Filter it
Ig = nanconv(act_map,G);
end




%% create output variable
out.image = Ig;
out.x = xvec;
out.y = yvec;


%% plotting
if plotflag
    if isempty(fig_Handle)
        h.f = figure();
    else
        h.f = fig_Handle;
        figure(h.f);
    end
    if ~isempty(axes_Handle)
        subplot(axes_Handle)
        imagesc(xvec,yvec,Ig);
%         set(gca,'Ydir','normal')
        h.im = axes_Handle;
    else
        h.im = imagesc(xvec,yvec,Ig);
    end

    hold on
    plot([0 0],ylim,'w--')
    plot(xlim,[0 0],'w--')
else h=NaN;
end