% 
% 
% function [out, h] = inferTuning(x,y,fr,varargin)
% 
% % varargin default values (varargin is a struct with the following possible fields)
% xwidth      = max(abs(x))+1;
% ywidth      = max(abs(y))+1;
% filtsize    = [3 3];
% filtsigma   = 0.5;
% plotflag    = 1;
% 
% Pfields = {'xwidth','ywidth','filtsize','filtsigma','plotflag'};
% for i = 1:length(Pfields) % if a params structure was provided as an input, change the requested fields
%     if ~isempty(varargin)&&isfield(varargin{1}, Pfields{i}), eval(sprintf('%s = varargin{1}.(Pfields{%d});', Pfields{i}, i)); end
% end
% if ~isempty(varargin)  % if there is a params input
%     fnames = fieldnames(varargin{1}); % cycle through field names and make sure all are recognized
%     for i = 1:length(fnames)
%         recognized = max(strcmp(fnames{i},Pfields));
%         if recognized == 0, fprintf('fieldname %s not recognized\n',fnames{i}); end
%     end
% end
% 
% 
% %% build the image matrix
% % make sure the data fit within the bounds
% bound = y>-ywidth & y<ywidth & x>-xwidth & x<xwidth;
% y=y(bound);
% x=x(bound);
% fr=fr(bound);
% 
% % create x and y vectors
% xvec = -xwidth:xwidth;
% yvec = -ywidth:ywidth;
% 
% % create matrix of saccade/firing rates
% rallx=round(x)+round(xwidth);
% rally=round(y)+round(ywidth);
% act_map=zeros(length(yvec),length(xvec));
% act_map(sub2ind(size(act_map),rally,rallx))=ones(1,length(rallx)).*fr;
% 
% %% filtering
% % Create the gaussian filter for smoothing
% G = fspecial('gaussian',filtsize,filtsigma);
% % Filter it
% Ig = imfilter(act_map,G,'same');
% 
% %% create output variable
% out.image = Ig;
% out.x = xvec;
% out.y = yvec;
% 
% 
% %% plotting
% if plotflag
% h.f = figure();
% h.im = imagesc(xvec,yvec,Ig);
% else h=NaN;
% end




% this script is still unfinished


ind=1;
for trial = 1:length(withSaccs)
    curtrial=withSaccs(trial);
    for saccade_num = 1:length(Trials(curtrial).Saccades)
       
        sx1 = Trials(curtrial).Saccades(saccade_num).horizontal_position_start;
        sy1 = Trials(curtrial).Saccades(saccade_num).vertical_position_start;
        sx2 = Trials(curtrial).Saccades(saccade_num).horizontal_position_end;
        sy2 = Trials(curtrial).Saccades(saccade_num).vertical_position_end;

        
        cx2 = sx2-sx1;
        cy2 = sy2-sy1;

        SExy(ind,:) = [cx2 cy2];
        ind = ind + 1;
    end
end


figure(1)
binsize=5;
xb=-55:binsize:55;
yb=xb;
n=length(xb);

% 
xr = interp1(xb, 0.5:numel(xb)-0.5, SExy(:,1), 'nearest');
yr = interp1(xb, 0.5:numel(yb)-0.5, SExy(:,2), 'nearest');

% remove any NaNs (for accumarray)
bad = find(isnan(yr) | isnan(xr));
xr(bad)=[];yr(bad)=[];

% build the histogram
Z = accumarray([yr xr] + 0.5, 1, [n n]);

% plotting
surf(xb, yb, Z);
xlabel('x');
ylabel('y');
zlabel('count');
view(2)
