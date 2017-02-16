function PSTH = radial_PSacTH(Trials,channels,plotflag)

if ispc
    opengl software % use openGL software rather than hardware (since you are using alpha transparency and this isn't compatible with openGL hardware currently)
end

if nargin<3
    plotflag=1;
end

% temporal parameters
dt = 10;
time_before = 500;
time_after = 200;


% plot only rewarded saccades?
onlyRewFlag = 1;


trialvec = 1:length(Trials);
statsctr = 1;
for electrode = channels
    
    if length(Trials(1).Electrodes(electrode).Units)==1
        unitvec=1:1;
        unitsub=0;
    else
        unitvec=2:length(Trials(1).Electrodes(electrode).Units);
        unitsub=1;
    end
    
    for unit = unitvec % first unit is unsorted spikes
        ind = 1;
        
        for trial = 1:length(trialvec)
            curtrial=trialvec(trial);
            ecodes = [Trials(curtrial).Events(:).Code];
            
            tx = Trials(curtrial).Target.x;
            ty = Trials(curtrial).Target.y;
            
            % get saccade angle
            sacc_angle = atan2d(ty,tx);
            sacc_angle(sacc_angle<0)=sacc_angle(sacc_angle<0)+360;
            
            
            % find rewarded saccade
            saccade_num = find([Trials(curtrial).Saccades.rewarded]);
            
            
            if ~isempty(saccade_num)
                
                data_struct(ind).angle = sacc_angle;
                
                % get saccade onset time
                start_time = Trials(curtrial).Saccades(saccade_num).t_start_sacc;
                
                % get neural data
                temp = {[Trials(curtrial).Electrodes(electrode).Units(unit).Times] - double(start_time)};
                data_cells{ind}=temp{1}';
                
                ind = ind + 1;
            end
                
                
                
            
%             for saccade_num = 1:length(Trials(curtrial).Saccades)
%                 
%                 if onlyRewFlag && ~Trials(curtrial).Saccades(saccade_num).rewarded
%                     continue
%                 end
%                 
%                 % get last saccade (assuming it is the one that was rewarded)
%                 sx1 = Trials(curtrial).Saccades(saccade_num).x_sacc_start;
%                 sy1 = Trials(curtrial).Saccades(saccade_num).y_sacc_start;
%                 sx2 = Trials(curtrial).Saccades(saccade_num).x_sacc_end;
%                 sy2 = Trials(curtrial).Saccades(saccade_num).y_sacc_end;
%                 
%                 % center by x1 and y1
%                 cx2 = sx2-sx1;
%                 cy2 = sy2-sy1;
%                 
%                 % get saccade angle
%                 sacc_angle = atan2d(cy2,cx2);
%                 sacc_angle(sacc_angle<0)=sacc_angle(sacc_angle<0)+360;
%                 data_struct(ind).angle = sacc_angle;
%                 
%                 % get saccade onset time
%                 start_time = Trials(curtrial).Saccades(saccade_num).t_start_sacc;
%                 
%                 % get neural data
%                 temp = {[Trials(curtrial).Electrodes(electrode).Units(unit).Times] - double(start_time)};
%                 data_cells{ind}=temp{1}';
% 
%                 ind = ind + 1;
%             end
        end
        
        
        
        % bin by angle
        numbins=8;
        bin_edges=linspace(0,360,numbins+1);
        [n_ang,bin_edges,ang_bins]=histcounts([data_struct.angle],bin_edges);
                
        time_params(1).zero_time=0;
        time_params(1).start_time=-time_before;
        time_params(1).end_time=time_after;
        time_params(1).dt=5;
        time_params(2:numbins)=time_params(1);
        
        other_params.errBars=1;
        other_params.useSEs=1;
        other_params.smoothflag=1;
        other_params.smoothtype='gauss';
        other_params.gauss_sigma = 10;
        other_params.plotLoc = [0.55 0 .43 1];
        other_params.figHand = figure(unit-unitsub);
        other_params.plotflag = plotflag;
        
        for tb = 1:numbins
            other_params.names{tb}=[num2str(bin_edges(tb)) '-' num2str(bin_edges(tb+1))];
            all_spikes2{tb}=data_cells(ang_bins==tb);
        end
        
        if plotflag
            subplot(1,2,2)
            title('Saccade locked')
            set(gcf,'position',[206         415        1578         547])
            set(gcf,'Name',['unit ' num2str(unit-unitsub)],'NumberTitle','off')
        end
        PSTH.electrode(electrode).unit(unit).data = nda_PSTH(all_spikes2,time_params,other_params);

        
%         % create subplots
%         fh(unit)=figure();
%         other_params.figHand=fh(unit);
%         heights=1/2;
%         widths=1/(numbins/2);
%         binlefts=0:widths:(1-widths);
%         lefts=repmat(binlefts,1,2);
%         bottoms=[.5*ones(1,length(binlefts)) zeros(1,length(binlefts))];
%         for tb = 1:numbins
%             all_spikes={data_cells(ang_bins==tb)};
%             
%             other_params.plotLoc = [lefts(tb) bottoms(tb) widths heights];
%             other_params.figTitle = [num2str(bin_edges(tb)) '-' num2str(bin_edges(tb+1))];
%             PSTH_rast(all_spikes,time_params,other_params);
%         end
    end
end


