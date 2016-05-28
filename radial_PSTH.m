function radial_PSTH(Trials,channels)

if ispc
    opengl software % use openGL software rather than hardware (since you are using alpha transparency and this isn't compatible with openGL hardware currently)
end

% temporal parameters
dt = 10;
time_before = 200;
time_after = 200;



trialvec = 1:length(Trials);
statsctr = 1;

for electrode = channels
    
    for unit = 2:length(Trials(1).Electrodes(electrode).Units) % first unit is unsorted spikes
        ind = 1;
        
        for trial = 1:length(trialvec)
            curtrial=trialvec(trial);
            ecodes = [Trials(curtrial).Events(:).Code];
            
            for saccade_num = 1:length(Trials(curtrial).Saccades)
                % get last saccade (assuming it is the one that was rewarded)
                sx1 = Trials(curtrial).Saccades(saccade_num).x_sacc_start;
                sy1 = Trials(curtrial).Saccades(saccade_num).y_sacc_start;
                sx2 = Trials(curtrial).Saccades(saccade_num).x_sacc_end;
                sy2 = Trials(curtrial).Saccades(saccade_num).y_sacc_end;
                
                % center by x1 and y1
                cx2 = sx2-sx1;
                cy2 = sy2-sy1;
                
                % get saccade angle
                sacc_angle = atan2d(cy2,cx2);
                sacc_angle(sacc_angle<0)=sacc_angle(sacc_angle<0)+360;
                data_struct(ind).angle = sacc_angle;
                
                % get saccade onset time
                start_time = Trials(curtrial).Saccades(saccade_num).t_start_sacc;
                
                % get neural data
                temp = {[Trials(curtrial).Electrodes(electrode).Units(unit).Times] - double(start_time)};
                data_cells{ind}=temp{1}';

                ind = ind + 1;
            end
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
        
        for tb = 1:numbins
            other_params.names{tb}=[num2str(bin_edges(tb)) '-' num2str(bin_edges(tb+1))];
            all_spikes2{tb}=data_cells(ang_bins==tb);
        end
        PSTH_rast(all_spikes2,time_params,other_params);
        
        
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


