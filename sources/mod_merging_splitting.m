function mod_merging_splitting(name)
%mod_merging_splitting(name)
%
%  Resolve merging and splitting interaction from trajectories
%  and apply a time filtering. Struture of tracks2 is the same as tracks
%  with addtiional flag about interaction, merging and splitting.
%
%  Fill indices of interaction eddy as 'interaction' field
%  and flag the interaction as 'merge' (=1) or/and 'split' (=1)
%  Remove eddies associated which are shorter than cut_off (save in 'short').

%  - cut_off (days) are the minimum duration tracks recorded
%       use cut_off=0 to use the explicit turnover time (tau) for each eddy
%       use cut_off=1 to keep all eddies
%  - Dt (days) is the tolerance (in time steps) to recover the detection
%       after an eddy disapeared
%  - T (seconds) is the period for one day
%  - dps is the number of days between 2 time step (dps=1/8 gives 3h)
%
%  For a description of the input parameters see param_eddy_tracking.m
%
%-------------------------
%   Ver. 3.2 Jul.2016 Briac Le Vu
%-------------------------
%
%=========================

% load key_source and parameters (use mod_eddy_params.m first)
%----------------------------------------------

load('param_eddy_tracking')

% begin the log file
%----------------------------------------------
diary([path_out,'log_eddy_merging_splitting',name,'.txt']);

%----------------------------------------------------------
% load unfiltrated eddy tracks

disp('Load tracks ...')

load([path_out,'eddy_tracks',name],'tracks','warn_tracks')

stepI = min(tracks(1).step);
stepF = max(tracks(end).step);

%----------------------------------------------------------
% set name indice belonging to shapes2
tracks_name = fieldnames(tracks);
n1 = find(strcmp(tracks_name,'x2'));
n2 = find(strcmp(tracks_name,'shapes2'));

%----------------------------------------------------------
% flag merging & spliting period of the interacting eddies

disp(['Find merging and splitting of the interacting eddies '...
    'among the ',num2str(length(tracks)),' eddies ...'])

for j=1:length(tracks)
   
    % initialize tracks interaction detection and type flag
    tracks(j).split = nan(length(tracks(j).step),1);
    tracks(j).merge = nan(length(tracks(j).step),1);
    tracks(j).split2 = nan(length(tracks(j).step),1);
    tracks(j).merge2 = nan(length(tracks(j).step),1);
    
    % flag first possible interaction
    ind  = ~isnan(tracks(j).interaction);
    indm = find(ind==1,1,'last');
    inds = find(ind==1,1,'first');

    % flag as splitting if eddy starts after step(i)-Dt/2
    tracks(j).split(ind & tracks(j).step - Dt/2 <= tracks(j).step(1)) = 0;
    if tracks(j).step(inds) - Dt/2 <= tracks(j).step(1)
        tracks(j).split(inds) = 1;
    end
    
    % flag as merging if eddy finishes before step(i)+Dt/2
    tracks(j).merge(ind & tracks(j).step + Dt/2 >= tracks(j).step(end)) = 0;
    if tracks(j).step(indm) + Dt/2 >= tracks(j).step(end)
        tracks(j).merge(indm) = 1;
    end
    
    % flag second possible interaction
    ind2 = ~isnan(tracks(j).interaction2);
    indm = find(ind2==1,1,'last');
    inds = find(ind2==1,1,'first');

    % flag as splitting if eddy starts after step(i)-Dt/2
    tracks(j).split2(ind2 & tracks(j).step - Dt/2 <= tracks(j).step(1)) = 0;
    if tracks(j).step(inds) - Dt/2 <= tracks(j).step(1)
        tracks(j).split2(inds) = 1;
    end
    % flag as merging if eddy finishes before step(i)+Dt/2
    tracks(j).merge2(ind2 & tracks(j).step + Dt/2 >= tracks(j).step(end)) = 0;
    if  tracks(j).step(indm) + Dt/2 >= tracks(j).step(end)
        tracks(j).merge2(indm) = 1;
    end
    
end

tracks_name = fieldnames(tracks);
warn_name = fieldnames(warn_tracks);

%----------------------------------------------------------
% Rearrange tracks to prevent lost of tracks due to simultaneaous
% splitting-merging process (may be due to a bad spatiotemporal satellite
% coverage)

disp('Extend tracks when merging with its child eddy just splitted')

for j=1:length(tracks)

    trymoreconcat = 1;
    
    while trymoreconcat
        
        ind = tracks(j).merge==1;
        Nind = tracks(j).interaction(ind);
        
        trymoreconcat = 0;
        
        if ~isempty(Nind) && ~isnan(Nind)
            inds = tracks(Nind).split==1;
            Ninds = tracks(Nind).interaction(inds);

            % find case where splitted eddy is merging with its parent
            % inside the Dt period
            if ~isempty(Ninds) && ~isnan(Ninds)
                
                if Ninds==j && tracks(j).step(ind) < tracks(Nind).step(inds) + Dt &&...
                    tracks(j).step(ind) >= tracks(Nind).step(inds)

                    Tind = tracks(j).step(end);
                    Tind1 = tracks(Nind).step > Tind;
                    Tind2 = tracks(Nind).step <= Tind;

                    if ~isempty(Tind1)

                        % record warn_flag
                        warn_tracks(j).concat(tracks(j).merge==1) = 1;
                        warn_tracks(j).concat(tracks(j).merge==0) = 0;

                        % remove merging flag
                        Mstp = tracks(j).step(tracks(j).merge==1);
                        tracks(j).merge(tracks(j).merge==1)=0;

                        % concatenate the part of the splitted eddy following the
                        % merging
                        for n=1:length(tracks_name)
                            tracks(j).(tracks_name{n}) =...
                                [tracks(j).(tracks_name{n}); tracks(Nind).(tracks_name{n})(Tind1)];
                        end

                        % modify the reference to the splitted eddy in interaction
                        % with other eddies
                        Iind = ~isnan(tracks(Nind).interaction);
                        Iind = unique(tracks(Nind).interaction(Iind));
                        for n=1:length(Iind)
                            tracks(Iind(n)).interaction(tracks(Iind(n)).interaction==Nind &...
                                tracks(Iind(n)).step > Tind)=j;
                        end
                        
                        % add merging flag for the child
                        %tracks(Nind).merge(tracks(Nind).step==Mstp) = 1;

                        % remove this part from the splitted eddy
                        for n=1:length(tracks_name)
                            tracks(Nind).(tracks_name{n}) =...
                                tracks(Nind).(tracks_name{n})(Tind2);
                        end

                        % the concatenate eddy may merge at the end of its life
                        trymoreconcat=1;
                    end
                end
            end
        end
    end
    
    %         tracks(j) = [tracks(j) tracks(Nind)(Tind)];
    %     ind2 = ~isnan(tracks(j).interaction2);
    %     Nind2 = tracks(j).interaction2(ind2);
    

end

%----------------------------------------------------------
% Record robust eddies as eddies tracked long enough to rotate typical
% twice (cut_off=0).

if cut_off==0
    disp('Filter eddies shorter than their turnover time ...')
else
    disp(['Filter eddies shorter than ',num2str(cut_off),' days ...'])
end

short = false(1,length(tracks));

for j=1:length(tracks)

    % first and last detection of the eddy j
    delta = [tracks(j).step(1) tracks(j).step(end)];
    tau = tracks(j).tau1;
    
    % taking merging and spltting flag eddies
    % shorter than twice their turnover time or the cut_off time
    if isnan(nanmean(tau))
        short(j) = true;
    elseif cut_off==0
        short(j) = (diff(delta) + 1)*dps < 2 * nanmean(tau);
    elseif cut_off>0
        short(j) = (diff(delta) + 1)*dps < cut_off;
    end
    
    % and keep the last and first tracks for NRT analysis
    if ~isnan(nanmean(tau)) && nrt
        if cut_off==0 && ( stepF - delta(1) + 1 )*dps < 2 * nanmean(tau) ||...
                cut_off==0 && ( delta(end) - stepI + 1 )*dps < 2 * nanmean(tau) ||...
                cut_off>0 && ( stepF - delta(1) +1 )*dps < cut_off ||...
                cut_off>0 && ( delta(end) - stepI + 1 )*dps < cut_off
            short(j) = false;
        end 
    end

end
        
disp(' ')

tracks(short) = [];
Cshort = cumsum(short);

disp([' ',num2str(length(tracks)),' eddies remaining'])

%----------------------------------------------------------
% adjust interaction indice after filtrering

disp(' adjust interaction indice and clean shapes2 ...')

moved = false(1,length(tracks));
moved2 = false(1,length(tracks));

for j=1:length(tracks)
    
    % find interaction indices
    ind = find(~isnan(tracks(j).interaction));
    
    if any(ind)
        % interacing eddies
        Nind = tracks(j).interaction(ind);

        % adjust tracking indices
        tracks(j).interaction(ind) = Nind - Cshort(Nind)';
        
        if any(short(Nind))
            % remove interaction if interacting eddy is missing
            tracks(j).interaction(ind(short(Nind))) = NaN;
            tracks(j).split(ind(short(Nind))) = NaN;
            tracks(j).merge(ind(short(Nind))) = NaN;
            
            Uind = unique(Nind(short(Nind))' - Cshort(Nind(short(Nind))));
            disp(['  the interacting eddies ',num2str(Uind),...
                ' has been filtred for eddy ',num2str(j),...
                ' on step ', num2str(tracks(j).step(ind(short(Nind)))')])
            
            moved(j) = 1;
        end
    end
    
    % find interaction2 indices
    ind2 = find(~isnan(tracks(j).interaction2));
    
    if any(ind2)
        % interacting eddies
        Nind2 = tracks(j).interaction2(ind2);

        % adjust tracking indices
        tracks(j).interaction2(ind2) = Nind2 - Cshort(Nind2)';
    
        if any(short(Nind2))
            % remove interaction if interacting eddy is missing
            tracks(j).interaction2(ind2(short(Nind2))) = NaN;
            tracks(j).split2(ind2(short(Nind2))) = NaN;
            tracks(j).merge2(ind2(short(Nind2))) = NaN;
            
            Uind2 = unique(Nind2(short(Nind2))' - Cshort(Nind2(short(Nind2))));
            disp(['  the interacting eddies ',num2str(Uind2),...
                ' has been filtred for eddy ',num2str(j),...
                ' on step ', num2str(tracks(j).step(ind2(short(Nind2)))')])

            moved2(j) = 1;
        end
    end
    
    % find interaction2 with no interaction
    ind3 = ~isnan(tracks(j).interaction2) & isnan(tracks(j).interaction);
    
    if any(ind3)
        % move interaction2 to interaction
        tracks(j).interaction(ind3) = tracks(j).interaction2(ind3);
        tracks(j).split(ind3) = tracks(j).split2(ind3);
        tracks(j).merge(ind3) = tracks(j).merge2(ind3);
    end
    
    % check interaction2 empty
    ind4 = find(~isnan(tracks(j).interaction2));
    
    if any(ind4)
        for jj=1:length(ind4)
            disp(['  eddy ',num2str(j),' is still interacting with 2 eddies ',...
                num2str([tracks(j).interaction(ind4(jj)) tracks(j).interaction2(ind4(jj))]),...
                ' on step ', num2str(tracks(j).step(ind4(jj))')])
        end
    end
    
end

disp(' ')

%----------------------------------------------------------
% clean double shape or double center no involved in interaction

disp(' then clean double shapes with no interaction ...')

for j=1:length(tracks)
    
    if moved(j)
        for i=1:length(tracks(j).step)

            % check the interacting eddy of eddy j and the step stp(i)
            if isnan(tracks(j).interaction(i))
                % no interaction
                tracks(j).shapes2{i} = NaN;
                for n=[n1:n1+3 n2+1:n2+8]
                    tracks(j).(tracks_name{n})(i) = NaN;
                end
            end
            
        end
    end
    
    if moved2(j)
        for i=1:length(tracks(j).step)
            
            % check the interacting eddy of eddy j and the step stp(i)
            if isnan(tracks(j).interaction2(i))
                % no interaction
                tracks(j).shapes2{i} = NaN;
                for n=[n1:n1+3 n2+1:n2+8]
                    tracks(j).(tracks_name{n})(i) = NaN;
                end
            end
            
        end
    end
        

end % end loop on tracks j

disp(' ')

%----------------------------------------------------------
% Reflag merging & spliting period of the interacting eddies

disp(['Flag merging and splitting of the interacting eddies '...
    'among the ',num2str(length(tracks)),' remaining eddies ...'])

for j=1:length(tracks)
   
    % initialize tracks interaction detection and type flag
    tracks(j).split = nan(length(tracks(j).step),1);
    tracks(j).merge = nan(length(tracks(j).step),1);
    tracks(j).split2 = nan(length(tracks(j).step),1);
    tracks(j).merge2 = nan(length(tracks(j).step),1);
    
    % flag first possible interaction
    ind  = ~isnan(tracks(j).interaction);
    indm = find(ind==1,1,'last');
    inds = find(ind==1,1,'first');

    % flag as splitting if eddy starts after step(i)-Dt/2
    tracks(j).split(ind & tracks(j).step - Dt/2 <= tracks(j).step(1)) = 0;
    if tracks(j).step(inds) - Dt/2 <= tracks(j).step(1)
        tracks(j).split(inds) = 1;
    end
    
    % flag as merging if eddy finishes before step(i)+Dt/2
    tracks(j).merge(ind & tracks(j).step + Dt/2 >= tracks(j).step(end)) = 0;
    if tracks(j).step(indm) + Dt/2 >= tracks(j).step(end)
        tracks(j).merge(indm) = 1;
    end
    
    % flag second possible interaction
    ind2 = ~isnan(tracks(j).interaction2);
    indm = find(ind2==1,1,'last');
    inds = find(ind2==1,1,'first');

    % flag as splitting if eddy starts after step(i)-Dt/2
    tracks(j).split2(ind2 & tracks(j).step - Dt/2 <= tracks(j).step(1)) = 0;
    if tracks(j).step(inds) - Dt/2 <= tracks(j).step(1)
        tracks(j).split2(inds) = 1;
    end
    % flag as merging if eddy finishes before step(i)+Dt/2
    tracks(j).merge2(ind2 & tracks(j).step + Dt/2 >= tracks(j).step(end)) = 0;
    if  tracks(j).step(indm) + Dt/2 >= tracks(j).step(end)
        tracks(j).merge2(indm) = 1;
    end
    
end

%----------------------------------------------------------
% save tracks after merging spltting and cut_off

disp(['Save tracks2 in eddy_tracks2',name,'.mat ...'])

tracks2=tracks;

warn_tracks2=warn_tracks;

save([path_out,'eddy_tracks2',name],'tracks2','warn_tracks2','short','-v7.3')

diary off
