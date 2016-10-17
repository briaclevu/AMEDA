function mod_merging_splitting(name)
%mod_merging_splitting(tracks)
%
%  Resolve merging and splitting interaction from trajectories
%
%  Find indice of interaction eddy as 'interaction' field
%  Flag the inetraction as 'merge' (=1) or/and 'split' (=1)
%  Remove eddies associated with their next interaction
%  which are shorter than cut_off (save in 'short').

%  - cut_off (days) are the minimum duration tracks recorded
%       use cut_off=0 to use the explicit turnover time (tau) for each eddy
%       use cut_off=1 to keep all eddies
%  - Dt (days) is the tolerance of time steps after eddy disapears
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
diary([path_out,'log_eddy_merging_splitting_',name,'.txt']);

%----------------------------------------------------------
% load unfiltrated eddy tracks

disp('Load tracks ...')

load([path_out,'eddy_tracks_',name],'tracks')
stepF = max(tracks(end).step);

%----------------------------------------------------------
% set name indice belonging to shapes2
tracks_name = fieldnames(tracks);
n1 = find(strcmp(tracks_name,'x2'));
n2 = find(strcmp(tracks_name,'shapes2'));

%----------------------------------------------------------
% initialize tracks interaction detection and type flag
for j1=1:length(tracks)
    
    tracks(j1).split = nan(length(tracks(j1).step),1);
    tracks(j1).merge = nan(length(tracks(j1).step),1);
    tracks(j1).split2 = nan(length(tracks(j1).step),1);
    tracks(j1).merge2 = nan(length(tracks(j1).step),1);
    
end

%----------------------------------------------------------
% flag merging & spliting period of the interacting eddies

disp(['Find merging and splitting of the interacting eddies '...
    'among the ',num2str(length(tracks)),' eddies ...'])

short = false(1,length(tracks));

for j=1:length(tracks)
    
    % flag first interaction
    ind = ~isnan(tracks(j).interaction);
    Nind = tracks(j).interaction(ind);
    
    % flag as splitting if eddy starts after step(i)-Dt
    tracks(j).split(ind & tracks(j).step - Dt <= tracks(j).step(1)) = 1;
    tracks(j).split(ind & tracks(j).step - Dt > tracks(j).step(1)) = 0;
    % flag as merging if eddy finishes before step(i)+Dt
    tracks(j).merge(ind & tracks(j).step + Dt >= tracks(j).step(end)) = 1;
    tracks(j).merge(ind & tracks(j).step + Dt < tracks(j).step(end)) = 0;
    
    % flag second interaction
    ind2 = ~isnan(tracks(j).interaction2);
    Nind2 = tracks(j).interaction2(ind2);
    
    % flag as splitting if eddy starts after step(i)-Dt
    tracks(j).split2(ind2 & tracks(j).step - Dt <= tracks(j).step(1)) = 1;
    tracks(j).split2(ind2 & tracks(j).step - Dt > tracks(j).step(1)) = 0;
    % flag as merging if eddy finishes before step(i)+Dt
    tracks(j).merge2(ind2 & tracks(j).step + Dt >= tracks(j).step(end)) = 1;
    tracks(j).merge2(ind2 & tracks(j).step + Dt < tracks(j).step(end)) = 0;
    
    % first and last detection of the eddy j
    delta = [tracks(j).step(1) tracks(j).step(end)];
    tau = tracks(j).tau1;
    
    % increase detection limit with merged or spitted eddies
    if nanmax(tracks(j).merge)==1 || nanmax(tracks(j).split)==1 ||...
            nanmax(tracks(j).merge2)==1 || nanmax(tracks(j).split2)==1

        % change first detection in case of splitting
        if nanmax(tracks(j).split)==1
            delta(1) = min(tracks(Nind(1)).step(1),delta(1));
            if Nind(1)~=Nind(end)
                tau = [tracks(Nind(1)).tau1;tau];
            end
        end
        if nanmax(tracks(j).split2)==1
            delta(1) = min(tracks(Nind2(1)).step(1),delta(1));
            if Nind2(1)~=Nind2(end)
                tau = [tracks(Nind2(1)).tau1;tau];
            end
        end

        % change last detection in case of merging
        if nanmax(tracks(j).merge)==1
            delta(2) = max(tracks(Nind(end)).step(end),delta(2));
            tau = [tau;tracks(Nind(end)).tau1];
        end
        if nanmax(tracks(j).merge2)==1
            delta(2) = max(tracks(Nind2(end)).step(end),delta(2));
            tau = [tau;tracks(Nind2(end)).tau1];
        end
    end
    
    % taking merging and spltting flag eddies
    % shorter than their turnover time or the cut_off time
    if cut_off==0 && ( stepF - delta(1) )*dps >= mean(tau) + Dt
        
        short(j) = (diff(delta) + 1)*dps < mean(tau);
        
    elseif ( stepF - delta(1) )*dps >= cut_off + Dt
        
        short(j) = (diff(delta) + 1)*dps < cut_off;
        
    end
    
end
        
disp(' ')

%----------------------------------------------------------
% filter tracks shorter than eddy turnover time or cut_off days

if cut_off==0
    disp('Filter eddies shorter than their turnover time ...')
else
    disp(['Filter eddies shorter than ',num2str(cut_off),' days ...'])
end

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
        tracks(j).split(ind3) = tracks(j).interaction2(ind3);
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
% save tracks after merging spltting and cut_off

disp(['Save tracks2 in eddy_tracks_',name,'.mat ...'])

tracks2=tracks;

save([path_out,'eddy_tracks_',name],'tracks2','short','-append')

diary off
