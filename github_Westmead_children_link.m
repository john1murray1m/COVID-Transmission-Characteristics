%% estimate the cases linked to submissions by the Children's Hospital Westmead
%23/10/23
% 8/4/24 Modify to exclude prob related to age gender except proference for
% children

% % save('prob_distance_save','p_cum_geo','lab_names') % this is the prob of this dist
% or more extreme

% load('link_probs')
% 
% load('prob_distance_save')

ii1=(TG_Delta.Laboratory==lab_names(41));

TG_Children=TG_Delta(ii1,:);

% the  cases
ii1=T_Delta.STATE=="NSW";
T_Delta_NSW=T_Delta(ii1,:);

p_match_children=cell(height(TG_Children),2); % the first is the cell of N ID,
            % the second is the number of samples on that date so a max of
            % family size

p_match_Westmead=cell(height(TG_Children),2); % the first is the G ID,
            % the second is IDN 

tab_fam=tabulate(TG_Children.Date);

for i=1:height(TG_Children)
    T_Gst=TG_Children(i,:);
    p_match_Westmead(i,1)={T_Gst.IDG};

    ii2=T_Gst.Date==tab_fam(:,1);
    p_match_children(i,2)=tab_fam(ii2,2);
    idN=prob_link_child_noany(T_Gst,T_Delta_NSW,p_time_diff,p_cum_geo);
    if ~isempty(idN)
        p_match_children(i,1)={idN}; % 1: NID, 2: time_diff, 3: prob, 4: postcode
                                        % 5: ages, 6: distance to Westmead
                                        % sorted by descending prob,
        p_match_Westmead(i,2)={idN}; % 1: NID, 2: time_diff, 3: prob, 4: postcode
                                    % 5: ages, 6: distance to Westmead
    end
end


    