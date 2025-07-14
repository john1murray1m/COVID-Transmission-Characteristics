%% incorporate geo data into the NSW links


p_cum_geo(2,:)=[]; % get rid of the second 0 entry

% save('save_NSW_lab_locations',"lab_LL")
load('save_NSW_lab_locations')

% %   save('save_Westmead','p_match_Westmead',"families",'p_match_children')
% % load('save_Westmead')
load('save_Westmead_noany')

p_match_NSW=cell(height(p_match_NSW_nogeo),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 

p_no_geo=interp1(p_cum_geo(:,1),p_cum_geo(:,2),50); % if there is no geo info 
%                                   then set the prob for 50kms (=0.0358)
% % if testing the influence of Histopath lacking locations on NSW distances
% p_no_geo=p_no_geo/10;

ii1=TG_Delta.Division=='New South Wales';
TG_NSW=TG_Delta(ii1,:);

lab_NSW=unique(TG_NSW.Laboratory);

westG=cell2mat(p_match_Westmead(:,1)); % the idg of Westmead

for i=1:height(p_match_NSW_nogeo)
    idG=p_match_NSW_nogeo{i,1};
    p_match_NSW(i,1)={idG}; 

    T_Gst=TG_Delta(idG,:);
    labi=T_Gst.Laboratory;
    ilab=find(labi==lab_NSW);
    if ilab~=34 % Westmead
        idN=p_match_NSW_nogeo{i,2};
        if ~isempty(idN)      
            T_Nst=T_Delta(idN(:,1),:);
            if ~isempty(lab_LL{ilab,2}) % there is some geo info for this lab
                idN1=prob_link_NSW_geo(T_Gst,T_Nst,idN,p_cum_geo,lab_LL{ilab,2},lab_LL{ilab,3});
                p_match_NSW(i,2)={idN1}; 
            else
                idN(:,3)=idN(:,3)*p_no_geo; % reduce the prob to match with additional geo data changes
                if height(idN)>100
                    idN=idN(1:100,:); % restrict to the top 100, they were sorted in the nogeo link
                end
                p_match_NSW(i,2)={idN}; 
            end
        end

    else % swap with Westmead match
        ii2=find(westG==idG);
        p_match_NSW(i,2)=p_match_Westmead(ii2,2); % sorted and restricted to top 100
    end
end


save('save_links_NSW_noany',"p_match_NSW") 

            
