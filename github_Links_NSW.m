%% determine initial link sets for the NSW data without geo data
% The geo data and associated prob is determined later


%% the sex and age of many seqs are included for NSW as well as some geo info

ii0=(TG_Delta.Division=="New South Wales");

TG_NSW=TG_Delta(ii0,:);

% the NSW cases
ii1=T_Delta.STATE=="NSW";
T_Delta_NSW=T_Delta(ii1,:);

p_match_NSW_nogeo=cell(height(TG_NSW),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 


for i=1:height(TG_NSW)
    T_Gst=TG_NSW(i,:);
    idG=TG_NSW{i,'IDG'};
    % if eliminating age from linkage prob
    p_noage=p_age;
    p_noage(:,2)=1/height(p_age);
    idN=prob_link_NSW_noany(T_Gst,T_Delta_NSW,p_time_diff);

    p_match_NSW_nogeo(i,1)={idG}; 
    p_match_NSW_nogeo(i,2)={idN}; 
end

save('save_links_NSW_nogeo_noany',"p_match_NSW_nogeo")