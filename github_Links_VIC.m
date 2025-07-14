%% determine link sets for the VIC data


%% the sex and age of many seqs are included for VIC

ii0=(TG_Delta.Division=="Victoria");

TG_VIC=TG_Delta(ii0,:);

% the LQD cases
ii1=T_Delta.STATE=="VIC";
T_Delta_VIC=T_Delta(ii1,:);

p_match_VIC=cell(height(TG_VIC),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 

% if removing age from links prob
    p_noage=p_age;
    p_noage(:,2)=1/height(p_age);
% if removing gender from links prob
p_nogender=p_gender;
    p_nogender(:,2)=1/height(p_gender);

for i=1:height(TG_VIC)
    T_Gst=TG_VIC(i,:);
    idG=TG_VIC{i,'IDG'};
    idN=prob_link_VIC_noany(T_Gst,T_Delta_VIC,p_time_diff);

    p_match_VIC(i,1)={idG}; 
    p_match_VIC(i,2)={idN}; 
end

save('save_links_VIC_noany',"p_match_VIC")