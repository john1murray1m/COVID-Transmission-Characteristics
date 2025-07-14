%% create linkage sets for ACT

ii0=(TG_Delta.Division=="Australian Capital Territory");

TG_ACT=TG_Delta(ii0,:);

% the LQD cases
ii1=T_Delta.STATE=="ACT";
T_Delta_ACT=T_Delta(ii1,:);

gender_names={'Female','Male','X'};

mindiff=p_time_diff(1,1); % this is -2 at the moment


p_match_ACT=cell(height(TG_ACT),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 





%% case 6
% neither gender nor age are known-
% if removing age from links prob
    p_noage=p_age;
    p_noage(:,2)=1/height(p_age);
% if removing gender from links prob
p_nogender=p_gender;
    p_nogender(:,2)=1/height(p_gender);


for i=1:height(TG_ACT)
    T_Gst=TG_ACT(i,:);
    idG=TG_ACT{i,'IDG'};
    idN=prob_link_state(T_Gst,T_Delta_ACT,p_noage,p_nogender,p_time_diff);

    p_match_ACT(i,1)={idG}; 
    p_match_ACT(i,2)={idN}; 
end
save('save_links_ACT_noany',"p_match_ACT")