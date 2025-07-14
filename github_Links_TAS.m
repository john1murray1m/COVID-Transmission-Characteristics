%% link cases in TAS

ii0=(TG_Delta.Division=="Tasmania");

TG_TAS=TG_Delta(ii0,:);

% the LQD cases
ii1=T_Delta.STATE=="TAS";
T_Delta_TAS=T_Delta(ii1,:);

gender_names={'Female','Male','X'};

mindiff=p_time_diff(1,1); % this is -2 at the moment


p_match_TAS=cell(height(TG_TAS),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 





%% case 6
% neither gender nor age are known-
% if removing age from links prob
    p_noage=p_age;
    p_noage(:,2)=1/height(p_age);
% if removing gender from links prob
p_nogender=p_gender;
    p_nogender(:,2)=1/height(p_gender);

for i=1:height(TG_TAS)
    T_Gst=TG_TAS(i,:);
    idG=TG_TAS{i,'IDG'};
    idN=prob_link_state(T_Gst,T_Delta_TAS,p_noage,p_nogender,p_time_diff);

    p_match_TAS(i,1)={idG}; 
    p_match_TAS(i,2)={idN}; 
end

save('save_links_TAS_noany',"p_match_TAS")