%% determine link sets for the SA data


%% the sex of all cases is included for SA

ii0=(TG_Delta.Division=="South Australia");

TG_SA=TG_Delta(ii0,:);

% the LQD cases
ii1=T_Delta.STATE=="SA";
T_Delta_SA=T_Delta(ii1,:);

gender_names={'Female','Male','X'};

mindiff=p_time_diff(1,1); % this is -2 at the moment


p_match_SA=cell(height(TG_SA),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 





%% case 6
% neither gender nor age are known-
% if removing age from links prob
    p_noage=p_age;
    p_noage(:,2)=1/height(p_age);

    
for i=1:height(TG_SA)
    T_Gst=TG_SA(i,:);
    idG=TG_SA{i,'IDG'};
    idN=prob_link_state_SA(T_Gst,T_Delta_SA,p_noage,p_time_diff);

    p_match_SA(i,1)={idG}; 
    p_match_SA(i,2)={idN}; 
end

save('save_links_SA_noany',"p_match_SA")