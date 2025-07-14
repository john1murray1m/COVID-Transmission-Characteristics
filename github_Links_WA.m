%% Link WA cases

ii0=(TG_Delta.Division=="Western Australia");

TG_WA=TG_Delta(ii0,:);

% the LQD cases
ii1=T_Delta.STATE=="WA";
T_Delta_WA=T_Delta(ii1,:);

gender_names={'Female','Male','X'};

mindiff=p_time_diff(1,1); % this is -2 at the moment


p_match_WA=cell(height(TG_WA),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 





%% case 6
% neither gender nor age are known-
% if removing age from links prob
    p_noage=p_age;
    p_noage(:,2)=1/height(p_age);
% if removing gender from links prob
p_nogender=p_gender;
    p_nogender(:,2)=1/height(p_gender);

for i=1:height(TG_WA)
    T_Gst=TG_WA(i,:);
    idG=TG_WA{i,'IDG'};
    idN=prob_link_state(T_Gst,T_Delta_WA,p_noage,p_nogender,p_time_diff);

    p_match_WA(i,1)={idG}; 
    p_match_WA(i,2)={idN}; 
end

save('save_links_WA_noany',"p_match_WA")