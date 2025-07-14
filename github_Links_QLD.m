%% determine link sets for the QLD data
% 2/11/23
% run NNDSS_data_delta_input.m and Linkage_start.m first
% also run specimen_date_diffs to determine Links and probs

% 18/9/23 This differs from specimen_date_diffs since it collects N cases
% within 2 days of Date and uses p_time_diff from the non-dupl matches
% 6/10/23 Modify to call a function to extract the possible N and their
% probs dependent on match type

% save('link_probs',"p_time_diff",'p_age','p_gender',"gender_code",...
%     'p_age_grp','age_groups','match_age_sex','iidup','non_iidup','iisanon',...
%     'Link_sa')
load('link_probs')

ii0=(TG_Delta.Division=="Queensland");

TG_QLD=TG_Delta(ii0,:);

% the LQD cases
ii1=T_Delta.STATE=="QLD";
T_Delta_QLD=T_Delta(ii1,:);

gender_names={'Female','Male','X'};

mindiff=p_time_diff(1,1); % this is -2 at the moment


p_match_QLD=cell(height(TG_QLD),2); % 1) the G ID, 2) the [N.IDs, tdiffs, probs] 





%% case 6
% neither gender nor age are known-
% if removing age from links prob
    p_noage=p_age;
    p_noage(:,2)=1/height(p_age);
% if removing gender from links prob
p_nogender=p_gender;
    p_nogender(:,2)=1/height(p_gender);

for i=1:height(TG_QLD)
    T_Gst=TG_QLD(i,:);
    idG=TG_QLD{i,'IDG'};
%     idN=prob_link_state(T_Gst,T_Delta_QLD,p_age,p_gender,p_time_diff);
    idN=prob_link_state(T_Gst,T_Delta_QLD,p_noage,p_nogender,p_time_diff);

    p_match_QLD(i,1)={idG}; 
    p_match_QLD(i,2)={idN}; 
end

% save('save_links_QLD',"p_match_QLD")
% save('save_links_QLD_noage',"p_match_QLD")
save('save_links_QLD_noany',"p_match_QLD")