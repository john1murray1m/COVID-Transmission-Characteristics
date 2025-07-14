%% determine sample-time matches for individual South Australian GISAID seqs 
% these all include sex
% and also return the prob dependent on match

% this only matches by time since other data is missing from most states
function idN=prob_link_state_SA(T_Gst,T_Nst,p_age,p_time_diff)

% T_Nst is the case data restricted to the particular state
% now idN has 3 columns: 1) the N ids, 2) the time diffs, 3) the probs of
% the match


gender_names={'Female','Male','X'};


if height(T_Gst)>1
    error('There should only be a single G seq')
end

mindiff=p_time_diff(1,1); % this is -2 at the moment
maxdiff=p_time_diff(end,1); % this is 1 at the moment


dateG=T_Gst.Date;
ddays=days(T_Nst.SPECIMEN_DATE-dateG);
ii1=ddays>=mindiff & ddays<=maxdiff;
ddays=ddays(ii1);

T_Nst1=T_Nst(ii1,:);

% the same gender
ii1=T_Nst1.SEX==T_Gst.Sex;
if any(ii1)
    T_Nst1=T_Nst1(ii1,:);
    ddays=ddays(ii1);
end

idN=NaN(height(T_Nst1),3); % the matrix of N ids and their date diffs
idN(:,1)=T_Nst1.ID;
idN(:,2)=ddays;

% calculate their probs dependent on the matching mechanism
p_time=p_time_diff(idN(:,2)-mindiff+1,2);

ages=T_Nst1.AGE;

  % if any of the ages are above max in prob then replace with max
 age_max=p_age(end,1);
 ii1=ages>age_max;
 if any(ii1)
     ages(ii1)=age_max;
 end
    p_age_state=p_age(ages+1,2);

idN(:,3)=p_time.*p_age_state;  

[~,isort]=sortrows(idN,3,"descend");
idN=idN(isort,:); % sort in order of decreasing prob

