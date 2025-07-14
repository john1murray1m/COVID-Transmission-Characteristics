%% determine sample-time matches for individual GISAID seqs within 2 of the Date 
% and also return the prob dependent on match

% this only matches by time since other data is missing from most states
function idN=prob_link_state(T_Gst,T_Nst,p_age,p_gender,p_time_diff)

% T_Nst is the case data restricted to the particular state
% now idN has 3 columns: 1) the N ids, 2) the time diffs, 3) the probs of
% the match

% set a min alloed p value (also gets rid of cases where SEX is unknown
% (128 cases) or age is unknown (2)
minp=1e-8;

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

ss=T_Nst1.SEX;
psex=zeros(size(p_time)); 
for k=1:3
    iis=ss==gender_names{k};
    if any(iis)
        psex(iis)=p_gender(k,2);
    end
end
idN(:,3)=p_time.*p_age_state.*psex;  

ii1=idN(:,3)>=minp;
idN=idN(ii1,:);

[~,isort]=sortrows(idN,3,"descend");
idN=idN(isort,:); % sort in order of decreasing prob

