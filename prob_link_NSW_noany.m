%% determine sample-time matches for individual NSW GISAID seqs 
% using gender and age data when available but not geo


function idN=prob_link_NSW_noany(T_Gst,T_Nst,p_time_diff)

% now idN has 3 columns: 1) the N ids, 2) the time diffs, 3) the probs of
% the match

% set a min alloed p value (also gets rid of cases where SEX is unknown
% (128 cases) or age is unknown (2)
minp=1e-8;

if height(T_Gst)>1
    error('There should only be a single G seq')
end

minpdiff=min(p_time_diff(:,2));

mindiff=p_time_diff(1,1); % this is -2 at the moment
maxdiff=p_time_diff(end,1); % this is 1 at the moment

dateG=T_Gst.Date;
ddays=days(T_Nst.SPECIMEN_DATE-dateG);
ii1=ddays>=mindiff & ddays<=maxdiff;
ddays=ddays(ii1);


T_Nst1=T_Nst(ii1,:);

% check if sex info available for this seq
sex=T_Gst.Sex;
age=T_Gst.Age;



% restrict to the same sex if known

if sex~='unknown'
       ii1=T_Nst1.SEX==sex; % the same sex
       if any(ii1)
            T_Nst1=T_Nst1(ii1,:);
            ddays=ddays(ii1);
            psex=ones(height(T_Nst1),1); 
       else 
            psex=[];
            idN=[];
       end
else
    psex=0.5*ones(height(T_Nst1),1); 
end

if ~isempty(psex)
    
    % the case data
    ages=T_Nst1.AGE;
    sexes=T_Nst1.SEX;

    p_age_NSW=zeros(height(ages),1);     
    
    if ~ismissing(age)
           iiN1=ages==age; % the same age (the same sex enforced above if known)
           if any(iiN1)
                p_age_NSW(iiN1)=1;
           end
    else
        p_age_NSW=minpdiff;
    end
    
    % calculate their probs dependent on the matching mechanism
    p_time=p_time_diff(ddays-mindiff+1,2);

    prob=p_time.*psex.*p_age_NSW;
    
    idN(:,1)=T_Nst1.ID;
    idN(:,2)=ddays;
    idN(:,3)=prob;  

    ii1=idN(:,3)>=minp;
    idN=idN(ii1,:);
    
    [~,isort]=sortrows(idN,3,"descend");
    idN=idN(isort,:); % sort in order of decreasing prob
else
    idN=[];

end
