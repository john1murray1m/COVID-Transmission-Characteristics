%% determine sample-time matches for individual GISAID seqs within 2 of the Date 
% and also return the prob dependent on match
% this only considers links to children sampled at Westmead Children's
% hospital
% The ages will be <18 for the children, plus possibly their carers
% also use a prob of distance from Westmead at


function idN=prob_link_child_noany(T_Gst,T_Nst,p_time_diff,p_cum_geo)

Lat_long_westmead=[-33.80139,150.99238];

% now idN has 3 columns: 1) the N ids, 2) the time diffs, 3) the probs of
% the match

% tdiff is the max time allowed between dates
% imatch determines the type of match (within a state): 1, sex and age; 10, sex and age_grp for the unmatched in 1; 
% 2, % age, 3, sex, 4) postcode, 5) age, 6) distance from Westmead




prob_min=1e-8;

if height(T_Gst)>1
    error('There should only be a single G seq')
end

% children are more likely to be tested and sequenced compared to
% adults
prob_child=0.9;

% check if sex info available for this seq
sex=T_Gst.Sex;
age=T_Gst.Age;

minpdiff=min(p_time_diff(:,2));

mindiff=p_time_diff(1,1); % this is -2 at the moment
maxdiff=p_time_diff(end,1); % this is 1 at the moment

dateG=T_Gst.Date;
ddays=days(T_Nst.SPECIMEN_DATE-dateG);
ii1=ddays>=mindiff & ddays<=maxdiff;
ddays=ddays(ii1);

T_Nst1=T_Nst(ii1,:); % within the 2 day limit
% idN=NaN(height(T_Nst1),6); % the matrix of N ids,their date diffs,...
%         probs, postcodes
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
    p_age_Westmead=zeros(length(ages),1);

    if ~ismissing(age)
           iiN1=ages==age; % the same age (the same sex enforced above if known)
           if any(iiN1)
                p_age_Westmead(iiN1)=1;
           end
    else

        ii1=ages<18;
        if any(ii1)
            p_age_Westmead(ii1)=prob_child*minpdiff;
        end
        if any(~ii1)
            p_age_Westmead(~ii1)=(1-prob_child)*minpdiff;
        end    
    end


    pcode=T_Nst1.POSTCODE;
    idN(:,4)=pcode;
    idN(:,5)=ages;

    LL=T_Nst1{:,{'LAT','LONG'}};
    dist_childG=deg2km(distance(LL(:,1),LL(:,2),Lat_long_westmead(1),Lat_long_westmead(2)));

    idN(:,6)=dist_childG;
    % get the prob of this dist or higher
    p_dist=zeros(length(ages),1);
    for i=1:length(ages)
        ii1=find(p_cum_geo(:,1)>dist_childG(i));
        if ~isempty(ii1)
            p_dist(i)=p_cum_geo(ii1(1),2);
        end
    end


    idN(:,1)=T_Nst1.ID;
    idN(:,2)=ddays;

    % calculate their probs dependent on the matching mechanism
    p_time=p_time_diff(ddays-mindiff+1,2);

    idN(:,3)=p_time.*psex.*p_age_Westmead.*p_dist;  
        % drop low or zero probs
    ii1=idN(:,3)>=prob_min;
    idN=idN(ii1,:);

    [~,isort]=sortrows(idN,3,"descend");
    idN=idN(isort,:); % sort in order of prob
    % restrict to the top 100
    if height(idN)>100
        idN=idN(1:100,:);
    end

else
    idN=[];
end
