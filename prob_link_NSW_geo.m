%% determine sample-time matches for individual GISAID seqs 
% when there are multiple sites to which the individual is compared
% the nogeo match has already been processed


function idN=prob_link_NSW_geo(T_Gst,T_Nst,idN,p_cum_geo,GLL,remote)

% GLL is an array of Lats, Longs for all sites


% now idN has 3 columns: 1) the N ids, 2) the time diffs, 3) the probs of
% the match; 4) distance from pathology site

prob_min=1e-8;



if height(T_Gst)>1
    error('There should only be a single G seq')
end

    LL=T_Nst{:,{'LAT','LONG'}};
    dist0=NaN(height(LL),height(GLL));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    scale_rural0=1+(remote'-1)/3; % this scales remotefrom a number between 1 (city) 
    % to 7 (very remote and small popn) to an number between 1 and 3. This
    % will give a max distance of 100kms for city links to 300 for the most
    % remote.
    scale_rural=repmat(scale_rural0,height(LL),1);

    for k=1:height(GLL)
        % will need to scale the distances by factor_rural
        dist0(:,k)=deg2km(distance(LL(:,1),LL(:,2),GLL(k,1),GLL(k,2)));
    end
    [dist_Gs,ik]=nanmin(dist0./scale_rural,[],2); % min for each N over all path sites
                        % scaling by rural or not
     dist_G=dist_Gs.*scale_rural(1,ik)'; % the real distance

    idN(:,4)=dist_G;
    % get the prob of this dist or higher
    p_dist=zeros(height(idN),1);
    
    for i=1:height(idN)
        ii1=find(p_cum_geo(:,1)>dist_Gs(i)); % use the rural scaled distance for the prob
        if ~isempty(ii1)
            p_dist(i)=p_cum_geo(ii1(1),2);
        end
    end

    idN(:,3)=idN(:,3).*p_dist;  

    % drop low or zero probs
    ii1=idN(:,3)>=prob_min;
    idN=idN(ii1,:);

    [~,isort]=sortrows(idN,3,"descend");
    idN=idN(isort,:); % sort in order of decreasing prob

    % restrict to top 100
    if height(idN)>100
        idN=idN(1:100,:);
    end
