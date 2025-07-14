%% DEtermine the full network of clones and seqs within each clone, link in the cases and run an MC analysis

%% if running without prior values loaded in Matlab, then uncomment these load files
% load('save_all_matches_noany') 
% 
% 
% % save('aligned_input_save','seqa','big_dels','del_grp','Disj_vals','Subs_vals','Bet_vals')
% load('aligned_input_save')
% 
% % save('aligned_clade_del_grps',"ngrps_lin_del",'lineages','del_grp_big')
% load('aligned_clade_del_grps')
% 
% % save('aligned_clone_del_grp_clade_opt',"CLE_adj_clone_grps",'maxnt_dist',...
% %     'ngrps_lin_del_big','min_grp_size')
% load('aligned_clone_del_grp_clade_opt')
% 
% 
% % load in the NSW LHD names - don't need the other things at the moment
% % save('save_LHD_trans_MC','LHD_cell','LHD_names')
% load('save_LHD_trans_MC')

%% specify the weighting for distance between cases lamd, and ages lama
% with weighting exp(-lamd*dist-lama*agediff);
rng default % For reproducibility
lamamax=0.5; 
lamdmax=0.5;
nmc=10; % the number of MC simulations
Lam=lhsdesign(nmc,2);
Lam(:,1)=Lam(:,1)*lamamax;
Lam(:,2)=Lam(:,2)*lamdmax;

%% run the MC calculations
MC_out_tabs=cell(nmc,3); % Tcell  % the first is age_transmissions, 2nd is distances,...
% 3rd is trans between NSW LHD

%% some new results related to numbers of transmissions
T_outd_NSW=table; % the top 10 trans in each sim for NSW
T_outd_VIC=table;
T_outd_AUS=table;

T_outd_edge_NSW=cell(nmc,1); % the edges for the highest trans in each sim
T_outd_edge_VIC=cell(nmc,1);

Out0_perc=NaN(length(age_groups),nmc); % the percent in each age group that have no trans

for imc=1:nmc
    lama=Lam(imc,1);
    lamd=Lam(imc,2);
    %%
    figcount=0;
    Tcell=cell(39,2); % T_GN, T_SET
    
    % % GK is igrp=1
    
    
    %% 
    % igrp=1;
    % jk=10;
    % jjk=1;
    for igrp=1:3
        CLE_adj_grps=CLE_adj_clone_grps{igrp};
        jrow=find(cellfun(@(x) ~isempty(x), CLE_adj_grps(:,1)));
        
        for jk0=1:length(jrow)
            jk=jrow(jk0);
    
            lin_name=CLE_adj_grps{jk,1};
            % choose the particular connected component within this
                Cclone=ngrps_lin_del_big{igrp}{jk,4}; % the vector of del grps
                sc_matz=sc_mat_gen(Cclone); % generate the vector linking seq to clones
    
            for jjk=1:length(CLE_adj_grps{jk,2})
            
                % update the figure counter
                figcount=figcount+1;
    
                Date_val=datestr(CLE_adj_grps{jk,2}(jjk));
                
                
                Gopt=CLE_adj_grps{jk,4}{jjk};
                
                 % fix some of the graphs that have a cycle because KickOut did
                % not complete with a forest
                [cyc,cycedge]=allcycles(Gopt);
                if ~isempty(cyc)
                    edgesout=zeros(1,length(cyc)); % the edges that will be kicked out from the cycles
                    for i=1:length(cyc)
                        [~,iworst]=max(Gopt.Edges.Weight(cycedge{i}));
                        edgesout(i)=cycedge{i}(iworst);
                    end
                    Gopt=rmedge(Gopt,edgesout);
                end
                
                iclone=str2num(char(Gopt.Nodes.Name));
                iseq=[Cclone{iclone}];
                seqs_all=seqa(iseq);
                
                            %% get the estimated case data for these seqs
                ns=length(iseq);
                T_GN=table('Size',[ns,17],'VariableTypes',{'double','double','categorical','categorical','double','categorical',...
                    'categorical','double','datetime','datetime','datetime','double','double','double','categorical,...' ...
                    'datetime','categorical'},'VariableNames',...
                    {'IDG','IDN','Division','State','Age','Age_Grp','Sex','Pcode','Onset_Date','GSpec_Date',...
                    'NSpec_Date','Lat','Long','OutDegree','Died','Date_Died','Place_Acq'});
                T_GN.Age_Grp=categorical(T_GN.Age_Grp,age_groups,'Ordinal',true);
                T_GN.IDG=iseq';
                for j=1:length(iseq)
                    ii1=opt_match_all(:,1)==iseq(j);
                    if any(ii1)
                        T_G=TG_Delta(opt_match_all(ii1,1),:);
                        T_GN{j,'GSpec_Date'}=(T_G.Date);
                        T_GN{j,'Division'}=(T_G.Division);
                    
                        T_N=T_Delta(opt_match_all(ii1,2),:);
                        T_GN{j,'IDN'}=opt_match_all(ii1,2);
                        T_GN{j,'State'}=T_N.STATE;
                        T_GN{j,'Age'}=T_N.AGE;
                        T_GN{j,'Age_Grp'}=T_N.AGE_GRP;
                        T_GN{j,'Sex'}=T_N.SEX;
                        T_GN{j,'Pcode'}=T_N.POSTCODE;
                        T_GN{j,'Onset_Date'}=(T_N.TRUE_ONSET_DATE);
                        T_GN{j,'NSpec_Date'}=(T_N.SPECIMEN_DATE);
                        T_GN{j,'Lat'}=T_N.LAT;
                        T_GN{j,'Long'}=T_N.LONG;
                        T_GN{j,'Died'}=T_N.DIED;
                        T_GN{j,'Date_Died'}=T_N.CV_DATE_DIED_;
                        T_GN{j,'Place_Acq'}=T_N.PLACE_OF_ACQUISITION;
                    else
                        T_GN{j,'IDN'}=NaN;
                        T_GN{j,'Age'}=NaN;
                        T_GN{j,'Pcode'}=NaN;
                        T_GN{j,'Lat'}=NaN;
                        T_GN{j,'Long'}=NaN;
                  
                    end
                end
    
                %% calculate the network of seqs
                Gopt_seqs=clone2seqs_geo(Gopt,seqa,Cclone,iclone,T_GN,lama,lamd);
                outd=outdegree(Gopt_seqs);
                harry=str2double(string(Gopt_seqs.Nodes.Name));
                for i=1:height(T_GN)               
                    ii2=harry==iseq(i);
                    if any(ii2)
                        T_GN{i,'OutDegree'}=outd(ii2);
                    else
                        T_GN{i,'OutDegree'}=NaN;
                    end
                end
                        
    
    
                % add in the clone names of each of the endnodes
                C10=str2num(char(Gopt_seqs.Edges.EndNodes(:,1)));
                C20=str2num(char(Gopt_seqs.Edges.EndNodes(:,2)));
                C1=sc_matz(C10);
                C2=sc_matz(C20);
                Gopt_seqs.Edges.C1=C1;
                Gopt_seqs.Edges.C2=C2;
                
                % remove the edges that have been doubly included between clones
                ii1=find(Gopt_seqs.Edges.Within==0); % the edges between clones, their should be only one between any two
                [~,iid]=unique(Gopt_seqs.Edges{ii1,{'C1','C2'}},'rows');
                Gopt_seqs=rmedge(Gopt_seqs,setdiff(ii1,ii1(iid)));
    
                % figure(figcount+100)
                % clf
                % popt=plot(Gopt_seqs,'NodeLabel',Gopt_seqs.Nodes.Name,'EdgeLabel',mod(Gopt_seqs.Edges.Mut_Dist,1000));
                % title([lin_name,', grp = ',num2str([igrp jk jjk]),', no. clones = ',num2str(height(Gopt.Nodes)),...
                %     ', no. seqs = ',num2str(sum(Gopt.Nodes.Size)),', ',Date_val])
                % if height(Gopt.Nodes)>50
                %              layout(popt,'force','UseGravity',true)
                % else
                %              layout(popt,'layered')
                % end
    
                % now extract the edge data
                ne=numedges(Gopt_seqs);
                seq_ends=NaN(ne,4); % the first is the start node the second the end node
                                                        % with indices being the seq no.
                                                        % The third and fourth are the T_GN
                                                        % indices
                seq_ends(:,1)=str2double(Gopt_seqs.Edges.EndNodes(:,1));
                seq_ends(:,2)=str2double(Gopt_seqs.Edges.EndNodes(:,2));
                
                se=seq_ends(:,[1 2]);
                se2=NaN(size(se(:))); % the translated indices to T_GN
                se1=unique(se(:));
                for j=1:length(se1)
                    ii1=(se(:)==se1(j));
                    ii2=find(T_GN.IDG==se1(j)); % the row in the table T_GN
                    if ~isempty(ii2)
                        se2(ii1)=ii2;
                    end
                end
                seq_ends(:,3)=se2(1:ne);
                seq_ends(:,4)=se2((ne+1):end);
                T_SE.State1=T_GN{seq_ends(:,3),'State'};
                T_SE.State2=T_GN{seq_ends(:,4),'State'};
                T_SE.Odate1=T_GN{seq_ends(:,3),'Onset_Date'};
                T_SE.Odate2=T_GN{seq_ends(:,4),'Onset_Date'};
                T_SE.Sdate1=T_GN{seq_ends(:,3),'NSpec_Date'};
                T_SE.Sdate2=T_GN{seq_ends(:,4),'NSpec_Date'};
                T_SE.Age1=T_GN{seq_ends(:,3),'Age'};
                T_SE.Age2=T_GN{seq_ends(:,4),'Age'};
                T_SE.Age_Grp1=T_GN{seq_ends(:,3),'Age_Grp'};
                T_SE.Age_Grp2=T_GN{seq_ends(:,4),'Age_Grp'};
                T_SE.Sex1=T_GN{seq_ends(:,3),'Sex'};
                T_SE.Sex2=T_GN{seq_ends(:,4),'Sex'};
                T_SE.OutD1=T_GN{seq_ends(:,3),'OutDegree'};
                T_SE.OutD2=T_GN{seq_ends(:,4),'OutDegree'};
                T_SE.Pcode1=T_GN{seq_ends(:,3),'Pcode'};
                T_SE.Pcode2=T_GN{seq_ends(:,4),'Pcode'};
                T_SE.Died1=T_GN{seq_ends(:,3),'Died'};
                T_SE.Died2=T_GN{seq_ends(:,4),'Died'};
                T_SE.Date_Died1=T_GN{seq_ends(:,3),'Date_Died'};
                T_SE.Date_Died2=T_GN{seq_ends(:,4),'Date_Died'};
                T_SE.Place_Acq1=T_GN{seq_ends(:,3),'Place_Acq'};
                T_SE.Place_Acq2=T_GN{seq_ends(:,4),'Place_Acq'};
                T_SE.Dist=deg2km(distance(T_GN{seq_ends(:,3),{'Lat','Long'}},...
                    T_GN{seq_ends(:,4),{'Lat','Long'}}));
                dist_edges=[0 5 10 20 50 100 200 500 1000 1e6];
                dist_groups={'0-<5','5-<10','10-<20','20-<50','50-<100',...
                    '100-<200','200-<500','500-<1000','1000+'};
                T_SE.Dist_Grp=discretize(T_SE.Dist,dist_edges,'categorical',dist_groups);
                
                T_SET=struct2table(T_SE);
                
                Tcell(figcount,1)={T_GN};
                Tcell(figcount,2)={T_SET};
            end
            
        end
    end

    %% only save the tables etc that are to be analysed

        for i=1:height(Tcell)
        if i==1
            T_GN_all=Tcell{i,1};
            T_SET_all=Tcell{i,2};
        else
            T_GN_all=[T_GN_all; Tcell{i,1}];
            T_SET_all=[T_SET_all; Tcell{i,2}];
        end
    end
 
    % change the full names of the states to the shorter versions
     scats=categories(T_GN_all.Division);       
    T_GN_all.Division=renamecats(T_GN_all.Division,state_names([1 3 2 8 7 4 5 6]));
    
    % there were 544 seqs that had no links, so remove these
    ii1=~ismissing(T_GN_all.Division);
    T_GN_all=T_GN_all(ii1,:);
    
    % there are 920 0f 15,702 edges with at least one missing endnode definition
    % remove them
    ii1=~ismissing(T_SET_all.Age_Grp1) & ~ismissing(T_SET_all.Age_Grp2);
    T_SET_all=T_SET_all(ii1,:);



%% tables of transmission aspects
% age groups
    [table_Age_Grp,chi2a,p,labels]=crosstab(T_SET_all.Age_Grp1,T_SET_all.Age_Grp2);
    % trans_age_grps*rec_age_grps

    % Transmissions to each age group from other ages
age_trans=table_Age_Grp./repmat(sum(table_Age_Grp),length(age_groups),1)*100;
% each column considers the break down of trans_agegrps to that rec_age_grp
% column
   
%% repeat the calculations but restricted to NSW where the matches should be better
ii1=T_SET_all.State1=='NSW' & T_SET_all.State2=='NSW';
    [table_Age_Grp_NSW,chi2a,pNSW,labels]=crosstab(T_SET_all{ii1,'Age_Grp1'},T_SET_all{ii1,'Age_Grp2'});
age_trans_NSW=table_Age_Grp_NSW./repmat(sum(table_Age_Grp_NSW),length(age_groups),1)*100;


% 
%% repeat the calculations but restricted to VIC where the matches should be better
ii1=T_SET_all.State1=='VIC' & T_SET_all.State2=='VIC';
    [table_Age_Grp_VIC,chi2a,pVIC,labels]=crosstab(T_SET_all{ii1,'Age_Grp1'},T_SET_all{ii1,'Age_Grp2'});
age_trans_VIC=table_Age_Grp_VIC./repmat(sum(table_Age_Grp_VIC),length(age_groups),1)*100;

MC_out_tabs(imc,1)={{age_trans,age_trans_NSW,age_trans_VIC}};

  % distances
    tab_dist=tabulate(T_SET_all.Dist_Grp);
ii1=T_SET_all.State1=='NSW' & T_SET_all.State2=='NSW';
tab_dist_NSW=tabulate(T_SET_all{ii1,'Dist_Grp'});

ii1=T_SET_all.State1=='VIC' & T_SET_all.State2=='VIC';
tab_dist_VIC=tabulate(T_SET_all{ii1,'Dist_Grp'});


%% how does distance change over time

% NSW
st='NSW';
ii1=T_SET_all.State1==st & T_SET_all.State2==st & ~isnan(T_SET_all.Dist);
T_SET_st=T_SET_all(ii1,:);
mindate=min(T_SET_st.Odate1);
maxdate=max(T_SET_st.Odate2);
dd=days(maxdate-mindate);
dbreak=mindate+floor((0:10)*dd/10);
dmid=mindate+floor(dd/20+(0:9)*dd/10);
dmid_NSW=dbreak;

% break into 10  intervals
Dset1=cell(10,4); % 1) time, 2) Dists, 3) Odate1

for i=1:10
    Dset1(i,1)={dmid(i)};
    ii2=T_SET_st.Odate1>=dbreak(i) & T_SET_st.Odate1<dbreak(i+1);
    Dset1(i,2)={T_SET_st{ii2,'Dist'}};
    Dset1(i,3)={T_SET_st{ii2,'Odate1'}};
    Dset1(i,4)={T_SET_st{ii2,'Age_Grp1'}};
end


% the median distances over all time periods
dist_NSW=NaN(10,1);
iqr_NSW=NaN(10,1);
for i=1:10
    dist_NSW(i)=nanmedian(Dset1{i,2});
    iqr_NSW(i)=iqr(Dset1{i,2});
end



% VIC
st='VIC';
ii1=T_SET_all.State1==st & T_SET_all.State2==st & ~isnan(T_SET_all.Dist);
T_SET_st=T_SET_all(ii1,:);
mindate=min(T_SET_st.Odate1);
maxdate=max(T_SET_st.Odate2);
dd=days(maxdate-mindate);
dbreak=mindate+floor((0:10)*dd/10);
dmid=mindate+floor(dd/20+(0:9)*dd/10);
dmid_VIC=dbreak;
% break into 10  intervals
Dset2=cell(10,4); % 1) time, 2) Dists, 3) Odate1

for i=1:10
    Dset2(i,1)={dmid(i)};
    ii2=T_SET_st.Odate1>=dbreak(i) & T_SET_st.Odate1<dbreak(i+1);
    Dset2(i,2)={T_SET_st{ii2,'Dist'}};
    Dset2(i,3)={T_SET_st{ii2,'Odate1'}};
    Dset2(i,4)={T_SET_st{ii2,'Age_Grp1'}};
end


% the median distances over all time periods
dist_VIC=NaN(10,1);
iqr_VIC=NaN(10,1);
for i=1:10
    dist_VIC(i)=nanmedian(Dset2{i,2});
    iqr_VIC(i)=iqr(Dset2{i,2});
end




%% ACT

st='ACT';
ii1=T_SET_all.State1==st & T_SET_all.State2==st & ~isnan(T_SET_all.Dist);
T_SET_st=T_SET_all(ii1,:);
ngrp=4; %the number of groups in the plot
mindate=min(T_SET_st.Odate1);
maxdate=max(T_SET_st.Odate2);
dd=days(maxdate-mindate);
dbreak=mindate+floor((0:ngrp)*dd/ngrp);
dmid=mindate+floor(dd/(ngrp*2)+(0:ngrp-1)*dd/ngrp);
dmid_ACT=dbreak;
% break into ngrp  intervals
Dset3=cell(ngrp,4); % 1) time, 2) Dists, 3) Odate1, 4) Age_Grp1

for i=1:ngrp
    Dset3(i,1)={dmid(i)};
    ii2=T_SET_st.Odate1>=dbreak(i) & T_SET_st.Odate1<dbreak(i+1);
    Dset3(i,2)={T_SET_st{ii2,'Dist'}};
    Dset3(i,3)={T_SET_st{ii2,'Odate1'}};
    Dset3(i,4)={T_SET_st{ii2,'Age_Grp1'}};
end


% the median distances over all time periods
dist_ACT=NaN(ngrp,1);
iqr_ACT=NaN(ngrp,1);
for i=1:ngrp
    dist_ACT(i)=nanmedian(Dset3{i,2});
    iqr_ACT(i)=iqr(Dset3{i,2});
end



MC_out_tabs(imc,2)={{tab_dist,tab_dist_NSW,tab_dist_VIC, dist_NSW,dmid_NSW,...
    dist_VIC,dmid_VIC,dist_ACT,dmid_ACT}};

%% transmissions between LHD for NSW
trans_LHD=zeros(length(LHD_names));
for i=1:length(LHD_names)
    p1=ismember(categorical(T_SET_all.Pcode1),LHD_cell{i,2});
    for j=1:length(LHD_names)
        p2=ismember(categorical(T_SET_all.Pcode2),LHD_cell{j,2});
        trans_LHD(i,j)=sum(p1 & p2);
    end
end
MC_out_tabs(imc,3)={trans_LHD};


%% number of transmissions

 % get the max outdegree for NSW and VIC separately
    ii1=T_GN_all.State=='NSW' & ~isnan(T_GN_all.OutDegree);
    TT=T_GN_all(ii1,:);
    [~,itt]=sort(TT.OutDegree,'descend');
    outm_NSW=TT(itt(1:10),{'IDN','OutDegree','Age_Grp','Sex','Lat','Long'});
    % get the edges for the highest
    ii1=T_SET_all.State1=='NSW' & T_SET_all.OutD1==TT{itt(1),'OutDegree'} & ...
        T_SET_all.Age_Grp1==TT{itt(1),'Age_Grp'} & T_SET_all.Sex1==TT{itt(1),'Sex'} & ...
        T_SET_all.Pcode1==TT{itt(1),'Pcode'} ;
    T_outd_edge_NSW(imc)={T_SET_all(ii1,:)};
   
 ii1=T_GN_all.State=='VIC' & ~isnan(T_GN_all.OutDegree);
       TT=T_GN_all(ii1,:);
    [~,itt]=sort(TT.OutDegree,'descend');
    outm_VIC=TT(itt(1:10),{'IDN','OutDegree','Age_Grp','Sex','Lat','Long'});
    % get the edges for the highest
    ii1=T_SET_all.State1=='VIC' & T_SET_all.OutD1==TT{itt(1),'OutDegree'} & ...
        T_SET_all.Age_Grp1==TT{itt(1),'Age_Grp'} & T_SET_all.Sex1==TT{itt(1),'Sex'} & ...
        T_SET_all.Pcode1==TT{itt(1),'Pcode'} ;
    T_outd_edge_VIC(imc)={T_SET_all(ii1,:)};

    % no all of Austr
     ii1=~isnan(T_GN_all.OutDegree);
       TT=T_GN_all(ii1,:);
    [~,itt]=sort(TT.OutDegree,'descend');
    outm_AUS=TT(itt(1:10),{'IDN','OutDegree','Age_Grp','Sex','Lat','Long'});

    % % tab the outdegrees - this should be neg bin dist
    % tab_out=tabulate(T_GN_all{:,'OutDegree'});
    T_outd_NSW=[T_outd_NSW; outm_NSW]; % the top 10 trans for each sim in NSW
    T_outd_VIC=[T_outd_VIC; outm_VIC];
    T_outd_AUS=[T_outd_AUS; outm_AUS];

    %% what percentage of each age group have no transmissions
     ii1=~isnan(T_GN_all.OutDegree);
    TT=T_GN_all(ii1,:);
    NT=groupcounts(TT,'Age_Grp');
    ii1=TT.OutDegree==0;
    NT0=groupcounts(TT(ii1,:),'Age_Grp');
    Out0_perc(:,imc)=NT0.GroupCount./NT.GroupCount*100;
end



save('save_align_link_all_noany_geo_shortb_MC','Lam',"MC_out_tabs",...
    'T_outd_NSW','T_outd_VIC','T_outd_AUS',"Out0_perc",'T_outd_edge_VIC',...
    'T_outd_edge_NSW')
