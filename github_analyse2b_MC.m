%% analyse the MC results

%% uncomment these save files if needed
% % % save('save_align_link_all_noany_geo_shortb_MC','Lam',"MC_out_tabs",...
% % %     'T_outd_NSW','T_outd_VIC','T_outd_AUS',"Out0_perc",'T_outd_edge_VIC',...
% % %     'T_outd_edge_NSW')
% load('save_align_link_all_noany_geo_shortb_MC')
% 
% % save('save_LHD_trans_MC','LHD_cell','LHD_names')
% load('save_LHD_trans_MC')


% % save('save_background',"state_names_full","state_names","age_groups","gender_names","dist_groups")
% load('save_background')

% create shortened LHD names
LHD_names_short={'Albury','C Coast', 'Far W', 'Hunter N E', 'Ill-Shoal',...
    'Interstate','Mid N Coast','Murrumbidgee','Nepean B M','Northern','N Sydney',...
    'SE Sydney','SW Sydney','Southern','Sydney','Western','W Sydney'};

nmc=length(MC_out_tabs);

% collect all the age_trans 
age_trans_MC=zeros(length(age_groups),length(age_groups),nmc);
age_trans_NSW_MC=zeros(length(age_groups),length(age_groups),nmc);
age_trans_VIC_MC=zeros(length(age_groups),length(age_groups),nmc);

% also the tab_dist values
tab_dist_MC=NaN(length(dist_groups),nmc); % percentages
tab_dist_NSW_MC=NaN(length(dist_groups),nmc); % percentages NSW
tab_dist_VIC_MC=NaN(length(dist_groups),nmc); % percentages VIC

% the changes in median distance over time
% get the dates for changes in dist 
    fred=MC_out_tabs{1,2};
    dmid_NSW_MC=fred{5};
    dmid_VIC_MC=fred{7};
    dmid_ACT_MC=fred{9};

    % drop off the final time
    dmid_NSW_MC=dmid_NSW_MC(1:end-1);
    dmid_VIC_MC=dmid_VIC_MC(1:end-1);
    dmid_ACT_MC=dmid_ACT_MC(1:end-1);

dist_NSW_MC=NaN(length(dmid_NSW_MC),nmc); % the median distances for NSW
dist_VIC_MC=NaN(length(dmid_VIC_MC),nmc); % the median distances for VIC
dist_ACT_MC=NaN(length(dmid_ACT_MC),nmc); % the median distances for ACT

trans_LHD_MC=zeros(length(LHD_names),length(LHD_names),nmc); % the transmissions 
        % between LHD in NSW


for imc=1:nmc
    fred=MC_out_tabs{imc,1};
    age_trans_MC(:,:,imc)=fred{1};
    age_trans_NSW_MC(:,:,imc)=fred{2};
    age_trans_VIC_MC(:,:,imc)=fred{2};

    fred=MC_out_tabs{imc,2};
    tab_dist_MC(:,imc)=cell2mat(fred{1}(:,3));
    tab_dist_NSW_MC(:,imc)=cell2mat(fred{2}(:,3));
    tab_dist_VIC_MC(:,imc)=cell2mat(fred{3}(:,3));

    dist_NSW_MC(:,imc)=fred{4};
    dist_VIC_MC(:,imc)=fred{6};
    dist_ACT_MC(:,imc)=fred{8};
    
    trans_LHD_MC(:,:,imc)=MC_out_tabs{imc,3};
   

end


%% plot the % transmissions to each age group from other ages
xx=(1:length(age_groups));
yym=median(age_trans_MC,3);
yyp=prctile(age_trans_MC,97.5,3);
yyn=prctile(age_trans_MC,2.5,3);

xconf=[xx xx(end:-1:1)];

% restrict the x axis age groups to every second one
ixage=1:2:17;

figure(1)
clf
for ia=1:length(age_groups)
    subplot(4,5,ia)
    yconf=[yyp(:,ia)' yyn(end:-1:1,ia)'];
    p=fill(xconf,yconf,'red');
    p.FaceColor=[1 0.8 0.8];
    p.EdgeColor='none';
    
    hold on
    plot(xx,yym(:,ia)','ro')
    hold off
    set(gca,'XLim',[1,length(age_groups)],'XTick',xx(ixage),...
        'XTickLabel',age_groups(ixage),'YLim',[0,25])
    title(age_groups(ia))
    fontname('Arial')
end

% NSW
yym=median(age_trans_NSW_MC,3);
yyp=prctile(age_trans_NSW_MC,97.5,3);
yyn=prctile(age_trans_NSW_MC,2.5,3);

xconf=[xx xx(end:-1:1)];
figure(2)
clf
for ia=1:length(age_groups)
    subplot(4,5,ia)
    yconf=[yyp(:,ia)' yyn(end:-1:1,ia)'];
    p=fill(xconf,yconf,'blue');
    p.FaceColor=[0.8 0.8 1];
    p.EdgeColor='none';
    
    hold on
    plot(xx,yym(:,ia)','bo')
    hold off
    set(gca,'XLim',[1,length(age_groups)],'XTick',xx(ixage),...
        'XTickLabel',age_groups(ixage),'YLim',[0,25])
    title(age_groups(ia))
        fontname('Arial')

end

% VIC
yym=median(age_trans_VIC_MC,3);
yyp=prctile(age_trans_VIC_MC,97.5,3);
yyn=prctile(age_trans_VIC_MC,2.5,3);

xconf=[xx xx(end:-1:1)];
figure(3)
clf
for ia=1:length(age_groups)
    subplot(4,5,ia)
    yconf=[yyp(:,ia)' yyn(end:-1:1,ia)'];
    p=fill(xconf,yconf,'green');
    p.FaceColor=[0.8 1 0.8];
    p.EdgeColor='none';
    
    hold on
    plot(xx,yym(:,ia)','go')
    hold off
    set(gca,'XLim',[1,length(age_groups)],'XTick',xx(ixage),...
        'XTickLabel',age_groups(ixage),'YLim',[0,25])
    title(age_groups(ia))
        fontname('Arial')

end

%% plot distance distributions
xdist=repmat((1:length(dist_groups))',1,nmc);

figure(4)
clf
    fontname('Arial')

subplot(1,3,1) % all
p=boxchart(xdist(:),tab_dist_MC(:),'Notch','on','MarkerStyle','none');
p.BoxFaceColor=[1 0.8 0.8];
p.BoxEdgeColor='red';
xlabel('Distance','FontWeight','bold')
ylabel('Percentage')
set(gca,'XTick',1:length(dist_groups),'XTickLabel',dist_groups)
title('Australia')
set(gca,'YLim',[0, 35])


subplot(1,3,2) % NSW
p=boxchart(xdist(:),tab_dist_NSW_MC(:),'Notch','on','MarkerStyle','none');
xlabel('Distance','FontWeight','bold')
ylabel('Percentage')
set(gca,'XTick',1:length(dist_groups),'XTickLabel',dist_groups)
title('NSW')
set(gca,'YLim',[0, 35])

subplot(1,3,3) % VIC
p=boxchart(xdist(:),tab_dist_VIC_MC(:),'Notch','on','MarkerStyle','none');
p.BoxFaceColor=[0.8 1 0.8];
p.BoxEdgeColor='green';
xlabel('Distance','FontWeight','bold')
ylabel('Percentage')
set(gca,'XTick',1:length(dist_groups),'XTickLabel',dist_groups)
title('VIC')
set(gca,'YLim',[0, 35])


%% plot the changes in median distances over time
figure(5)
clf
 % NSW
 subplot(1,3,1)

xx=(1:length(dmid_NSW_MC))';
yym=median(dist_NSW_MC,2);
yyp=prctile(dist_NSW_MC,97.5,2);
yyn=prctile(dist_NSW_MC,2.5,2);

xconf=[xx; xx(end:-1:1)];


    yconf=[yyp; yyn(end:-1:1)];
    p=fill(xconf,yconf,'blue');
    p.FaceColor=[0.8 0.8 1];
    p.EdgeColor='none';
    
    hold on
    plot(xx,yym,'bo')
    hold off
    set(gca,'XLim',[1,length(xx)],'XTick',xx,...
        'XTickLabel',datestr(dmid_NSW_MC))
    title('NSW')
    ylabel('Kilometers')

     % VIC
 subplot(1,3,2)

xx=(1:length(dmid_VIC_MC))';
yym=median(dist_VIC_MC,2);
yyp=prctile(dist_VIC_MC,97.5,2);
yyn=prctile(dist_VIC_MC,2.5,2);

xconf=[xx; xx(end:-1:1)];


    yconf=[yyp; yyn(end:-1:1)];
    p=fill(xconf,yconf,'green');
    p.FaceColor=[0.8 1 0.8];
    p.EdgeColor='none';
    
    hold on
    plot(xx,yym,'go')
    hold off
    set(gca,'XLim',[1,length(xx)],'XTick',xx,...
        'XTickLabel',datestr(dmid_VIC_MC))
        set(gca,'YLim',[0 100])

    title('VIC')
    ylabel('Kilometers')

         % ACT
 subplot(1,3,3)

xx=(1:length(dmid_ACT_MC))';
yym=median(dist_ACT_MC,2);
yyp=prctile(dist_ACT_MC,97.5,2);
yyn=prctile(dist_ACT_MC,2.5,2);

xconf=[xx; xx(end:-1:1)];


    yconf=[yyp; yyn(end:-1:1)];
    p=fill(xconf,yconf,'magenta');
    p.FaceColor=[1 0.8 1];
    p.EdgeColor='none';
    
    hold on
    plot(xx,yym,'mo')
    hold off
    set(gca,'XLim',[1,length(xx)],'XTick',xx,...
        'XTickLabel',datestr(dmid_ACT_MC))
    set(gca,'YLim',[0 100])
    title('ACT')
    ylabel('Kilometers')
        fontname('Arial')

        %% plot new graphs related to number of transmissions
       figure(6)
clf
geoscatter(T_outd_NSW.Lat,T_outd_NSW.Long,'ro')
geobasemap streets
 geolimits([-38, -28],[141,154])
title('NSW high transmitters')

figure(7)
clf
geoscatter(T_outd_VIC.Lat,T_outd_VIC.Long,'ro')
geobasemap streets
geolimits([-37, -36],[140.5,150.5])
title('VIC high transmitters')

figure(8)
clf

% plot bar graphs of the highest 10 transmitters from each sim but removing
% duplicates across sim
[~,ii1,~]=unique(T_outd_AUS.IDN);
TT=T_outd_AUS(ii1,:);
P=pivot(TT,Columns='Sex',Rows='Age_Grp');

subplot(1,3,1)
bar(P.Age_Grp,[P.Male, P.Female])
legend('Male','Female')
title('High Transmitters')
xlabel('Age Group')

TT.Age_Grp=categorical(TT.Age_Grp,Ordinal=true);
subplot(1,3,2)
gscatter(TT.Age_Grp,TT.OutDegree,TT.Sex,'rb')
xlabel('Age Group')
ylabel('No. Transmissions')
title('High Transmitters')

subplot(1,3,3)
boxplot(Out0_perc','Labels',age_groups,'Symbol','')
ylabel('Percentage')
title('% of non transmitters in each age')
xlabel('Age Group')
 

    %% calculate median and 95%ci of LHD trans
    % but remove interstate and Albury
    trans_LHD_MC_nonI=trans_LHD_MC;
    trans_LHD_MC_nonI=trans_LHD_MC_nonI(:,[2:5, 7:17],:);
    trans_LHD_MC_nonI=trans_LHD_MC_nonI([2:5, 7:17],:,:);
    LHD_names_short_nonI=LHD_names_short([2:5, 7:17]);

    % turn the number of transmissions to %
    tt=sum(sum(trans_LHD_MC_nonI,1),2);
    trans_LHD_perc_MC=trans_LHD_MC_nonI./repmat(tt,15,15,1)*100;
    trans_LHD_med=median(trans_LHD_perc_MC,3);
    trans_LHD_med_p=prctile(trans_LHD_perc_MC,97.5,3);
    trans_LHD_med_n=prctile(trans_LHD_perc_MC,2.5,3);

    fred=string(round(trans_LHD_med,2));
    fred=fred+' ('+string(round(trans_LHD_med_n,2))+','+string(round(trans_LHD_med_p,2))+')';

    Ttrans=array2table(fred,'VariableNames',string(LHD_names_short_nonI));
    Ttrans=addvars(Ttrans,string(LHD_names_short_nonI)','Before',1,...
        'NewVariableNames','Trans\Rec');

    % reorder from largest self
    [ord,iord]=sort(diag(trans_LHD_med),'descend');
    Rorder=iord;
    Corder=[1 iord'+1];



%% output LHD transmissions to table
      writetable(Ttrans(Rorder,Corder),'Trans_LHD_short_NSW_100.xlsx')
