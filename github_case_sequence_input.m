%% This file reads in the Australian Delta case data and also the sequence data. It sets it up for linkages

%% Input files

% read in Delta data 
T_Delta=readtable('J:Attachment B 550_2022 Data.XLSX','Sheet','Delta');

% read in the postcode data obtained from https://www.matthewproctor.com/australian_postcodes
T_postcodes=readtable("australian_postcodes.xlsx");

% read in the Delta sequences from previous network calculations
% Aligned Delta sequences
% save('aligned_input_save','seqa','big_dels','del_grp','Disj_vals','Subs_vals','Bet_vals')
load('../aligned_input_save')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify behaviour groups
age_groups={'0-4' '5-9' '10-14' '15-19' '20-24' '25-29' '30-34' '35-39' ...
            '40-44' '45-49' '50-54' '55-59' '60-64' '65-69' '70-74' ...
            '75-79' '80-84' '85+'};

% the full state names in NNDSS and GISAID in order
state_names={'NSW', 'VIC', 'QLD', 'SA', 'WA', 'TAS', 'NT', 'ACT'}; % NNDSS
state_names_full={'New South Wales', 'Victoria', 'Queensland', 'South Australia', ...
    'Western Australia', 'Tasmania', 'Northern Territory', 'Australian Capital Territory'}; % GISAID

gender_code={'F'; 'M'; 'X'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reformat the data

%% postcodes
T_postcodes=T_postcodes(:,{'Postcode','Locality','State','Long','Lat',...
    'SA3Name','SA4Name','MMM2019'});

%% cases
T_postcodes=convertvars(T_postcodes,1:width(T_postcodes),'string');
T_postcodes=convertvars(T_postcodes,{'Postcode','Long','Lat','MMM2019'},'double');

%convert some to categorical
varid=[1 3:7 9 15 16 18];
T_Delta=convertvars(T_Delta,varid,'categorical');

% %% Exclude some that are subtyped as omicron or are untyped (or Alpha)
type_drop={'B.1.1.7' 'B.1.1.529' 'BA.1' 'BA.2' 'BA.4' 'BA.4.6' 'BA.5' };
%     'No Data', 'UNTYPED'};
ii1=ismember(T_Delta.MAPPED_SEROGROUP_SUBTYPE,type_drop);
T_Delta(ii1,:)=[];

% the sequences were sampled up to 24 October 2021 so restrict to this date
% and where a SPECIMEN_DATE is recorded
ii1=(T_Delta.TRUE_ONSET_DATE<='24-Oct-2021' & ~ismissing(T_Delta.SPECIMEN_DATE));
T_Delta=T_Delta(ii1,:);

% EXCLUDE some variables
var_drop={'COVIDWAVE' 'ORGANISM_NAME' 'STATE_SEROGROUP_SUBTYPE'};
T_Delta(:,var_drop)=[];

% shorten some of the names
T_Delta.Properties.VariableNames("MAPPED_SEROGROUP_SUBTYPE") = "SUBTYPE";
T_Delta.Properties.VariableNames("LAB_DIAGNOSIS_METHOD") = "TEST";
T_Delta.Properties.VariableNames("RESIDENT_POSTCODE") = "POSTCODE";
T_Delta.Properties.VariableNames("RESIDENT_LOCATION") = "LOCATION";
T_Delta.Properties.VariableNames("AGE_AT_ONSET") = "AGE";

% throw in a sequential ID number
ID=(1:height(T_Delta))';
T_Delta=addvars(T_Delta,ID,'Before',1);

age_edges=[0:5:85 110];
AGE_GRP=discretize(T_Delta.AGE,age_edges,'categorical',age_groups);
T_Delta=addvars(T_Delta,AGE_GRP,'after','AGE');

% attach postcode Lat and Long
long_lat=NaN(height(T_Delta),2);
for i=1:height(T_Delta)
    ii1=find(T_postcodes.Postcode==T_Delta.POSTCODE(i));
    if ~isempty(ii1)
        long_lat(i,:)=T_postcodes{ii1(1),{'Long','Lat'}};
    end
end
T_Delta.LONG=long_lat(:,1);
T_Delta.LAT=long_lat(:,2);

%% sequence table
stf=categorical(state_names_full);

TG_Delta=struct2table(seqa);
TG_Delta.Sequence=[];
varnames={'Sex','Lineage','Clade','Division','Location','Laboratory',...
    'LaboratorySub'};
TG_Delta = convertvars(TG_Delta,varnames,'categorical');
TG_Delta=convertvars(TG_Delta,{'Header','ISL'},'string');

% convert strings in cells to nums
TG_Delta.Age=standardizeMissing(TG_Delta.Age,{'-6','NaN','unknown'});
TG_Delta.Age=str2double(string(TG_Delta.Age));

% replace U in Sex with unknown
ii1=TG_Delta.Sex=='U';
TG_Delta{ii1,'Sex'}=categorical("unknown");


% replace Other in Sex with X to be consistent with the NNDSS
ii1=TG_Delta.Sex=='Other';
TG_Delta{ii1,'Sex'}=categorical("X");
% % then the sequences are
% seqs=seqa(fred(:,1));

% replace 'Victoria?' with 'Victoria'
ii1=TG_Delta.Division=='Victoria?';
TG_Delta{ii1,'Division'}=categorical("Victoria");

% Some of the labs have variations, so fix to the most common
ii1=contains(string(TG_Delta.LaboratorySub),'MDU-PHL') | contains(string(TG_Delta.LaboratorySub),'Microbiological Diagnostics Unit');
TG_Delta{ii1,'LaboratorySub'}=categorical("Microbiological Diagnostic Unit - Public Health Laboratory (MDU-PHL)");

ii1=contains(string(TG_Delta.LaboratorySub),'Westmead');
TG_Delta{ii1,'LaboratorySub'}=categorical("NSW Health Pathology - Institute of Clinical Pathology and Medical Research; Westmead Hospital; University of Sydney");

ii1=contains(string(TG_Delta.LaboratorySub),'PHV-FSS') | ...
    contains(string(TG_Delta.LaboratorySub),'Forensic');
TG_Delta{ii1,'LaboratorySub'}=categorical("Queensland Health Forensic and Scientific Services");

ii1=contains(string(TG_Delta.LaboratorySub),'PathWest');
TG_Delta{ii1,'LaboratorySub'}=categorical("PathWest Laboratory Medicine WA Microbial Surveillance Unit");

% now fix some of the lab names for those that submit for sequencing
ii1=contains(string(TG_Delta.Laboratory),'4Cyte');
TG_Delta{ii1,'Laboratory'}=categorical("4Cyte Pathology");

ii1=contains(string(TG_Delta.Laboratory),'Australian Clinical Labs');
TG_Delta{ii1,'Laboratory'}=categorical("Australian Clinical Labs (formerly Healthscope Pathology)");

ii1=contains(string(TG_Delta.Laboratory),'Westmead');
TG_Delta{ii1,'Laboratory'}=categorical("The Children's Hospital at Westmead");

ii1=contains(string(TG_Delta.Laboratory),'HISTOPATH');
TG_Delta{ii1,'Laboratory'}=categorical("Histopath");

ii1=contains(string(TG_Delta.Laboratory),'Laverty');
TG_Delta{ii1,'Laboratory'}=categorical("Laverty Pathology");

ii1=contains(string(TG_Delta.Laboratory),'MedLab');
TG_Delta{ii1,'Laboratory'}=categorical("Medlab Pathology");

ii1=contains(string(TG_Delta.Laboratory),'Microbiological');
TG_Delta{ii1,'Laboratory'}=categorical("Microbiological Diagnostic Unit - Public Health Laboratory (MDU-PHL)");

ii1=contains(string(TG_Delta.Laboratory),'PathWest Laboratory Medicine');
TG_Delta{ii1,'Laboratory'}=categorical("PathWest Laboratory Medicine WA");

ii1=contains(string(TG_Delta.Laboratory),'Southern');
TG_Delta{ii1,'Laboratory'}=categorical("Southern. IML Pathology");

ii1=contains(string(TG_Delta.Laboratory),'St Vincent');
TG_Delta{ii1,'Laboratory'}=categorical("St Vincents Pathology (SydPath)");

% change dates to same format as NNDSS
TG_Delta.Date.Format='dd-MMM-yyyy';
TG_Delta.DateSub.Format='dd-MMM-yyyy';

% throw in a sequential ID number
IDG=(1:height(TG_Delta))';
TG_Delta=addvars(TG_Delta,IDG,'Before',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some probabilities
% the time difference probabilities
p_time_diff = [-2.0000    0.0075; -1.0000    0.0225; 0    0.9595; 1.0000    0.0105];

% prob for age groups
tabg=tabulate(T_Delta.AGE_GRP);
p_age_grp=[(1:length(age_groups))' cell2mat(tabg(:,3))/100];

% % prob for ages - not used
% restrict ages from 0 to 100 and make uniform
p_age=[(0:100)' 1/101*ones(101,1)];

tab_gender=tabulate(T_Delta{~ii3,'SEX'});
p_gender=[(1:3)' cell2mat(tab_gender([1,2,4],3))/100];

%% use the unique age-gender matches for NSW SAVID to determine distance
% probs
ii3=TG_Delta.Sex~='unknown' & ~ismissing(TG_Delta.Age);

TG_sa=TG_Delta(ii3,:);

Link_sa=cell(2,3);

for istate=1:2
        stG=state_names_full(istate);
    stN=state_names(istate);
    ii1=TG_sa.Division==stG;
    TG_state=TG_sa(ii1,:); % the GISAID sequences for this state
    ii2=T_Delta.STATE==stN;
    TN_state=T_Delta(ii2,:); % the NNDSS cases  for this state

    [P_Link,idG,idN]=prob_link_sa(TG_state,TN_state);

    Link_sa(istate,1)={P_Link}; % this is the smallest time diff in sample times
    Link_sa(istate,2)={idG}; % this is the ID in the GISAID date
    Link_sa(istate,3)={idN}; % this is the cell of ids and time diffs


end

% all together
idN1=Link_sa{1,3};
idN2=Link_sa{2,3};
idNN=[cell2mat((idN1(:,2))); cell2mat((idN2(:,2)))];
idGG=cell2mat(Link_sa(:,2)); % all G matches
idNN0=[Link_sa{1,3}; Link_sa{2,3}];
% the unique matches within 2 of the G Date
idNN1=find(cell2mat(cellfun(@(x) length(x)==1,idNN0(:,2),'UniformOutput',false)));
fred=cell2mat(idNN0(idNN1,:));


match_age_sex=[idGG(idNN1) fred]; % the links [G.ID N.Id N.time_diff]
% any duplication for the N links?
tab2=tabulate(match_age_sex(:,2));
idup=tab2(:,2)>1;
iidup=tab2(idup,1); % the N.IDs of the duplicates. See if any are closer in time
dup_vals=cell(length(iidup),2); % the iidup value from N, the [G vals, the time diffs]
for i=1:length(iidup)
    fred=iidup(i);
    ii1=match_age_sex(:,2)==fred;
    dup_vals(i,1)={fred};
    dup_vals(i,2)={match_age_sex(ii1,[1 3])};
end

iisanon=arrayfun(@(x) ismember(x,iidup),match_age_sex(:,2));
umas=match_age_sex(~iisanon,:);

lab_names=unique(TG_Delta.Laboratory);

ilab=find(contains(string(lab_names),'Randwick'));
ii2=find(T_postcodes.Locality=="RANDWICK");

TT=TG_Delta(umas(:,1),:);

lab_name=lab_names(ilab);
ii1=TT.Laboratory==lab_name;
TT=TT(ii1,:);
um=umas(ii1,:);

% also restrict to time_diff range
ii1=um(:,3)>=p_time_diff(1,1) & um(:,3)<=p_time_diff(end,1);
um=um(ii1,:);
TT=TT(ii1,:);

% also don't want to include infections from overseas
TN=T_Delta(um(:,2),:);
ii1=TN.PLACE_OF_ACQUISITION=="Australia";
um=um(ii1,:);
TT=TT(ii1,:);

pcode=T_postcodes.Postcode(ii2(1));
Long=T_postcodes.Long(ii2(1));
Lat=T_postcodes.Lat(ii2(1));
G_LL=[Long, Lat];

exact_dist_time=NaN(height(TT),3); % the dist to Randwick and time diff
% overseas (NaN or 2999) will have NaN distance
for i=1:height(TT)

    % the distance in kms
    fred=T_Delta{um(i,2),{'LONG','LAT','POSTCODE'}};
    dist_2_G=deg2km(distance(fred(:,2),fred(:,1),G_LL(2),G_LL(1)));
    exact_dist_time(i,:)=[dist_2_G,um(i,3) fred(:,3)];
end

fred=exact_dist_time(:,1);
ii1=~isnan(fred);
fred=fred(ii1);
[f,x]=ecdf(fred);
p_cum_geo=[x,1-f];