%% load in the NSW estimated lab sites that are listed in the specified files or are specified directly

% restrict to postcodes to NSW
ii1=T_postcodes.State=='NSW';
T_NSW_postcodes=T_postcodes(ii1,:);

lab_LL=cell(35,3); % 1 is the lab name, 2) is a cell of Lat,Long, 3) is the 
% Modified Monash Model 2019 measure of remoteness and population size 

% get the lab names for NSW
ii1=TG_Delta.Division=='New South Wales';
TG_NSW=TG_Delta(ii1,:);

lab_LL(:,1)=cellstr(unique(TG_NSW.Laboratory));

% 1: 4Cyte
fred=readtable('SHARED\4Cyte Pathology_Pathology_Coordinates_googleAPI_217_NSW_JM.xlsx');
fred=convertvars(fred,"Town",'string');
lab_LL(1,2)={fred{:,{'Latitude','Longitude'}}};
% get remoteness
remote=ones(height(fred),1);
for i=1:height(fred)
    town=fred{i,'Town'};
    ii1=strcmpi(T_NSW_postcodes{:,'Locality'},town);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(1,3)={remote};

% 2: ACT pathology
lab_LL(2,2)={[-35.34653, 249.10257]}; % ACT centre (Garran)
remote=1;
lab_LL(2,3)={remote};

% 3: SAVID, Randwick, use postcodes for South East Sydney LHD
fred=readtable('NSW_postcodes_LHD.xlsx');
fred=convertvars(fred,"LHD",'categorical');
ii1=fred.LHD=='South East Sydney';

pcodes=unique(fred{ii1,'Postcode'});
lat_long=NaN(height(pcodes),2);
for i=1:height(pcodes)
    ii1=find(T_postcodes{:,'Postcode'}==pcodes(i));
    if ~isempty(ii1)
        lat_long(i,:)=T_postcodes{ii1(1),{'Lat','Long'}};
    end
end
lab_LL(3,2)={lat_long};
remote=ones(height(pcodes),1);
for i=1:height(pcodes)
    ii1=T_NSW_postcodes{:,'Postcode'}==pcodes(i);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(3,3)={remote};

% 4: Austech
fred=readtable('SHARED\Austech Medical Laboratories_Coordinates_googleAPI_5_NSW.xlsx');
fred=convertvars(fred,"Town",'string');
lab_LL(4,2)={fred{:,{'Latitude','Longitude'}}};
remote=ones(height(fred),1);
for i=1:height(fred)
    town=fred{i,'Town'};
    ii1=strcmpi(T_NSW_postcodes{:,'Locality'},town);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(4,3)={remote};

% 5: Australian Clinical Labs
fred=readtable('SHARED\Australian_Clinical Labs_Coordinates_googleAPI_NSW_377_withDups_JM.xlsx');
fred=convertvars(fred,"Town",'string');
lab_LL(5,2)={fred{:,{'Latitude','Longitude'}}};
remote=ones(height(fred),1);
for i=1:height(fred)
    town=fred{i,'Town'};
    ii1=strcmpi(T_NSW_postcodes{:,'Locality'},town);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(5,3)={remote};

% 6: Capital Pathology
Capital_LL=[-35.34653, 249.10257; % ACT centre (Garran)
    -36.41849, 148.61602; % Jindabyne
     -34.73977, 149.71905; % Goulburn
     -36.67946, 149.83611; % Bega
    -36.23239, 149.12535]; % Cooma
lab_LL(6,2)={Capital_LL};
remote=[1; 5; 3; 5; 4];
lab_LL(6,3)={remote};


% 7: Douglass Hanly Moir
fred=readtable('SHARED\Douglass_Hanly_Moir_Pathology_Coordinates_googleAPI_208_JM.xlsx');
fred=convertvars(fred,"Town",'string');
lab_LL(7,2)={fred{:,{'Latitude','Longitude'}}};
remote=ones(height(fred),1);
for i=1:height(fred)
    town=fred{i,'Town'};
    ii1=strcmpi(T_NSW_postcodes{:,'Locality'},town);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(7,3)={remote};

% 8: Histopath
% use all of NSW so replace with the nogeo match in the link program

% 9: Laverty
fred=readtable('SHARED\Laverty_Pathology_Coordinates_googleAPI_370_towns_with_dups.xlsx');
fred=convertvars(fred,"Town",'string');
towns=unique(fred.Town);
itowns=zeros(size(towns)); % extract the first entry for each town/suburb
for i=1:length(towns)
    ii1=find(fred.Town==towns(i));
    itowns(i)=ii1(1);
end
fred=fred(itowns,:);
remote=ones(height(fred),1);
for i=1:height(fred)
    town=fred{i,'Town'};
    ii1=strcmpi(T_NSW_postcodes{:,'Locality'},town);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(9,3)={remote};

lab_LL(9,2)={fred{:,{'Latitude','Longitude'}}};

% 10: Liverpool Hospital
lab_LL(10,2)={[-33.92044, 150.93211]};
remote=1;
lab_LL(10,3)={remote};

% 11: Medihealth: use Penrith where they had a Drive thru centre at 243
% High Street
lab_LL(11,2)={[-33.75477, 150.70472]};
remote=1;
lab_LL(11,3)={remote};

% 12: Medlab: now owned by Australian Clinical Labs
lab_LL(12,2)=lab_LL(5,2);
lab_LL(12,3)=lab_LL(5,3);

% 13: MDU-PHL: Victorian path lab so use nogeo

% 14: Nepean Hospital
lab_LL(14,2)={[-33.75929, 150.71385]};
remote=1;
lab_LL(14,3)={remote};

% 15: RPA
lab_LL(15,2)={[-33.88937, 151.1830]}; % SSPWS RPA 
remote=1;
lab_LL(15,3)={remote};

% 16: Other

% 17: Pathwest: WA

% 18:24 Pathology North 
Path_North_LL=[-33.26027, 151.47891; % Wyong Public Hospital
    -33.41940, 151.33945; % Gosford Hospital
    -32.92231, 151.69248; % John Hunter Hospital, Newcastle
    -31.07243, 150.92532; % Tamworth Pathology North
    -28.80822, 153.29143; % Lismore Pathology North
    -30.31612, 153.09335; % Coffs Harbour Pathology North
    -33.82160, 151.19154]; % Royal North Shore Hospital
remote=[1; 1; 1; 3; 3; 3; 1];

lab_LL(18,2)={Path_North_LL(1,:)};
lab_LL(19,2)={Path_North_LL(2,:)};
lab_LL(20,2)={Path_North_LL(3,:)};
lab_LL(21,2)={Path_North_LL(4,:)};
lab_LL(22,2)={Path_North_LL(5,:)};
lab_LL(23,2)={Path_North_LL(6,:)};
lab_LL(24,2)={Path_North_LL(7,:)};

lab_LL(18,3)={remote(1)};
lab_LL(19,3)={remote(2)};
lab_LL(20,3)={remote(3)};
lab_LL(21,3)={remote(4)};
lab_LL(22,3)={remote(5)};
lab_LL(23,3)={remote(6)};
lab_LL(24,3)={remote(7)};


% 25: Pathology West
% read in the LHD suburbs
postcodes_LHD=readtable('NSW_postcodes_LHD.xlsx');
postcodes_LHD=convertvars(postcodes_LHD,{'Suburb','LHD'},'string');
ii1=ismember(postcodes_LHD.LHD,{'Far West NSW','Western NSW'});
postcodes_West=postcodes_LHD(ii1,:);

% read in Pathology NSW sites
pathology_NSW=readtable('SHARED\NSW_Health_Pathology_googleAPI_NSW_157.xlsx');
pathology_NSW=convertvars(pathology_NSW,'Town','string');
iwest=zeros(height(pathology_NSW),1);
for i=1:height(pathology_NSW)
    town=pathology_NSW{i,'Town'};
    iwest(i)=any(arrayfun(@(x) strcmpi(town,x), postcodes_West.Suburb));
end
fred=pathology_NSW(find(iwest),:);
lab_LL(25,2)={fred{:,{'Latitude','Longitude'}}};
remote=ones(height(fred),1);
for i=1:height(fred)
    town=fred{i,'Town'};
    ii1=strcmpi(T_NSW_postcodes{:,'Locality'},town);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(25,3)={remote};


% 26: QML pathology
% Ballina, Byron, Lismore, Murwillumbah, Tweed Heads, Nimbin
QML_LL=[-28.86645, 153.56509; % Ballina
        -28.64553, 153.61659; % Byron Bay
        -28.33020, 153.39583; % Murwillumbah
        -28.18262, 153.53007; % Tweed Heads
        -28.60205, 153.22065]; %Nimbin
lab_LL(26,2)={QML_LL};
remote=[3; 4; 1; 1; 5];
lab_LL(26,3)={remote};

% 27: Safework Laboraties
 % use the location of their lab in Pymble though could collect samples
 % from anywhere as part of Safework NSw but most likely Sydney
 lab_LL(27,2)={[-33.75017, 151.14467]};
remote=1;
lab_LL(27,3)={remote};

 % 28: SEALS
 seals=readtable('SHARED\SEALS_Coordinates_googleAPI_46.xlsx');
 pcodes=cellfun(@(x) str2num(x((end-3):end)),seals.LocationName);
 lab_LL(28,2)={seals{:,{'Latitude','Longitude'}}};

remote=ones(height(pcodes),1);
for i=1:height(pcodes)
    ii1=T_NSW_postcodes{:,'Postcode'}==pcodes(i);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(28,3)={remote};

 % 29: Southern IML Pathology
 % along the Princes highway: Engadine Thirroul, Wollongong, Kiama
 % Berry, Nowra, Vincentia, Ulladulla, Bateman's Bay, Moruya, Narooma
 
 % the distance between Ulladulla and Bateman's is 55kms. So maybe make the
 % prob associated with this dist or half of it, to be the background
 % prob_geo for sites throughout NSW?

IML_LL=[-34.06994, 151.00445; % Engadine
        -34.31718, 150.92199; % Thirroul
        -34.42372, 150.88542; % Wollongong
        -34.67027, 150.84851; % Kiama
        -34.77772, 150.68719; % Berry
        -34.88237, 150.60220; % Nowra
        -35.07052, 150.65227; % Vincentia
        -35.35928, 150.47000; % Ulladulla
        -35.71417, 150.17845; % Batemans Bay
        -35.92044, 150.07985; % Moruya
        -36.24004, 150.11520]; % Narooma
lab_LL(29,2)={IML_LL};
remote=[1; 1; 1; 2; 3; 3; 4; 4; 4; 5; 5];
lab_LL(29,3)={remote};


% 30: St vincent's Sydpath
 sydpath=readtable('SHARED\SYDPATH_googleAPI_Australia_24.xlsx');
 sydpath=convertvars(sydpath,"Town",'string');

 sydpath_LL=sydpath{:,{'Latitude','Longitude'}};
 lab_LL(30,2)={sydpath_LL};
remote=ones(height(sydpath),1);
for i=1:height(sydpath)
    town=sydpath{i,'Town'};
    ii1=strcmpi(T_NSW_postcodes{:,'Locality'},town);
    if any(ii1)
        remote(i)=nanmean(T_NSW_postcodes{ii1,'MMM2019'});
    end
end
lab_LL(30,3)={remote};

 % 31:33
 SSWPS_LL=[-33.83724, 151.09415; % Concord Repat
    -33.92044, 150.93211; % Liverpool Hospital
    -33.88937, 151.1830]; % RPA - these are the sequences that were collected at RPA


 % 31: Concord Hospital
 lab_LL(31,2)={SSWPS_LL(1,:)};
  remote=1;
 lab_LL(31,3)={remote};


 % 32: Liverpool Hospital
 lab_LL(32,2)={SSWPS_LL(2,:)};
   remote=1;
 lab_LL(32,3)={remote};

  % 33: SSWPS RPA (same as #15 but lab record entered differently
 lab_LL(33,2)={SSWPS_LL(3,:)};
  remote=1;
 lab_LL(33,3)={remote};

 % 34: Children's Hospital Westmead. This is calculated separately in
 % Westmead_children_link incorporating p_cum_geo and saved in
 % save_Westmead.mat in p_match_Westmead

 % 35: Virtus Diagnostics -IVF
 % Greenwich, York St CBD, Alexandria, Dee Why, Wahroonga, Kogarah,
 % Westmead, Baulkham Hills, Miranda, Liverpool, Gosford, Auburn St
 % Wollongong, New Lambton Heights
 IVF_LL=[-33.83338, 151.18584; % Greenwich
        -33.86821, 151.20608; % CBD
        -33.90802, 151.19994; % Alexandria
        -33.75072, 151.28605; % Dee Why
        -33.72154, 151.11666; % Wahroonga
        -33.96824, 151.13419; % Kogarah
        -33.80293, 150.98696; % Westmead
        -33.76472, 150.99165; % Baulkham Hills
        -34.03110, 151.10221; % Miranda
        -33.92569, 150.91951; % Liverpool
        -33.42040, 151.33838; % Gosford
        -34.43305, 150.88838; % Wollongong
        -32.92244, 151.69309]; % Nw Lambron Hts (Newcastle)
 lab_LL(35,2)={IVF_LL};
remote=ones(height(IVF_LL),1);
 lab_LL(35,3)={remote};


save('save_NSW_lab_locations',"lab_LL")