%% collect the transmissions between and within NSW LHD

% load in the postcodes for the NSW LHD
LHD_postcodes=readtable('NSW_postcodes_LHD.xlsx');
LHD_postcodes(:,{'Var3','Var4'})=[];
LHD_postcodes=convertvars(LHD_postcodes,{'Suburb','LHD'},'categorical');

LHD_names=categorical(unique(LHD_postcodes.LHD));
LHD_cell=cell(length(LHD_names),3); % the name and postcodes for each LHD
for i=1:length(LHD_names)
    ii1=LHD_postcodes.LHD==LHD_names(i);
    LHD_cell(i,1)={LHD_names(i)};
    LHD_cell(i,2)={unique(LHD_postcodes{ii1,"Postcode"})};
    LHD_cell(i,3)={cell2table(tabulate(categorical(LHD_postcodes{ii1,"Postcode"})))};
    % the last one determines the freq of postcodes in an LHD to help
    % eliminate overlap
end

% now remove dupl postcodes
for i=1:(length(LHD_names)-1)
    T1=LHD_cell{i,3};
    for j=i+1:length(LHD_names)
        T2=LHD_cell{j,3};
        [ii1,ii2]=ismember(T1.Var1,T2.Var1);
        if any(ii1)
            ii3=find(ii1);
            ii4=find(ii2);
            idel=ones(length(ii3),1); % which of the tables has that pcode deleted
            for k=1:length(ii3)
                n1=T1{ii3(k),'Var2'};
                n2=T2{ii2(ii4(k)),'Var2'};
                if n1>=n2
                    idel(k)=2;
                end
            end
            ii5=idel==1;
            if any(ii5)
                T1(ii3(ii5),:)=[];
            end
            ii5=idel==2;
            if any(ii5)
                T2(ii2(ii4(ii5)),:)=[];
            end
            LHD_cell(j,3)={T2};
           
        end
    end
    LHD_cell(i,3)={T1};
end
for i=1:length(LHD_names)
    T1=LHD_cell{i,3};
    LHD_cell(i,2)={categorical(string(T1.Var1))};
end
                


save('save_LHD_trans_MC','LHD_cell','LHD_names')