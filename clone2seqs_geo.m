%% generate the network of seqs from network of clones


function Gopt_seqs=clone2seqs_geo(Gopt,seqa,Cclone,iclone,T_GN,lama,lamd)

% GG=Gopt;

%% GG will be the digraph over the sequences and constructed from an EdgeTable
% Nodes will have [SeqName, CloneName, Date]
% Edges will have [EndNodes (of SeqName), MutDist]
iseq=[Cclone{iclone}];
seqs=seqa(iseq);
% Dates=[seqs.Date]';
Dates=T_GN.Onset_Date;
mindate=min(Dates);


% get the edges for Gopt between Nodes of size 1

% also keep track of the clone edges already included
ii1=Gopt.Nodes.Size==1;

Gopts=Gopt; % now names are seqs
% change names so the seq names do not conflict with the clone names
fred=char(Gopts.Nodes.Name);
fred=strcat(fred,repmat('s',height(fred),1)); % add an s to the end of each name
Gopts.Nodes.Name=cellstr(fred);


Gopts.Edges.Within=zeros(numedges(Gopts),1);
Gopts.Edges.Dist=zeros(numedges(Gopts),1); % the dist between nodes
GGedgetable=[];
if any(ii1)
   Gopts.Nodes.Name(ii1)=cellstr(string([Cclone{iclone(ii1)}]'));

    g1=string(Gopts.Nodes.Name(ii1));
    ii2=cellfun(@(x) any(x==g1),Gopts.Edges.EndNodes(:,1)) & ...
        cellfun(@(x) any(x==g1),Gopts.Edges.EndNodes(:,2));
    GGedgetable=Gopts.Edges(ii2,{'EndNodes','Mut_Dist','Within'});
end

% for each of the clones>1 construct the optimal network
% after modifying the dates
imult=find(~ii1);
if ~isempty(imult) % there are some clones of size > 1
    for i=1:length(imult)
           Etable=within_clone_network_geo(Gopt,seqa,Cclone,iclone(imult(i)),...
               T_GN,lama,lamd);
           % check whether any of these clone connections have already been
           % added
           GGedgetable=[GGedgetable; Etable];
    end
end

Gopt_seqs=digraph(GGedgetable);