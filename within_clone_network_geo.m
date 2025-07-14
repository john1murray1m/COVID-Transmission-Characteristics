%% generate the network within a single clone
function Etable=within_clone_network_geo(Gopt,seqa,Cclone,iiclone,T_GN,lama,lamd)

Clone_node_name=string(iiclone); % the name of the node for this clone
iiclone0=find(Gopt.Nodes.Name==Clone_node_name);
Date0_within=Gopt.Nodes.Day(iiclone0); % the  date of the within seq root from the clone graph

seq_idx=Cclone{iiclone}; % the seq indices

num_seq=length(seq_idx);
seqs=seqa(seq_idx);
Dates=datenum([seqs.Date]);
[mind,ind_root]=min(Dates); % keep the orig dates for the root node otherwise can generate loops

% if want to modify the dates around their integer values. 
logm=1;
logs=1;
logmu=log(logm^2/sqrt(logs+logm^2));
logsig=sqrt(log(logs/logm^2+1));
Dates00=Dates-lognrnd(logmu,logsig,1,length(Dates)); % log normal dist sd 1 around the original dates

% % remove the logrnd variation in dates if desired. Replace with above if want
% % the logrnd variation
% Dates00=Dates;


dmin=min(Dates00);
dates=Dates00-dmin;

dates_orig=Dates-min(Dates); % this will be used to decide a variety of outgoing nodes


A=CLE_adjacency_within_clone(seqs,dates);

% remove any edges coming into the root
A(:,ind_root)=0;

%% get the nonzero entries and modify these probs by the dist and age diffs
% get the indices within T_GN for these seqs
iign=arrayfun(@(x) find(T_GN.IDG==x), seq_idx);
TT=T_GN(iign,:);
Age=TT.Age;
Long=TT.Long;
Lat=TT.Lat;
[i,j]=find(A);
Age_diff=abs(Age(i)-Age(j));
Dist_diff=deg2km(distance(Lat(i),Long(i),...
                Lat(j),Long(j)));
pp=lama*Age_diff+lamd*Dist_diff;
if any(ismissing(pp)) % some seqs have no case linked to them. For them use either the
                        % original A value or the median
    if any(~ismissing(pp))
        for ij=1:length(i)
            if ismissing(pp(ij)) % otherwise all are missing and leave with A value
                 A(i(ij),j(ij))=A(i(ij),j(ij))+median(pp,"omitmissing");
            else
                 A(i(ij),j(ij))=A(i(ij),j(ij))+pp(ij);
            end
        end
    end
        
else
    for ij=1:length(i)
             A(i(ij),j(ij))=A(i(ij),j(ij))+pp(ij);
    end
end

%% set up the initial digraph
Nodes=string(seq_idx); % label the nodes to be their indices within the whole group
NodeTable0=table(Nodes', 'VariableNames',{'Name'});

G=digraph(A,NodeTable0);

% number the edges
G.Edges.Number=(1:numedges(G))';

% the weights will be the values in A which are the prob
% calculations modified by nt distance
% modify the weights so they are unique
G.Edges.Weight=G.Edges.Weight+(1:numedges(G))'/(numedges(G)*1e4); % so no edges have same weight

%% calculate the optimal network
[Gopt_within,errorflag]=CLE_optimal(G,ind_root);
if errorflag==0 % CLE calculation successful
    Gopt_within.Nodes.Day=dates';
    opt_within={Gopt_within};

                % fix some of the graphs that have a cycle because KickOut did
            % not complete with a forest
            [cyc,cycedge]=allcycles(Gopt_within);
            if ~isempty(cyc)
                edgesout=zeros(1,length(cyc)); % the edges that will be kicked out from the cycles
                for i=1:length(cyc)
                    [~,iworst]=max(Gopt_within.Edges.Weight(cycedge{i}));
                    edgesout(i)=cycedge{i}(iworst);
                end
                Gopt_within=rmedge(Gopt_within,edgesout);
            end

end

%% add in the connecting clone nodes
% get the Clone node of the clone graph coming in to the root
fred=Gopt.Edges.EndNodes;
ii2in=find(cellfun(@(x) strcmp(x,Clone_node_name), fred(:,2)));
in_node=fred(ii2in,1); % the name of the node connecting into this clone at the root
if ~isempty(in_node) % this is not a root node
    in_node_mut=Gopt.Edges.Mut_Dist(ii2in); % the mut dist between these nodes +1000
    in_node_idx=Gopt.Edges.Pre_Idx(ii2in); % the mut dist between these nodes +1000
    in_idx=Cclone{str2num(in_node{1})}(in_node_idx);
    in_clone_date0=seqa(in_idx).Date;
    in_clone_date=datenum(in_clone_date0)-mind;
    in_clone_name=num2str(in_idx);
    % determine the name of the root node within this clone
    ii5=find(indegree(Gopt_within)==0 & outdegree(Gopt_within)>0);
    ii5a=find(indegree(Gopt_within)==0);
    if ~isempty(ii5)
        [~,iia]=min(Gopt_within.Nodes.Day(ii5));
        ii5=ii5(iia); % choose the one with the lowest day value
    elseif ~isempty(ii5a)
        [~,iia]=min(Gopt_within.Nodes.Day(ii5a));
        ii5=ii5a(iia); % choose the one with the lowest day value
    end
    root_within=Gopt_within.Nodes{ii5,1};
    
    % add in the zero (mod(1000) mut dist within the clone)
    Gopt_within.Edges.Mut_Dist=1000*ones(numedges(Gopt_within),1);

    % the edges so far are within the clone. 
    Gopt_within.Edges.Within=ones(numedges(Gopt_within),1);
    
    %% add in an edge coming into this clone
    % get the max of edge numbers
    in_max=max(Gopt_within.Edges.Number);
    if isempty(in_max)
        in_max=0;
    end
    
    NewNode=table({in_clone_name},in_clone_date,'VariableNames', {'Name', 'Day'});
    H=addnode(Gopt_within,NewNode);
    % if height(in_node_mut)~=height(in_max+1)
    %     [in_clone_name root_within]
    %     in_node_mut
    %     in_max+1
    %     root_within
    % end
    NewEdge=table([in_clone_name root_within],in_node_mut,in_max+1,in_node_mut,zeros(size(in_node_mut)),...
        'VariableNames',{'EndNodes','Weight','Number','Mut_Dist','Within'});
    H=addedge(H,NewEdge);
else
    in_max=max(Gopt_within.Edges.Number);
    Gopt_within.Edges.Mut_Dist=1000*ones(numedges(Gopt_within),1);
    Gopt_within.Edges.Within=ones(numedges(Gopt_within),1);
    H=Gopt_within;
end

%% add in the outgoing edges to other clones

% get the out edges going to the roots of other clones
ii2out=find(cellfun(@(x) strcmp(x,Clone_node_name), fred(:,1)));
out_node=fred(ii2out,2); % the name of the node connected to this clone 
% get the root seq name for this clone
out_node_num=str2num(char(out_node));
out_node_seq=out_node_num;
for i=1:length(out_node_num)
   iseq=Cclone{out_node_num(i)};
    [~,idx]=min([seqa(iseq).Date]);
    out_node_seq(i)=iseq(idx);
end
out_node_mut=Gopt.Edges.Mut_Dist(ii2out); % the mut dist between these nodes +1000
out_node_idx=Gopt.Edges.Pre_Idx(ii2out); % the index within this clone of the startfor each edge
out_node_idx_name=Gopt_within.Nodes{out_node_idx,1};
% get the unique idx and distribute the number of times this is
out_node_idx_nameu=out_node_idx_name; % a randomly drawn set of indices with the same day value

%allow nodes with same dates to be possible outgoing nodes rather than a
%single one
% if flag_out_permute % this is set to 1 if don't want all out edges originating from a single node
%                     % for each connection from this clone
    [out_nameu,ia,ic]=unique(out_node_idx_name);
    for k=1:length(out_nameu)
        iia=(ic==k); % the number of times this index is used
        dayu=dates_orig(out_node_idx(ia(k))); % the day value
        iij=dates_orig==dayu; % the nodes with the same day
        iijj=randi(sum(iij),sum(iia),1); % indices drawn randomly amoung 
        % the same day indices for the number of times used in the
        % graph for out edges
        iij0=find(iij);
        out_node_idx_nameu(iia)=Gopt_within.Nodes.Name(iij0(iijj));
    end
% end

out_clones_dates=zeros(length(out_node),1);
 out_clone_name=cell(length(out_node),1);
for j=1:length(out_node)
    iij=find(strcmp(Gopt.Nodes.Name,out_node(j)));
    out_clones_dates(j)=Gopt.Nodes.Day(iij)-Date0_within;
     out_clone_name(j)={num2str(out_node_seq(j))};
end

  NewNodes=table(out_clone_name,out_clones_dates,'VariableNames', {'Name' 'Day'});
 H=addnode(H,NewNodes);
 if ~isempty(in_max) % Gopt_within is not empty
        NewEdges=table([out_node_idx_nameu out_clone_name],out_node_mut,in_max+1+(1:length(ii2out))',...
            out_node_mut,zeros(size(out_node_mut)),...
            'VariableNames',{'EndNodes','Weight','Number','Mut_Dist','Within'});
        H=addedge(H,NewEdges);
 end

%% export the EdgeTable but only with {'EndNodes','Mut_Dist','Within'}
Etable=H.Edges(:,{'EndNodes','Mut_Dist','Within'});
