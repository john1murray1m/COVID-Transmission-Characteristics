%% Determine the problem structure for the best match between a set of cases and sequences
function [iGG,allN,A0,A,b,ww,v0]=bipartite_prob(p_match)

% Input:
% p_match is height of the seq being matched: 1) is the G ID, 2) is cell
% for each of (at least) [N.IDs, tdiffs, probs]
% some p_match will have additional columns [pcode, ages, dist to site]

% Output:
%   iGG= the seq being matched (and have some poss cases from p_match
%   allN = the N cases of poss matches over all seqa
%   A = inequality constraint matrix A*v<=b
%   b = rhs
%   ww = the weight coeff
%   v0 = a greedy start value

% the number of variables will be the length of the union of cases that are
% linked

idG=cell2mat(p_match(:,1)); % the seq id
ii2=cellfun(@(x) width(x)~=0,p_match(:,2));

fred=cellfun(@(x) x(:,1),p_match(ii2,2),'UniformOutput',false); % the idN
george=cellfun(@(x) x(:,3),p_match(ii2,2),'UniformOutput',false); % the probs
allN=unique(cell2mat(fred)); % all the N ID that were linked to at 
                            % least one seq
num_cases=length(allN);

% were there any seqs that could not be linked?
% ii0=cellfun(@isempty , fred);

G_unmatched=idG(~ii2);

% the seqs that are left to be matched
iGG=idG(ii2);
num_seqs=length(iGG);


inum=0; % the counter for the variables
iseq_case=cell(num_seqs,1); % the cases for this seq and the p
for i=1:num_seqs
    [~,ia,ib]=intersect(allN,fred{i}); % ia is the case index within allN that 
                                        % are present in this row/seq
    iseq_case(i)={[inum+(1:length(ia))', ia, george{i}(ib), i*ones(size(ia))]}; % the x id, the case it reps, 
                                        % the prob, and the row/seq
    inum=inum+length(ia);
end

A0=cell2mat(iseq_case); 

num_vars=A0(end,1); % the last x index
w=-log(A0(:,3)); % the -log(p) values for these vars and 
maxp=nanmax(w);
wX=2*maxp; %the penalty term for X

% a sequence is matched with at most one real case
As=sparse(num_seqs,num_vars);
bs=ones(num_seqs,1);
for i=1:num_seqs
    idx=iseq_case{i}(:,1); % or idx=A0(:,4)==i;
    As(i,idx)=1;
end

% a real case is assigned to at most one sequence
Ac=sparse(num_cases,num_vars);
bc=ones(num_cases,1);
for i=1:num_cases
    ii1=A0(:,2)==i;
    if any(ii1)
        idx=A0(ii1,1);
        Ac(i,idx)=1;
    end
end

% the starting guess
v0=sparse(num_vars,1);
A0sort=sortrows(A0,3,"descend"); % arrange in decreasing p
iseq=[]; % the seqs assigned so far
icase=[]; % the cases assigned so far
iK=0; % vars assigned so far
i=1;
K=min(num_seqs,num_cases); % this will be the max no. nonzero vars
while iK<K && i<=height(A0sort)
    ivar=A0sort(i,1);
    iseq0=A0sort(i,4);
    icase0=A0sort(i,2);
    if ~ismember(iseq0,iseq) && ~ismember(icase0,icase) 
        v0(ivar)=1;
        iseq=[iseq; iseq0];
        icase=[icase; icase0];
        iK=iK+1;
    end
    i=i+1;
end

% the combined ineq constrains
A=[As; Ac];
b=[bs; bc];

% the penalised weights
ww=w-wX;
