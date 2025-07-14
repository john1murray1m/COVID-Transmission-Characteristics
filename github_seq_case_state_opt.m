%% Determine the optimal match between cases and sequences for each state


%   iGG= the seq being matched (and have some poss cases from p_match
%   allN = the N cases of poss matches over all seqa
%   A0 = description of the variables in terms of seq/cases
%   A = inequality constraint matrix A*v<=b
%   b = rhs

%   ww = the weight coeff
%   v0 = a greedy start value

%% run the optimization routine for each value of st
% st='QLD';
st='VIC';
% st='NSW';
% st='SA';
% st='ACT';
% st='WA';
% st='NT';
% st='TAS';

    load(['save_links_',st,'_noany']) % if ignoring age and gender prob (use

p_match=eval(['p_match_',st]);
ii1=cellfun(@(x) ~isempty(x),p_match(:,2));
p_match=p_match(ii1,:);

% restrict to the top (plim) possible case matches
    plim=100;

for i=1:height(p_match)
    fred=p_match{i,2};
    if height(fred)>plim
        fred=fred(1:plim,:);
        p_match(i,2)={fred};
    end
end

[iGG,allN,A0,A,b,ww,v0]=bipartite_prob(p_match);

% run the opt routine
intcon=1:length(v0); % all vars are integer
lb=zeros(size(v0));
ub=ones(size(v0));

% try LP without integer constraints
options = optimoptions('linprog','Display','iter');
[vv,fval,exitflag,output]=linprog(ww,A,b,[],[],lb,ub,options)

% get the seq, case matches
idxv=vv>1e-5;

% check the integrality of soln
disp(['Min of vv>1e-5 is ', num2str(min(vv(idxv)))]);

disp(['nnz(vv) is ',num2str(nnz(vv))]);
disp(['no. of vv>1e-5 is ',num2str(sum(idxv))]);
disp(['Max of vv<=1e-5 is ',num2str(max(vv(~idxv)))]);

A00=A0(idxv,:);

opt_match=[iGG(A00(:,4)), allN(A00(:,2))];

save(['opt_state_soln',st,'_noany'],'opt_match','vv')
