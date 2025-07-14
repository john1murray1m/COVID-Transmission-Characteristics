%% place all state matches into one 2 col matrix
state_names={'NSW', 'VIC', 'QLD', 'SA', 'WA', 'TAS', 'NT', 'ACT'}; % NNDSS

opt_match_all=[];
for i=1:length(state_names)
    st=state_names{i};
    load(['opt_state_soln',st,'_noany']);
    opt_match_all=[opt_match_all; opt_match];
end

save('save_all_matches_noany','opt_match_all') % with no geo labs at

