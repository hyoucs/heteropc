% input
dir_str = '~/Dropbox/glasso/concord_results';
incov_str = 'chosen_incov.csv';

% create proximate mask
node_num = 83;
prox_mask = t_edge_adjacency_mask(node_num);

% apply mask on each incov
dir_list = dir(dir_str);
incov_mat = {};
fincov_mat = {};
for i = 3:length(dir_list)
	incov = importdata([dir_str, '/', dir_list(i).name, '/', incov_str]);
	incov_mat{i-2} = incov.data;
	fincov_mat{i-2} = incov.data .* prox_mask;
end
save([dir_str, '/fincov.mat'], 'fincov_mat');

% count preserved nonzero elements
nnz_list = [];
for i = 1:length(incov_mat)
	nnz_list(i,:) = [nnz(incov_mat{i}), nnz(fincov_mat{i})];
end
% 			origin_nnz	filtered_nnz
% EMOTION    	1228421      157533
% GAMBLING    	1124777      156877
% LANGUAGE    	1009365      144487
% MOTOR     	968757       133635
% RELATIONAL    1085051      155887
% SOCIAL    	1055517      156433
% WM    		1052533      147095

% plot masked/filtered incovs
spy(incov_mat{1});
spy(fincov_mat{1});


% compute eigen centrality on edge-to-edge network
fincov_ec = {};
for i = 1:length(fincov_mat)
	fincov_ec{i} = t_eigenvector_centrality_und(fincov_mat{i});
end
task_str = {'EMOTION','GAMBLING','LANGUAGE','MOTOR',...
    'RELATIONAL','SOCIAL','WM'};
for i = 1:length(fincov_mat)-1
    for j = i+1:length(fincov_mat)    
        close all;
        figure;
        semilogy(1:length(fincov_ec{1}),fincov_ec{i},'r+'); hold on;
        semilogy(1:length(fincov_ec{1}),fincov_ec{j},'y+'); 
        print([task_str{i},'-',task_str{j},'.png'], '-dpng', '-r600');
    end
end






