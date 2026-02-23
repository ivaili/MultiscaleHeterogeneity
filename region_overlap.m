clear

load('/project/3022037.01/final_analysis/results/figure_helpers/netNames2.mat','netNames');
load('/project/3022037.01/final_analysis/results/figure_helpers/networks2.mat','networks');
addpath(genpath('/home/mrstats/ivaili/Toolboxes/BCT/BCT'));


t = 20;
thr = 2.3;
nP = 10000;
nRoi = 390;

z_asd = dlmread('/home/devpsych/ivaili/3022037.01/Heterogeneity/GPR/MM/z.csv');
z_td = dlmread('/home/devpsych/ivaili/3022037.01/Heterogeneity/GPR/MM/TD/z.csv');

n_asd = size(z_asd,1);
n_td = size(z_td,1);

%% Compute degree of deviation per region
[degree_pos_asd, degree_neg_asd] = compute_deviation_degree(z_asd, thr, nRoi);
[degree_pos_td, degree_neg_td] = compute_deviation_degree(z_td, thr, nRoi);

%% Permutation testing
z_idx = logical([ones(n_asd,1); zeros(n_td,1)]);
K = nRoi;

degree_pos_orig = [degree_pos_asd - 1, degree_pos_td - 1]';
degree_neg_orig = [degree_neg_asd - 1, degree_neg_td - 1]';
Nsubs = size(degree_pos_orig, 1);

cnt_neg = zeros(1,K);
cnt_pos = zeros(1,K);

tic
for p = 1:nP

    if p == 1
        idx = 1:Nsubs;
    else
        idx = randperm(Nsubs);
    end

    asd_mask = z_idx(idx);
    td_mask = ~asd_mask;

    overlap_pos_asd = zeros(t, K);
    overlap_neg_asd = zeros(t, K);
    overlap_pos_td = zeros(t, K);
    overlap_neg_td = zeros(t, K);

    for j = 1:t
        degree_pos = degree_pos_orig >= j;
        degree_neg = degree_neg_orig >= j;

        overlap_pos_asd(j,:) = sum(degree_pos(asd_mask,:)) / n_asd;
        overlap_neg_asd(j,:) = sum(degree_neg(asd_mask,:)) / n_asd;
        overlap_pos_td(j,:) = sum(degree_pos(td_mask,:)) / n_td;
        overlap_neg_td(j,:) = sum(degree_neg(td_mask,:)) / n_td;
    end

    auc_pos_asd = trapz(overlap_pos_asd);
    auc_neg_asd = trapz(overlap_neg_asd);
    auc_pos_td = trapz(overlap_pos_td);
    auc_neg_td = trapz(overlap_neg_td);

    delta_auc_pos = auc_pos_asd - auc_pos_td;
    delta_auc_neg = auc_neg_asd - auc_neg_td;

    if p == 1
        delta_auc_pos1 = delta_auc_pos;
        delta_auc_neg1 = delta_auc_neg;
        perc_diff_pos = (delta_auc_pos ./ auc_pos_td) * 100;
        perc_diff_neg = (delta_auc_neg ./ auc_neg_td) * 100;
        overlap_pos_asd1 = overlap_pos_asd;
        overlap_neg_asd1 = overlap_neg_asd;
        overlap_pos_td1 = overlap_pos_td;
        overlap_neg_td1 = overlap_neg_td;
    else
        cnt_neg = cnt_neg + (delta_auc_neg >= delta_auc_neg1);
        cnt_pos = cnt_pos + (delta_auc_pos >= delta_auc_pos1);
    end

end
toc

punc_neg = cnt_neg / nP;
punc_pos = cnt_pos / nP;
[p_neg] = mafdr(punc_neg, 'BHFDR', true);
[p_pos] = mafdr(punc_pos, 'BHFDR', true);


%% Helper function
function [degree_pos, degree_neg] = compute_deviation_degree(z, thr, nRoi)
    z_log_pos = z > thr;
    z_log_neg = z < -thr;

    adj_pos = recon3dMat(z_log_pos', nRoi);
    adj_neg = recon3dMat(z_log_neg', nRoi);

    nSubs = size(adj_pos, 3);
    degree_pos = zeros(nRoi, nSubs);
    degree_neg = zeros(nRoi, nSubs);

    for i = 1:nSubs
        degree_pos(:,i) = degrees_und(adj_pos(:,:,i));
        degree_neg(:,i) = degrees_und(adj_neg(:,:,i));
    end
end