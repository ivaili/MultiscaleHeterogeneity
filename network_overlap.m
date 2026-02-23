clear

load('/project/3022037.01/final_analysis/results/figure_helpers/netNames2.mat','netNames');
load('/project/3022037.01/final_analysis/results/figure_helpers/networks2.mat','networks');
load('/project/3022037.01/Heterogeneity/GPR/MM/ASD/asd_beh.mat')
load('/project/3022037.01/Heterogeneity/GPR/MM/ASD/td_beh.mat')

thr = 2.3;
nP = 10000;
K = 55;
thresholds = 0.01:0.01:0.1;
nThr = length(thresholds);

z_asd = dlmread('/home/devpsych/ivaili/3022037.01/Heterogeneity/GPR/MM/z.csv');
z_td = dlmread('/home/devpsych/ivaili/3022037.01/Heterogeneity/GPR/MM/TD/z.csv');

n_asd = size(z_asd,1);
n_td = size(z_td,1);

ylabs = {"Visual","Somatomotor_Visual","Somatomotor","Dorsalattn_Visual","Dorsalattn_Somatomotor","Dorsalattn","Ventralattn_Visual","Ventralattn_Somatomotor","Ventralattn_Dorsalattn","Ventralattn","Limbic_Visual","Limbic_Somatomotor","Limbic_Dorsalattn","Limbic_Ventralattn","Limbic","Frontoparietal_Visual","Frontoparietal_Somatomotor","Frontoparietal_Dorsalattn","Frontoparietal_Ventralattn","Frontoparietal_Limbic","Frontoparietal","Defaultmode_Visual","Defaultmode_Somatomotor","Defaultmode_Dorsalattn","Defaultmode_Ventralattn","Defaultmode_Limbic","Defaultmode_Frontoparietal","Defaultmode","Thalamus_Visual","Thalamus_Somatomotor","Thalamus_Dorsalattn","Thalamus_Ventralattn","Thalamus_Limbic","Thalamus_Frontoparietal","Thalamus_Defaultmode","Thalamus","Striatum_Visual","Striatum_Somatomotor","Striatum_Dorsalattn","Striatum_Ventralattn","Striatum_Limbic","Striatum_Frontoparietal","Striatum_Defaultmode","Striatum_Thalamus","Striatum","AmygHypp_Visual","AmygHypp_Somatomotor","AmygHypp_Dorsalattn","AmygHypp_Ventralattn","AmygHypp_Limbic","AmygHypp_Frontoparietal","AmygHypp_Defaultmode","AmygHypp_Thalamus","AmygHypp_Striatum","AmygHypp"};

%% Classify edges into networks for ASD and TD
[z_pos_asd, z_neg_asd] = classify_edges(z_asd, thr, networks, netNames);
[z_pos_td, z_neg_td] = classify_edges(z_td, thr, networks, netNames);

%% Permutation testing
z_idx = logical([ones(n_asd,1); zeros(n_td,1)]);
z_all_pos = [z_pos_asd; z_pos_td];
z_all_neg = [z_neg_asd; z_neg_td];
Nsubs = size(z_all_pos,1);

cnt_neg = zeros(1,K);
cnt_pos = zeros(1,K);

for p = 1:nP

    if p == 1
        idx = 1:Nsubs;
    else
        idx = randperm(Nsubs);
    end

    asd_mask = z_idx(idx);
    td_mask = ~asd_mask;

    overlap_pos_asd = zeros(nThr, K);
    overlap_neg_asd = zeros(nThr, K);
    overlap_pos_td = zeros(nThr, K);
    overlap_neg_td = zeros(nThr, K);

    for counter = 1:nThr
        j = thresholds(counter);
        netdev_pos = z_all_pos >= j;
        netdev_neg = z_all_neg >= j;

        overlap_pos_asd(counter,:) = sum(netdev_pos(asd_mask,:)) / n_asd;
        overlap_neg_asd(counter,:) = sum(netdev_neg(asd_mask,:)) / n_asd;
        overlap_pos_td(counter,:) = sum(netdev_pos(td_mask,:)) / n_td;
        overlap_neg_td(counter,:) = sum(netdev_neg(td_mask,:)) / n_td;
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

punc_neg = cnt_neg / nP;
punc_pos = cnt_pos / nP;
[p_neg] = mafdr(punc_neg, 'BHFDR', true);
[p_pos] = mafdr(punc_pos, 'BHFDR', true);


%% Helper function
function [z_pos, z_neg] = classify_edges(z, thr, networks, netNames)
    nRoi = 390;
    z_log_pos = z > thr;
    z_log_neg = z < -thr;

    z_log_pos_sq = recon3dMat(z_log_pos', nRoi);
    z_log_neg_sq = recon3dMat(z_log_neg', nRoi);

    nSubs = size(z_log_pos_sq, 3);
    z_net_pos = zeros(10, 10, nSubs);
    z_net_neg = zeros(10, 10, nSubs);

    for i = 1:nSubs
        [~,~,z_net_pos(:,:,i)] = plotClassifiedEdges2(z_log_pos_sq(:,:,i), networks, 0, netNames);
        [~,~,z_net_neg(:,:,i)] = plotClassifiedEdges2(z_log_neg_sq(:,:,i), networks, 0, netNames);
    end

    z_pos = flatten3dMat(z_net_pos, 0)';
    z_neg = flatten3dMat(z_net_neg, 0)';
end