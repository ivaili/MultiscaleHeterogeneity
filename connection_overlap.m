clear
addpath(genpath('/home/mrstats/ivaili/scripts'))

z_asd = dlmread('/home/devpsych/ivaili/3022037.01/Heterogeneity/GPR/MM/z.csv');
z_td = dlmread('/home/devpsych/ivaili/3022037.01/Heterogeneity/GPR/MM/TD/z.csv');
load('/project/3022037.01/Heterogeneity/edge_level/z_all.mat', 'z_idx')

K = 75855;
nP = 10000;
thr = 2.3;
n_asd = 759;
n_td = 504;

z_idx = logical(z_idx);
z_all = [z_asd; z_td];
z_all_pos = z_all > thr;
z_all_neg = z_all < -thr;
Nsubs = size(z_all_pos, 1);

cnt_neg = zeros(1, K);
cnt_pos = zeros(1, K);

tic
for p = 1:nP

    if p == 1
        idx = 1:Nsubs;
    else
        idx = randperm(Nsubs);
    end

    asd_mask = z_idx(idx);
    td_mask = ~asd_mask;

    overlap_neg = sum(z_all_neg(td_mask,:)) / n_td - sum(z_all_neg(asd_mask,:)) / n_asd;
    overlap_pos = sum(z_all_pos(td_mask,:)) / n_td - sum(z_all_pos(asd_mask,:)) / n_asd;

    if p == 1
        overlap1_neg = overlap_neg;
        overlap1_pos = overlap_pos;
    else
        cnt_neg = cnt_neg + (overlap_neg >= overlap1_neg);
        cnt_pos = cnt_pos + (overlap_pos >= overlap1_pos);
    end

end
toc

punc_neg = cnt_neg / nP;
punc_pos = cnt_pos / nP;
[p_neg] = mafdr(punc_neg, 'BHFDR', true);
[p_pos] = mafdr(punc_pos, 'BHFDR', true);
p_neg_recon = reconMat(p_neg', 390);
p_pos_recon = reconMat(p_pos', 390);