function sm_ne_NEmember_thresh_set(nefile, thresh)

load(nefile, 'exp_site_nedata')
nedata = exp_site_nedata.nedata;

NEmembers = nedata.NEmember([nedata.NEmember.thresh] == thresh).member;

nedata.MemberThr = thresh;
nedata.NEmembers = NEmembers;
exp_site_nedata.nedata = nedata;
save(nefile, 'exp_site_nedata', '-append')