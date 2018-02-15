
name = 'VHbbPlots'
input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_runonskim_20180212/haddjobs/sum_%s.root'
weight = 'weight*puWeight*sign(genWeight)'
enable_reuse_step = True  # try to find output on disk and don't run a step if present


from main_samples import the_samples_dict, sample_colors
the_samples_dict = dict(
    (sample, the_samples_dict[sample])
    for sample in ['TT_powheg', 'ZH125']
)

# from main_selections import the_category_dict
from main_selections import *
the_category_dict = {
    'SR': [cats_SR, sr_sel, main_plotvariables.vars_used_in_selections],
}
