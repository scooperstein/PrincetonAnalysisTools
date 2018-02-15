
name = 'VHbbPlots'
input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_runonskim_20180121/haddjobs/sum_%s.root'
weight = 'weight'
enable_reuse_step = True  # try to find output on disk and don't run a step if present


from main_samples import the_samples_dict, sample_colors
from main_selections import the_category_dict
