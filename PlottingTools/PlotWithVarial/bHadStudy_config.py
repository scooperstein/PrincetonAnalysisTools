
name = 'VHbbPlotsbHadStudy'
input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_runonskim_20180116/haddjobs/sum_%s.root'
weight = 'weight'
enable_reuse_step = True  # try to find output on disk and don't run a step if present

# fetch the_samples_dict and sample_colors
from bHadStudy_samples import *

# fetch the_category_dict
from bHadStudy_selections import *
