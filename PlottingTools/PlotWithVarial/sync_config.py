
name = 'VHbbPlotsSync'
input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_sync_tuple_20180208/sum_%s.root'
weight = '1'
enable_reuse_step = True  # try to find output on disk and don't run a step if present


the_samples_dict = {
#   name                     index  scale  legend                  input tokens (sample name in the previous step)
    'DYToLL_HT400to600':    [None,  1.,    'DYToLL_HT400to600',    ['DYToLL_HT400to600']],
    'TT_powheg':            [None,  1.,    'TT_powheg',            ['TT_powheg']],
    'ZH125_powheg':         [None,  1.,    'ZH125_powheg',         ['ZH125_powheg']],
}


plot_vars = {
    'count':        ('1',           ';;event count',    1,  .5, 1.5),
    'puWeight':     ('puWeight',    ';puWeight;',       60, 0,  1.5),
    'Hpt':          ('H_pt',        ';p_{T}(jj)@GeV;',  21, 0,  350),
    'Vpt':          ('V_pt',        ';p_{T}(W)@[GeV];', 30, 0,  500),
}
import main_plotvariables as pv
plot_vars.update(pv.vars_used_in_selections)

plot_vars_2lep = plot_vars.copy()
plot_vars_2lep.update(pv.vars_2lepCR_only)


import main_selections as ms
sel = ms.no_sel

the_category_dict = {
#   region block     categories                     sel.    histograms
    'CR_0Lepton':   [ms.cats_CR_0Lepton,            sel,    plot_vars],
    'CR_2LeptLO':   [ms.cats_CR_2Lepton_lowVpt,     sel,    plot_vars_2lep],
    'CR_2LeptHI':   [ms.cats_CR_2Lepton_highVpt,    sel,    plot_vars_2lep],
}


from main_config import sample_colors
