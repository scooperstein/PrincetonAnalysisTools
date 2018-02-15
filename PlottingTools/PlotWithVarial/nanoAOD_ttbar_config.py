
name = 'VHbbPlotsSync'
input_pattern = '/eos/cms/store/user/scoopers/VHbbPostNano2017_V1/%s/*/*/*/*.root'
weight = '1'
enable_reuse_step = True  # try to find output on disk and don't run a step if present
treename = 'Events'

the_samples_dict = {
#   name                     index  scale  legend      input tokens (sample name in the previous step)
    'TT_powheg':            [None,  1.,    'Top pair', ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8']],
}


plot_vars = {
    'count':        ('1',           ';;event count',    1,  .5, 1.5),
    'puWeight':     ('puWeight',    ';puWeight;',       60, 0,  1.5),
    'Hpt':          ('H_pt',        ';p_{T}(jj)@GeV;',  21, 0,  350),
    'Vpt':          ('V_pt',        ';p_{T}(W)@[GeV];', 30, 0,  500),
}


categories_CR_2Lepton_TTbar_lowVpt = {
    'CR_2Lepton_lowVpt_TTbar_Mu': '(Vtype == 0 && V_pt >= 50 && V_pt < 150)',
    'CR_2Lepton_lowVpt_TTbar_El': '(Vtype == 1 && V_pt >= 50 && V_pt < 150)',
}

sel_CR_2Lepton_TTbar = [
    '(HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)',
    'nJet >= 2',
    'jet_Pt[0] > 20',
    'jet_Pt[1] > 20',
    '(V_mass > 10 && (V_mass < 75 || V_mass > 120))',
    '(Jet_btagCMVA[0] > Jet_btagCMVA[1] ? Jet_btagCMVA[0] : Jet_btagCMVA[1]) > 0.9432',  # CMVAmax > CMVAT
    '(Jet_btagCMVA[0] > Jet_btagCMVA[1] ? Jet_btagCMVA[1] : Jet_btagCMVA[0]) > -0.5884', # CMVAmin > CMVAL
]

the_category_dict = {
#   region block     categories                             sel.                     histograms
    'TTbar':         [categories_CR_2Lepton_TTbar_lowVpt,   sel_CR_2Lepton_TTbar,    plot_vars],
}


sample_colors = {'Top pair': 600}
