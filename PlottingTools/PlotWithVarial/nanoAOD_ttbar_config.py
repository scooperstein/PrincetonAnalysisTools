
name = 'VHbbPlotsNanoAOD'
input_pattern = '/eos/cms/store/user/scoopers/VHbbPostNano2017_V1/%s/*/*/*/*.root'
weight = '1'
enable_reuse_step = True  # try to find output on disk and don't run a step if present
treename = 'Events'

the_samples_dict = {
#   name           index  scale             legend      input tokens (goes into input_pattern)
    #'TT_powheg':  [None,  13245./35653.,    'Top pair', ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8']],
    'TT_powheg':  [None,  43714./50265.,    'Top pair', ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8']],
    'Data_Mu':    [None,  1.,               'Data',     ['SingleMuon']],
    # 'Data_El':    [None,  1.,               'Data',     ['DoubleEG']],
}


plot_vars = {
    'Hpt':      ('H_pt',                    ';p_{T}(jj) [GeV];',                    21, 0,      350 ),
    'Vpt':      ('V_pt',                    ';p_{T}(V) [GeV];',                     30, 0,      500 ),
    'mjj':      ('H_mass',                  ';m(jj) [GeV];',                        17, 0,      255 ),
    'CMVAmax':  ('Jet_btagCMVA[hJidx[0]]',  ';CMVA_{max};',                         20, -1,     1   ),
    'CMVAmin':  ('Jet_btagCMVA[hJidx[1]]',  ';CMVA_{min};',                         20, -1,     1   ),
    'jjWPtBal': ('H_pt/V_pt',               ';H/V p_{T} balance;',                  25, 0,      2.  ),
    # 'HVdPhi':   ('HVdPhi',                  ';HVdPhi [rad];',                       30, -3.2,   3.2 ),
    # 'drjj':     ('HJ1_HJ2_dR',              ';reg. Delta R(jj);',                   30, 0,      6   ),
}


categories_CR_2Lepton_TTbar_lowVpt = {
    'CR_2Lepton_lowVpt_TTbar_Mu': ['Vtype == 0', 'V_pt >= 50 && V_pt < 150'],
#    'CR_2Lepton_lowVpt_TTbar_El': '(Vtype == 1 && V_pt >= 50 && V_pt < 150)',
}

sel_CR_2Lepton_TTbar = [
    #'(HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL || HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)',
    #'(HLT_IsoMu27 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)',
    'HLT_IsoMu27',
    '(HLT_IsoMu27 ? Muon_pt[0] > 30 : 1)',
    'nJet >= 2',
    'Jet_Pt[hJidx[0]] > 20',
    'Jet_Pt[hJidx[1]] > 20',
    '(V_mass > 10 && (V_mass < 75 || V_mass > 120))',
    'Jet_btagCMVA[hJidx[0]] > 0.9432',  # CMVAmax > CMVAT
    'Jet_btagCMVA[hJidx[1]] > -0.5884', # CMVAmin > CMVAL
]

the_category_dict = {
#   region block     categories                            sel.                     histograms
    'TTbarCR':      [categories_CR_2Lepton_TTbar_lowVpt,   sel_CR_2Lepton_TTbar,    plot_vars],
    'TTbarCR_N-1':  [categories_CR_2Lepton_TTbar_lowVpt,   sel_CR_2Lepton_TTbar,    plot_vars],
}


sample_colors = {'Top pair': 600}

from varial_ext.treeprojector import TreeProjectorFileBased
TreeProjector = TreeProjectorFileBased
