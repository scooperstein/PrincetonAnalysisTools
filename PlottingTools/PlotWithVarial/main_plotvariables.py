

lep_pt_safe = 'lepInd1 > -1 ? selLeptons_pt[lepInd1] : -1.'
CMVAmax = 'Jet_btagCMVA[hJetInd1] > Jet_btagCMVA[hJetInd2] ? Jet_btagCMVA[hJetInd1] : Jet_btagCMVA[hJetInd2]'
CMVAmin = 'Jet_btagCMVA[hJetInd1] > Jet_btagCMVA[hJetInd2] ? Jet_btagCMVA[hJetInd2] : Jet_btagCMVA[hJetInd1]'
lepRelIso1_safe = 'isWmunu || isWenu || isZmm || isZee ? selLeptons_relIso_0 : -1'
lepRelIso2_safe = 'isZmm || isZee ? selLeptons_relIso_1 : -1'


vars_used_in_selections = {
#   histoname            variable                               title;x-axis;y-axis             n-bins   low     up
# signal region
    'Vpt':              ('V_pt',                                ';p_{T}(W)@[GeV];',             30,      0,      500     ),
    'Vmass':            ('V_mass',                              ';M_{ll}@GeV;',                 25,      0,      250     ),
    'lepPt':            (lep_pt_safe,                           ';Lepton@p_{T};',               30,      0,      300     ),
    'jetleadpt':        ('Jet_pt_reg[hJetInd1]',                ';Jet1@P_{T}@GeV;',             30,      0,      500     ),
    'jetsubleadpt':     ('Jet_pt_reg[hJetInd2]',                ';Jet2@P_{T}@GeV;',             25,      0,      250     ),
    'Hpt':              ('H_pt',                                ';p_{T}(jj)@GeV;',              21,      0,      350     ),
    'Hmass':            ('H_mass',                              ';M_{jj}@GeV;',                 25,      0,      250     ),
    'HVdPhi':           ('HVdPhi',                              ';HVdPhi;',                     16,      0,      3.2     ),
    'CMVAmax':          (CMVAmax,                               ';CMVA_{max};',                 60,      -1,     1       ),
    'CMVAmin':          (CMVAmin,                               ';CMVA_{min};',                 60,      -1,     1       ),
    'nAddJets':         ('nAddJets252p9_puid',                  ';nAddJets252p9_puid;',         11,      -.5,    10.5    ),
    'nAddLeptons':      ('nAddLeptons',                         ';nAddLeptons;',                11,      -.5,    10.5    ),
    'met_pt':           ('met_pt',                              ';MET@p_{T}@GeV;',              30,      0,      500     ),
    # MetjDPhi not yet in tree
    'MetTkMetDPhi':     ('MetTkMetDPhi',                        ';MetTkMetDPhi;',               16,      0,      3.2     ),
    'lepMetDPhi':       ('lepMetDPhi',                          ';lepMetDPhi;',                 16,      0,      3.2     ),
    'lepRelIso1':       (lepRelIso1_safe,                       ';1st Lep@Rel@Iso;',            30,      0,      0.3     ),
    'lepRelIso2':       (lepRelIso2_safe,                       ';2st Lep@Rel@Iso;',            30,      0,      0.3     ),
    'BDT':              ('CMS_vhbb_BDT_Wln_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),

# additional from 0-lepton CR selection
    'Vtype':            ('Vtype',                               ';Vtype;',                      10,      -.5,    9.5     ),
    'minMetjDPhi':      ('minMetjDPhi',                         ';MetjDPhi;',                   16,      0,      3.2     ),

# additional from 1-lepton CR selection
    'sigma_met_pt':     ('met_pt / sqrt(htJet30)',              ';#sigma(p_{T}^{miss});',       20,      0,      10      ),

# other
    'cutFlow':          ('cutFlow',                             ';cutFlow;',                    15,      -.5,    14.5    ),
    'controlSample':    ('controlSample',                       ';controlSample;',              30,      -1.5,   28.5    ),
}


vars_1lepCR_only = {  # copied from plotvariables_CS.dat
 # BJet variables
    'leadjetptNR':      ('Jet_pt[hJetInd1]',                    ';Jet1@P_{T}@GeV;',             30,      0,      500     ),
    'leadjetpt':        ('Jet_pt_reg[hJetInd1]',                ';Jet1@P_{T}@GeV;',             30,      0,      500     ),
    'leadjeteta':       ('Jet_eta[hJetInd1]',                   ';Jet1#eta;',                   30,      -3,     3       ),
    'leadjetCMVA':      ('Jet_btagCMVA[hJetInd1]',              ';Jet1@CMVA;',                  30,      -1,     1       ),
    'leadjetCMVAZoom':  ('Jet_btagCMVA[hJetInd1]',              ';Jet1@CMVA;',                  20,      0.94,   1       ),
    'subjetpt':         ('Jet_pt_reg[hJetInd2]',                ';Jet2@P_{T}@GeV;',             25,      0,      250     ),
    'subjetptNR':       ('Jet_pt[hJetInd2]',                    ';Jet2@P_{T}@GeV;',             25,      0,      250     ),
    'subjeteta':        ('Jet_eta[hJetInd2]',                   ';Jet2#eta;',                   30,      -3,     3       ),
    'subjetCMVA':       ('Jet_btagCMVA[hJetInd2]',              ';Jet2@CMVA;',                  30,      -1,     1       ),

# lepton
    'lepRelIso':        ('selLeptons_relIso_0',                 ';Lep@Rel@Iso;',                20,      0,      0.2     ),
    'lepPt':            ('selLeptons_pt[lepInd1]',              ';Lepton@p_{T};',               30,      0,      300     ),
    'lepPhi':           ('selLeptons_phi[lepInd1]',             ';Lepton@#phi;',                30,      -3.2,   3.2     ),
    'lepEta':           ('selLeptons_eta[lepInd1]',             ';Lepton@#eta;',                30,      -3,     3       ),

# DiJet Variables
    'mjj':              ('H_mass',                              ';M_{jj}@GeV;',                 25,      0,      250     ),
    'mjj_noreg':        ('H_mass_noreg',                        ';M_{jj}@NoReg;',               25,      0,      250     ),
    'Hpt':              ('H_pt',                                ';p_{T}(jj)@GeV;',              21,      0,      350     ),
    'dphijj':           ('HJ1_HJ2_dPhi',                        ';dPhi(jj);',                   20,      -3,     3       ),
    'detajj':           ('HJ1_HJ2_dEta',                        ';dEta(jj);',                   20,      0,      3       ),
    'drjj':             ('HJ1_HJ2_dR',                          ';dR(jj);',                     20,      0,      5       ),
    'jjWPtBal':         ('H_pt/V_pt',                           ';jj/W@P_{T}@Bal;',             30,      0,      2       ),
    'HVdPhi':           ('HVdPhi',                              ';HVdPhi;',                     16,      0,      3.2     ),

# MET
    'met_pt':           ('met_pt',                              ';MET@p_{T}@GeV;',              30,      0,      500     ),
    'met_phi':          ('met_phi',                             ';MET@Phi;',                    30,      -3.15,  3.15    ),
    'lepMetDPhi':       ('lepMetDPhi',                          ';Lep+Met@d#phi;',              30,      0,      3.14    ),
    'sqrtHTJet30':      ('sqrt(htJet30)',                       ';sqrt(htJet30);',              25,      0,      50      ),
    'MetSig':           ('met_pt/sqrt(htJet30)',                ';MET@Sig;',                    25,      0,      25      ),
    'min_met':          ('min(met_pt,mhtJet30)',                ';min(MET,MHT);',               25,      0,      400     ),

# BDT
    'BDT':              ('CMS_vhbb_BDT_Wln_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),

# other
    'TopMass':          ('Top1_mass_fromLepton_regPT_w4MET',    ';m_{top}@[GeV];',              20,      0,      600     ),
    'TopMassNoReg':     ('Top1_mass_fromLepton_w4MET',          ';TopMass@No@Reg;',             20,      0,      600     ),
    'nSAJets5_sel':     ('softActivityVH_njets5',               ';num.@SA5@Jets;',              10,      0,      10      ),
    'nAddJets':         ('nAddJets252p9_puid',                  ';nAddJets252p9_puid;',         11,      -.5,    10.5    ),
    'Vmt':              ('V_mt',                                ';W@m_{T};',                    20,      0,      200     ),
    'Vpt':              ('V_pt',                                ';p_{T}(W)@[GeV];',             30,      0,      500     ),
}
vars_1lepCR = vars_used_in_selections.copy()
vars_1lepCR.update(vars_1lepCR_only)

vars_2lepCR_only = {
    # these plot definitions have different binning and overwrite the above for the 2-lepton channel
    # 'BDT':              ('CMS_vhbb_BDT_Wln_13TeV',              ';BDT@Output;',                 20,      -1,     1       ),

    'Vpt':              ('V_pt',        ';p_{T} (V) [GeV];',                40,      0,      400    ),
    'Hpt':              ('H_pt',        ';Regressed p_{T} (jj) [GeV];',     40,      0,      400    ),
    'mjj':              ('H_mass',      ';Regressed m(jj) [GeV];',          17,      0,      255    ),
    'CMVAmax':          (CMVAmax,       ';CMVA_{max};',                     20,      -1,     1      ),
    'CMVAmin':          (CMVAmin,       ';CMVA_{min};',                     20,      -1,     1      ),
    'HVdPhi':           ('HVdPhi',      ';HVdPhi [rad];',                   30,      -3.2,   3.2    ),
    'jjWPtBal':         ('H_pt/V_pt',   ';p_{T} balance after regression;', 25,      0,      2.     ),
    'drjj':             ('HJ1_HJ2_dR',  ';reg. Delta R(jj);',               30,      0,      6      ),

    'lep1Pt':            ('selLeptons_pt[lepInd1]',              ';1st Lepton@p_{T};',               30,      0,      150     ),
    'lep1Phi':           ('selLeptons_phi[lepInd1]',             ';1st Lepton@#phi;',                30,      -3.2,   3.2     ),
    'lep1Eta':           ('selLeptons_eta[lepInd1]',             ';1st Lepton@#eta;',                30,      -3,     3       ),
    'lep1IDMu':          ('selLeptons_looseIdPOG[lepInd1]',      ';1st Lepton mu ID;',               6,      -1.5,    4.5     ),
    'lep1IDEl':          ('selLeptons_eleMVAIdSppring16GenPurp[lepInd1]', ';1st Lepton ele ID;',     6,      -1.5,    4.5     ),

    'lep2Pt':            ('selLeptons_pt[lepInd2]',              ';2nd Lepton@p_{T};',               30,      0,      150     ),
    'lep2Phi':           ('selLeptons_phi[lepInd2]',             ';2nd Lepton@#phi;',                30,      -3.2,   3.2     ),
    'lep2Eta':           ('selLeptons_eta[lepInd2]',             ';2nd Lepton@#eta;',                30,      -3,     3       ),
    'lep2IDMu':          ('selLeptons_looseIdPOG[lepInd2]',      ';2nd Lepton mu ID;',               6,      -1.5,    4.5     ),
    'lep2IDEl':          ('selLeptons_eleMVAIdSppring16GenPurp[lepInd2]', ';2nd Lepton ele ID;',     6,      -1.5,    4.5     ),
}
vars_2lepCR = vars_used_in_selections.copy()
vars_2lepCR.update(vars_2lepCR_only)
