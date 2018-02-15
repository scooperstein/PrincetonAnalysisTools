from main_selections import *

cats = {
    'CR_0Lepton_VpHF':     '(controlSample==3 && isZnn==1)',
    'CR_1Lepton_VpHF_Mu':  '(controlSample==13 && isWmunu==1 && Vtype==2)',
    'CR_2Lepton_lowVpt_VpHF_Mu':  '(controlSample==23 && isZmm==1%s)' % low_vpt,
    'CR_2Lepton_highVpt_VpHF_Mu':  '(controlSample==23 && isZmm==1%s)' % high_vpt,
}

plot_vars = {
    'GenStatus2bHad0_sorted_pt':    ('GenStatus2bHad_sorted_pt[0]',  ';GenStatus2bHad_sorted_pt[0];', 25, 0, 250),
    'GenStatus2bHad1_sorted_pt':    ('GenStatus2bHad_sorted_pt[1]',  ';GenStatus2bHad_sorted_pt[1];', 25, 0, 250),
    'GenStatus2bHad2_sorted_pt':    ('GenStatus2bHad_sorted_pt[2]',  ';GenStatus2bHad_sorted_pt[2];', 25, 0, 250),
    'GenStatus2bHad_pt':            ('GenStatus2bHad_pt',            ';GenStatus2bHad_pt;',           25, 0, 250),
    'GenStatus2bHad_eta':           ('GenStatus2bHad_eta',           ';GenStatus2bHad_eta;',          20, -5, 5),
}

the_category_dict = {
    'Testing':   [cats, [], plot_vars],
}
