
name = 'VHbbPlotsSyncEventsOnlyInOneFwk'
input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_sync_tuple_20180208/sum_%s.root'
weight = '1'
enable_reuse_step = True  # try to find output on disk and don't run a step if present
n_parallel_treeplotters = 7
from main_config import sample_colors
import ROOT


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
    'cutFlow':      ('cutFlow',     ';cutFlow;',        15,-.5,14.5),
}
import main_plotvariables as pv
plot_vars.update(pv.vars_used_in_selections)
plot_vars_2lep = plot_vars.copy()
plot_vars_2lep.update(pv.vars_2lepCR_only)

import main_selections as ms
cat_names = ms.cats_CR_2Lepton_lowVpt.keys() + ms.cats_CR_2Lepton_highVpt.keys()


def tree_prep(tree, lists, fname, region):
    import os
    all_region_lists = lists.get(os.path.basename(fname), {})
    mask_json = set(all_region_lists.get(region, []))

    if not mask_json:
        tree.SetEventList(ROOT.TEventList())
    else:
        mask_root = ROOT.TEventList()
        set(mask_root.Enter(e.GetReadEntry()) for e in tree if e.evt in mask_json)
        tree.SetEventList(mask_root)


def make_teventlists(evt_json):
    def get_teventlist(fname, single_file_json):
        print '='*10, 'starting', fname
        region_tuple = tuple(
            (rn, set(events), ROOT.TEventList())
            for rn, events in single_file_json.iteritems()
            if events
        )
        f = ROOT.TFile(input_pattern%fname)
        t = f.Get('tree')
        for e in t:
            for _, evt_set, root_mask in region_tuple:
                if e.evt in evt_set:
                    root_mask.Enter(e.GetReadEntry())

        for rn, _, root_mask in region_tuple:
            print rn, root_mask.GetN()

        return dict(
            (str(rn), root_mask)
            for rn, _, root_mask in region_tuple
        )

    return dict(
        (str(fn), get_teventlist(fn, single_file_json))
        for fn, single_file_json in evt_json.iteritems()
    )


def make_cat_dict(token, evt_json):
    # normalise file/sample names:
    evt_json = dict(
        (k.replace('sum_', '').replace('.root',''), v)
        for k,v in evt_json.iteritems()
    )

    def tree_prep_namespace(lists, fname, cname):
        return lambda tree: tree_prep(tree, lists, fname, cname)

    # teventlists = make_teventlists(evt_json)
    return dict(
        (
            token+'___'+fname+'___'+cname,
            [
                {token+'___'+cname:'1'},
                tree_prep_namespace(evt_json, fname, cname),
                plot_vars_2lep
            ]
        )
        for fname in the_samples_dict
        for cname in cat_names
        if evt_json[fname][cname]
    )


import json
with open('sync_eventlists.json') as f:
    all_events_in_AT = json.load(f)
with open('sync_events_only_in_AT.json') as f:
    events_only_in_AT = json.load(f)
with open('sync_events_only_in_Xb.json') as f:
    events_only_in_Xb = json.load(f)


the_category_dict = {}
the_category_dict.update(make_cat_dict('OnlyInXb', events_only_in_Xb))
#the_category_dict.update(make_cat_dict('OnlyInAT', events_only_in_AT))
#the_category_dict.update(make_cat_dict('AllInAT', all_events_in_AT))

