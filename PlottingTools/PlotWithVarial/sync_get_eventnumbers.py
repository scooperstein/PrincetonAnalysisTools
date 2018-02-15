#!/usr/bin/env python

import sync_config
import ROOT
import glob
import json
import uuid
import os

cats = {}
cats.update(sync_config.the_category_dict['CR_2LeptLO'][0])
cats.update(sync_config.the_category_dict['CR_2LeptHI'][0])
basic_sel = sync_config.the_category_dict['CR_2LeptHI'][1]
basic_sel = '&&'.join(basic_sel)

def make_eventlist_for_cat(filepath, selection):
    f = ROOT.TFile(filepath)
    tree = f.Get('tree')
    selection = '(%s)&&%s' % (selection, basic_sel)

    name = uuid.uuid4().hex
    tree.Draw('>>{0}'.format(name), selection, 'goff')
    eventlist = ROOT.gDirectory.Get(name)
    if not eventlist.GetN():
        return []

    tree.SetEventList(eventlist)
    evt_list = list(e.evt for e in tree if eventlist.Contains(e.GetReadEntry()))

    return evt_list


def make_eventlists_for_file(filepath):
    print filepath
    return dict(
        (catname, make_eventlist_for_cat(filepath, selection))
        for catname, selection in cats.iteritems()
    )


if __name__ == '__main__':
    syncfiles = glob.glob(sync_config.input_pattern%'*')
    eventlists = dict(
        (os.path.basename(f), make_eventlists_for_file(f))
        for f in syncfiles
    )
    with open('sync_eventlists.json', 'w') as f:
        json.dump(eventlists, f, indent=4, sort_keys=True)
