
from main_samples import *

# for plotting, samples are merged using the legend entry
# example: merge flavours in EWK samples
for samplename, items in the_samples_dict.iteritems():
    if samplename.startswith('W_'):
        items[2] = 'W+jets'   # set a new legend name
    if samplename.startswith('Z_'):
        items[2] = 'Z+jets'   # set a new legend name
