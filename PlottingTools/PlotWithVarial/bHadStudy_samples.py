from main_samples import *

the_samples_dict = dict(
    (k, v)
    for k, v in the_samples_dict.iteritems()
    if k.startswith('W_')
)
