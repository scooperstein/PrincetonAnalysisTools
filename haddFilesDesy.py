#!/bin/env python

import sys, os
join = os.path.join

if len(sys.argv) < 2:
    print 'I need the path to job output directory...'
    exit(-1)

base_dir = sys.argv[1]
out_dir = join(base_dir, 'haddjobs')
print '>    output goes to', out_dir
os.system('mkdir -p ' + out_dir)

_, sub_dirs, __ = next(os.walk(base_dir))
sub_dirs.remove('haddjobs')

for sub_dir in sub_dirs:
    cmd = 'hadd -f %s/sum_%s.root ' % (out_dir, sub_dir)
    cmd += join(base_dir, sub_dir, '*.root')
    print '>    issuing', cmd
    os.system(cmd)

cmd = 'hadd -f {od}/sum_Data.root {od}/sum_Run*.root'.format(od=out_dir)
print '>    issuing', cmd
os.system(cmd)
