#!/usr/bin/env python

# compare differences between pdb directories and use
# upload-single-to-box.sh to upload only different files to Box.

import sys
import os
import re
import subprocess
import time

def upload(f):
    tries = 0
    done = False
    while tries < 20 and not done:
        try:
            o = subprocess.check_output(['./upload-single-to-box.sh', f])
            print o
            done = True
        except:
            time.sleep(5)
            tries += 1

files = os.listdir('.')
weeks = [x for x in files if re.match(r'pdb.\d+', x)]
for i in range(len(weeks)):
    if i < len(weeks) - 1:
        old = weeks[i]
        new = weeks[i+1]

        print 'processing', old, new

        p = subprocess.Popen(['diff', '-r', '--brief', old, new], 
            stdout=subprocess.PIPE)
        for line in iter(p.stdout.readline, b''):
            m = re.match(r'Only in ([^\s]+): ([^\s]+)', line)
            if m:
                f = '{}/{}'.format(m.group(1), m.group(2))
                upload(f)
            else:
                m = re.match(r'Files ([^\s]+) and ([^\s]+) differ', line)
                if m:
                    new_f = m.group(2)
                    upload(new_f)
                else:
                    print 'WARNING: unmatched', line

        p.stdout.close()
        p.wait()
