#!/usr/bin/env python
import os
masses = [1675, 1685, 1695, 1705, 1715, 1725, 1735, 1745, 1755, 1765, 1775]
cmds = "source /cvmfs/cms.cern.ch/cmsset_default.sh;"
for m in masses:
    cmds += "combine -M MultiDimFit cardMorph.txt --redefineSignalPOIs r,MT --floatOtherPOIs 1 -S 0 -t 1000 --expectSignal 1 -s -1 -m %d --setParameters MT=%d;" % (m,m)

os.system(cmds)
