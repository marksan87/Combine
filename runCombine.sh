#!/bin/bash
combine -M MultiDimFit cardMorph.txt --redefineSignalPOIs r,MT -P MT --floatOtherPOIs 1 -S 1 --expectSignal 1 --algo grid --points 100 --setParameterRanges MT=1665,1785 --robustFit 1 --saveWorkspace --saveFitResult
