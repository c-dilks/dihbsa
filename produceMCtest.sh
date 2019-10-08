#!/bin/bash
loopAsym.sh outroot.MC.gen -m2 && mv spinFinal.root forMC/gen.${1}.root
loopAsym.sh outroot.MC.rec -m2 && mv spinFinal.root forMC/rec.${1}.root
