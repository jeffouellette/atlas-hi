#!/bin/bash

# usage: ./zgen SEED PTHATMIN MINZPT MAXZPT MINTRKPT MAXTRKPT NEVT FILENAMEOUT


echo "Universe = vanilla" > submit.job
echo "Notification = Complete" >> submit.job
echo "Executable = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/PullStudy/bin/gen" >> submit.job
echo "Priority = +1" >> submit.job
echo "GetEnv = True" >> submit.job
echo "InitialDir = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/PullStudy" >> submit.job
echo "accounting_group = group_atlas.boulder" >> submit.job

echo "" >> submit.job
echo "" >> submit.job

njobs=0

for iPtZ in {2..4}
do
  for iCent in {1..3}
  do
    for iSeed in {1001..2000}
    do
      seed=$((iSeed+10000*(3*(iPtZ-2)+(iCent-1))))
      echo "Arguments = ${seed} ${iPtZ} ${iCent} /atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/PullStudy/iSeed${iSeed}_iPtZ${iPtZ}_iCent${iCent}.root" >> submit.job
      echo "Log = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/PullStudy/logs/iSeed${iSeed}_iPtZ${iPtZ}_iCent${iCent}.log" >> submit.job
      echo "Output = /atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/outputs/PullStudy/iSeed${iSeed}_iPtZ${iPtZ}_iCent${iCent}.out" >> submit.job
      echo "Error = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/PullStudy/errors/iSeed${iSeed}_iPtZ${iPtZ}_iCent${iCent}.err" >> submit.job
      echo "Queue 1" >> submit.job
      echo "" >> submit.job

      njobs=$((njobs+1))
    done
  done
done

echo "${njobs} jobs setup for submission"
