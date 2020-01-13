#!/bin/bash

# usage: ./zgen SEED PTHATMIN MINZPT MAXZPT MINTRKPT MAXTRKPT NEVT FILENAMEOUT

#for iSeed in {80..159}
#do
#  time ./bin/gen ${iSeed} 5 15 30 657 output/seed${iSeed}_iPtZ2_iCent1.root > log/seed${iSeed}_iPtZ2_iCent1.log 2> errors/seed${iSeed}_iPtZ2_iCent1.err &
#  if [[ $((iSeed % 8)) -eq 7 ]]
#  then
#    wait
#  fi
#done

#for iSeed in {80..159}
#do
#  time ./bin/gen $((iSeed+1000)) 5 15 30 1501 output/seed${iSeed}_iPtZ2_iCent2.root > log/seed${iSeed}_iPtZ2_iCent2.log 2> errors/seed${iSeed}_iPtZ2_iCent2.err &
#  if [[ $((iSeed % 8)) -eq 7 ]]
#  then
#    wait
#  fi
#done

#for iSeed in {80..159}
#do
#  time ./bin/gen $((iSeed+2000)) 5 15 30 1587 output/seed${iSeed}_iPtZ2_iCent3.root > log/seed${iSeed}_iPtZ2_iCent3.log 2> errors/seed${iSeed}_iPtZ2_iCent3.err &
#  if [[ $((iSeed % 8)) -eq 7 ]]
#  then
#    wait
#  fi
#done

#for iSeed in {80..159}
#do
#  time ./bin/gen $((iSeed+3000)) 5 30 60 354 output/seed${iSeed}_iPtZ3_iCent1.root > log/seed${iSeed}_iPtZ3_iCent1.log 2> errors/seed${iSeed}_iPtZ3_iCent1.err &
#  if [[ $((iSeed % 8)) -eq 7 ]]
#  then
#    wait
#  fi
#done

for iSeed in {80..159}
do
  time ./bin/gen $((iSeed+4000)) 5 30 60 848 output/seed${iSeed}_iPtZ3_iCent2.root > log/seed${iSeed}_iPtZ3_iCent2.log 2> errors/seed${iSeed}_iPtZ3_iCent2.err &
  if [[ $((iSeed % 8)) -eq 7 ]]
  then
    wait
  fi
done

for iSeed in {80..159}
do
  time ./bin/gen $((iSeed+5000)) 5 30 60 849 output/seed${iSeed}_iPtZ3_iCent3.root > log/seed${iSeed}_iPtZ3_iCent3.log 2> errors/seed${iSeed}_iPtZ3_iCent3.err &
  if [[ $((iSeed % 8)) -eq 7 ]]
  then
    wait
  fi
done

for iSeed in {80..159}
do
  time ./bin/gen $((iSeed+6000)) 30 60 10000 141 output/seed${iSeed}_iPtZ4_iCent1.root > log/seed${iSeed}_iPtZ4_iCent1.log 2> errors/seed${iSeed}_iPtZ4_iCent1.err &
  if [[ $((iSeed % 8)) -eq 7 ]]
  then
    wait
  fi
done

for iSeed in {80..159}
do
  time ./bin/gen $((iSeed+7000)) 30 60 10000 318 output/seed${iSeed}_iPtZ4_iCent2.root > log/seed${iSeed}_iPtZ4_iCent2.log 2> errors/seed${iSeed}_iPtZ4_iCent2.err &
  if [[ $((iSeed % 8)) -eq 7 ]]
  then
    wait
  fi
done

for iSeed in {80..159}
do
  time ./bin/gen $((iSeed+8000)) 30 60 10000 318 output/seed${iSeed}_iPtZ4_iCent2.root > log/seed${iSeed}_iPtZ4_iCent3.log 2> errors/seed${iSeed}_iPtZ4_iCent3.err &
  if [[ $((iSeed % 8)) -eq 7 ]]
  then
    wait
  fi
done

