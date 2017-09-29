root -b -q 'jets_Q2.C (313063, 0.03, true)'&
root -b -q 'jets_Q2.C (313067, 1.24, true)'&
root -b -q 'jets_Q2.C (313100, 9.66, true)'&
root -b -q 'jets_Q2.C (313107, 11.92, true)'&
root -b -q 'jets_Q2.C (313136, 10.40, true)'&
root -b -q 'jets_Q2.C (313187, 3.67, true)'&
root -b -q 'jets_Q2.C (313259, 5.12, true)'&
root -b -q 'jets_Q2.C (313285, 4.74, true)'&
root -b -q 'jets_Q2.C (313295, 10.69, true)'&
root -b -q 'jets_Q2.C (313333, 4.13, true)'&
root -b -q 'jets_Q2.C (313435, 0.39, true)'&

wait

root -b -q 'jets_Q2.C (313572, 0.01, false)'&
root -b -q 'jets_Q2.C (313574, 1.33, false)'&
root -b -q 'jets_Q2.C (313575, 7.54, false)'&
root -b -q 'jets_Q2.C (313603, 8.69, false)'&
root -b -q 'jets_Q2.C (313629, 6.86, false)'&
root -b -q 'jets_Q2.C (313630, 7.9, false)'&
root -b -q 'jets_Q2.C (313695, 4.53, false)'&
root -b -q 'jets_Q2.C (313833, 5.11, false)'&
root -b -q 'jets_Q2.C (313878, 2.16, false)'&
root -b -q 'jets_Q2.C (313929, 0.63, false)'&
root -b -q 'jets_Q2.C (313935, 10.96, false)'&
root -b -q 'jets_Q2.C (313984, 2.40, false)'&
root -b -q 'jets_Q2.C (314014, 7.36, false)'&
root -b -q 'jets_Q2.C (314112, 10.49, false)'&
root -b -q 'jets_Q2.C (314157, 9.83, false)'&
root -b -q 'jets_Q2.C (314170, 4.92, false)'&

wait

root -b -q 'jets_Q2_hist.C ({313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314112, 314157, 314170})'
