root -b -q 'jets_Q2.C (313067, 0.929)'&
root -b -q 'jets_Q2.C (313107, 10.843)'&
root -b -q 'jets_Q2.C (313136, 9.42)'&
root -b -q 'jets_Q2.C (313187, 3.24)'&

root -b -q 'jets_Q2.C (313259, 4.634)'&
root -b -q 'jets_Q2.C (313572, 0.0051)'&
root -b -q 'jets_Q2.C (313574, 1.2)'&
root -b -q 'jets_Q2.C (313575, 7.059)'&

root -b -q 'jets_Q2.C (313603, 8.158)'&
root -b -q 'jets_Q2.C (313629, 6.251)'&
root -b -q 'jets_Q2.C (313630, 6.632)'&
root -b -q 'jets_Q2.C (313695, 4.139)'&

root -b -q 'jets_Q2.C (313833, 4.70)'&
root -b -q 'jets_Q2.C (313878, 1.784)'&
root -b -q 'jets_Q2.C (313929, 0.0506)'&
root -b -q 'jets_Q2.C (314014, 6.850)'&

root -b -q 'jets_Q2.C (314157, 9.153)'&

wait

root -b -q 'jets_Q2_hist.C ({313067, 313107, 313136, 313187, 313259, 313572, 313574, 313575, 313603, 313629, 313630, 313695, 313833, 313878, 313929, 314014, 314157})'
