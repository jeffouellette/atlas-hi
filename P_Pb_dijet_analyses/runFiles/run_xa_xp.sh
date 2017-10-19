#root -b -q '../jet_analyses/jets_xa_xp.C (313063, 0.02736, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313067, 1.243, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313100, 9.66, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313107, 11.92, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313136, 10.4, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313187, 3.67, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313259, 5.124, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313285, 4.744, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313295, 10.69, true)'&
#
#wait
#
#root -b -q '../jet_analyses/jets_xa_xp.C (313333, 4.133, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313435, 0.3886, true)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313572, 0.006006, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313574, 1.334, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313575, 7.538, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313603, 8.687, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313629, 6.856, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313630, 7.903, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313695, 4.527, false)'&
#
#wait
#
#root -b -q '../jet_analyses/jets_xa_xp.C (313833, 5.113, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313878, 2.163, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313929, 0.6315, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313935, 10.96, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (313984, 2.404, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (314014, 7.36, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (314112, 10.49, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (314157, 9.83, false)'&
#root -b -q '../jet_analyses/jets_xa_xp.C (314170, 4.919, false)'&
#
#wait

root -b -q '../jet_analyses/jets_xa_xp_hist.C ({313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314112, 314157, 314170})'
