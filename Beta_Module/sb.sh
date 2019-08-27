f90 $F90FLAGS -o run -s -w super_beta_muon.f $LINK_FNL
./run
rm run 
echo "Done!"
