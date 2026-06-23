cd src

idl << EOF 

.RESET_SESSION
.compile oisst_daily_nc2hdf5.pro
Resolve_All
save,filename="../oisst_daily_nc2hdf5.sav", /routines

exit
EOF

cd ../
chmod 774 oisst_daily_nc2hdf5.sav

