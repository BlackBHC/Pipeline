# For hdf5 to txt file:
python hdf52txt.py -id=$VMGIC/test_ic -od=$VMGIC/test_ic -if=diskhalo -of=model -ir=$INI

# For radial profile
python radial_profile.py -id=$VMGIC/Fig -od=$VMGIC/Fig -if=model -of=radial_profile -ir=$INI
