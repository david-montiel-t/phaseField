
# Works for periodic BCs ONLY
# Generates a list of subsections for the adaptivity subsection for "n"

first_nuc_ind = 92
last_nuc_ind = 141

fl = open("nuc_dep_list.dat",'w')

for j in range(first_nuc_ind,last_nuc_ind+1):
    fl.write("subsection Nucleation parameters: n"+str(j)+'\n')
    fl.write("set Nucleus semiaxes (x, y, z) = 0.7, 0.7, 0.7\n")
    fl.write("set Freeze zone semiaxes (x, y, z) = 2.0, 2.0, 2.0\n")
    fl.write("set Nucleus rotation in degrees (x, y, z) = 0, 0, 0\n")
    fl.write("set Freeze time following nucleation = 2\n")
    fl.write("set Nucleation-free border thickness = 0\n")
    fl.write("end\n\n")
fl.close()
