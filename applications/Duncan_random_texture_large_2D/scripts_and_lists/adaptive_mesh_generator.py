
# Works for periodic BCs ONLY
# Generates a list of subsections for the adaptivity subsection for "n"

ngrains = 141

fl = open("adaptive_mesh_text.dat",'w')

for j in range(1,ngrains+1):
    fl.write("subsection Refinement criterion: n"+str(j)+'\n')
    fl.write("# Select whether the mesh is refined based on the variable value (VALUE),\n")
    fl.write("# its gradient (GRADIENT), or both (VALUE_AND_GRADIENT)\n")
    fl.write("set Criterion type = VALUE\n")
    fl.write("# Set the lower and upper bounds for the value-based refinement window\n")
    fl.write("set Value lower bound = 0.001\n")
    fl.write("set Value upper bound = 0.999\n")
    fl.write("end\n\n")
fl.close()
