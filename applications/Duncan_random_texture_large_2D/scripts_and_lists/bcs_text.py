ngrains = 141

fl = open("bcs_text.dat",'w')
for j in range(1,ngrains+1):
    fl.write("set Boundary condition for variable n"+str(j)+' = NATURAL\n')
fl.close()
