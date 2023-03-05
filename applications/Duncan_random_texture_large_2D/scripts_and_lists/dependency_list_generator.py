# Generates a list of order parameters separated by commas and quotation signs "n"

# Number of order parameters
nops = 141

# Writing out positions
fl = open("list_of_ops.txt",'w')
fl.write(" List of order parameters"+'\n"')

for j in range(1,nops+1):
    fl.write('n'+str(j))
    if j==nops:
        fl.write('"')
    else:
        fl.write(', ')
fl.close()
