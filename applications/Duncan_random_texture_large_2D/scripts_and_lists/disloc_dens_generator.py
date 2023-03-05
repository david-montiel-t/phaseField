import numpy as np
N=141
N_grains=91

rho_ind = np.zeros(N_grains)
rho_units = np.zeros(N_grains)

f1 = open("new_graindds.txt")
for j in range(0,N_grains):
    rho_ind[j], rho_units[j] = f1.readline().split()
f1.close()

rho_mean = np.mean(rho_units)
rho_std_dev = np.std(rho_units)

# Rescaled dislocation density
rho_i = rho_units/rho_mean

rho_i = np.append(rho_i,np.zeros(N-N_grains))

# Writing the vector for densities
print ('Average dislocartion density: ' + str(rho_mean))
print ('Std dev: ' + str(rho_std_dev))
print ('Normalized Std_ dev: ' + str(rho_std_dev/rho_mean))

f1 = open("dislocation_densities.txt",'w')
f1.write('{')
for j in range(N):
    if j+1<N:
        if (j+1)%4:
            f1.write('%.5f' % rho_i[j] +', ')
        else:
            f1.write('%.5f' % rho_i[j]+ ', '+'\n')
    if j+1==N:
        f1.write('%.5f' % rho_i[j])
f1.write('};')
f1.close()
