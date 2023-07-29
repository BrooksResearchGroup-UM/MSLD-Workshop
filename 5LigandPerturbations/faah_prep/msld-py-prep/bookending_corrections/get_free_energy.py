import numpy as np
from pymbar import MBAR
import seaborn as sns
import matplotlib.pyplot as plt

temp = 298.15
beta = 1/(0.001987204259*temp)
 
# Assemble matrix in reduced energies
ffE1 = np.loadtxt('ff_e1.csv',dtype=float)
crnE1 = np.loadtxt('ff_e2.csv',dtype=float)

crnE2 = np.loadtxt('crn_e1.csv',dtype=float)
ffE2 = np.loadtxt('crn_e2.csv',dtype=float)

ffE = np.concatenate((ffE1, ffE2))
crnE = np.concatenate((crnE1, crnE2))

u_kn = np.stack((ffE,crnE),axis=0)*beta
N_k = np.array([len(ffE1),len(crnE2)])

# Instantiate MBAR object
mbar = MBAR(u_kn,N_k)

result = mbar.compute_free_energy_differences()

# Get free energy/uncertainty in kcal/mol
fe = result['Delta_f'][0,1]/beta
unc = result['dDelta_f'][0,1]/beta

# Write to file
with open('Result.txt','w') as f:
    f.write(f'{fe:.4f} +/- {unc:.4f}')

# Get overlap
overlap = mbar.compute_overlap()
overlap = overlap['matrix'][0,0]

# Write overlap to file
with open('overlap.txt','w') as f:
    f.write(f'{overlap*100}:.3f')

# Plot and save figure
sns.set_style('darkgrid')
sns.histplot(crnE1-ffE1,color='b',kde=True)
sns.histplot(crnE2-ffE2,color='r',kde=True)

plt.xlabel(r'$\Delta$U (kcal/mol)',fontsize=15)
plt.ylabel('Count',fontsize=15)
plt.title(f'{overlap*100:.2f}% Overlap',fontsize=20)
plt.legend(['FF->CRN','CRN->FF'])
plt.savefig('overlap.pdf')

