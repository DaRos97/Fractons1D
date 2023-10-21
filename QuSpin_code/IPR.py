import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path
import functions as fs

sps, N, Np = fs.input_routine()
#### Load Eigenstates
dirname = Path('.')
text_s = "spin-" + fs.spin_text[sps-2] + "_short_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
text_l = "spin-" + fs.spin_text[sps-2] + "_long_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
path_to_load = dirname / "Data" / "Evs" / text_s
ev_3 = np.load(path_to_load)
path_to_load = dirname / "Data" / "Evcs" / text_s
evc_3 = np.load(path_to_load)
path_to_load = dirname / "Data" / "Evs" / text_l
ev_34 = np.load(path_to_load)
path_to_load = dirname / "Data" / "Evcs" / text_l
evc_34 = np.load(path_to_load)
### Calculate IPR
IPR3 = np.zeros(len(ev_3))
IPR34 = np.zeros(len(ev_34))
for e in range(len(ev_3)):
    IPR3[e] = (np.abs(evc_3[:,e])**4).sum()
for e in range(len(ev_34)):
    IPR34[e] = (np.abs(evc_34[:,e])**4).sum()
#### Plot
plt.subplots(1,2,figsize=(12,4), sharey=True)
plt.subplot(1,2,1)
bins=(51,int(50*max(IPR3)))
plt.hist2d(ev_3,
           IPR3,
           bins,
           norm=colors.LogNorm())
plt.xlabel('E')
plt.ylabel('IPR')
plt.colorbar()
plt.title('IPR of $H_3$')
plt.ylim(0,1)
#
bins=(51,int(50*max(IPR34)))
plt.subplot(1,2,2)
plt.hist2d(ev_34,
           IPR34,
           bins,
           norm=colors.LogNorm())
plt.xlabel('E')
plt.ylabel('IPR')
plt.colorbar()
plt.title('IPR of $H_{34}$')
plt.ylim(0,1)

plt.show()
