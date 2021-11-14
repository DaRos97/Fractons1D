import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import functions as fs

sps, N, Np = fs.input_routine()

#### Comparison Functions for spacings and ratios
def P(x):       #Poisson spacing
    return np.exp(-x)
def Pr(x):      #Poisson ratio
    return 2/(1+x)**2
def WD(x):      #Wigner-Dyson spacing
    return np.pi/2*x*np.exp(-np.pi/4*x**2)
def WDr(x):     #Wigner-Dyson ratio
    return (27*(x+x**2))/(4*np.power((1+x+x**2),2.5))
def SP(x):      #Semi-Poisson spacing
    return 4*x*np.exp(-2*x)
def SPr(x):     #Semi-Poisson ratio
    return 12*x/(1+x)**4

dirname = Path('.')
text_s = "spin-" + fs.spin_text[sps-2] + "_short_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
text_l = "spin-" + fs.spin_text[sps-2] + "_long_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
path_to_load = dirname / "Data" / "Evs" / text_s
ev3 = np.load(path_to_load)
path_to_load = dirname / "Data" / "Evs" / text_l
ev34 = np.load(path_to_load)
#### Remove degeneracies
lim = 1e-10
rev3 = [ev3[0]]
for i in range(len(ev3)-1):
    if np.abs(ev3[i+1]-ev3[i])>lim:
        rev3 += [ev3[i+1]]
rev3 = np.array(rev3)
#
rev34 = [ev34[0]]
for i in range(len(ev34)-1):
    if np.abs(ev34[i+1]-ev34[i])>lim:
        rev34 += [ev34[i+1]]
rev34 = np.array(rev34)

#### Evaluate statistics inside a given range
## eps is the normalized spectrum between 0 and 1
eps3 = [(e-rev3[0])/(rev3[-1]-rev3[0]) for e in rev3]
eps34 = [(e-rev34[0])/(rev34[-1]-rev34[0]) for e in rev34]
## Steps to analyze inside eps
pts = 15            #Lower has less fluctuation
step = np.linspace(0.3,0.7,pts)         #this choses the range of the spectrum
steps = [(step[i+1]+step[i])/2 for i in range(pts-1)]
## s and r 
el3 = []
for e in eps3:
    if e > step[0] and e < step[-1]:
        el3 += [e]
s3 = [el3[i+1]-el3[i] for i in range(len(el3)-1)]
s3 = np.array(s3)
s3 /= (s3.sum()/len(s3))
r3 = [min(s3[i],s3[i+1])/max(s3[i],s3[i+1]) for i in range(len(s3)-1)]
#
el34 = []
for e in eps34:
    if e > step[0] and e < step[-1]:
        el34 += [e]
s34 = [el34[i+1]-el34[i] for i in range(len(el34)-1)]
s34 = np.array(s34)
s34 /= (s34.sum()/len(s34))
r34 = [min(s34[i],s34[i+1])/max(s34[i],s34[i+1]) for i in range(len(s34)-1)]
## r in step
rf3 = []
for i in range(pts-1):
    el = []
    for e in eps3:
        if e > step[i] and e < step[i+1]:
            el += [e]
    if len(el)==0:
        continue
    st = [el[i+1]-el[i] for i in range(len(el)-1)]
    st = np.array(st)
    meant = st.sum()/len(st)
    st /= meant
    rt = [min(st[i],st[i+1])/max(st[i],st[i+1]) for i in range(len(st)-1)]
    rt = np.array(rt)
    rf3 += [rt.sum()/len(rt)]
rf3 = np.array(rf3)
rf3m = rf3.sum()/len(rf3)
#
rf34 = []
for i in range(pts-1):
    el = []
    for e in eps34:
        if e > step[i] and e < step[i+1]:
            el += [e]
    if len(el)==0:
        continue
    st = [el[i+1]-el[i] for i in range(len(el)-1)]
    st = np.array(st)
    meant = st.sum()/len(st)
    st /= meant
    rt = [min(st[i],st[i+1])/max(st[i],st[i+1]) for i in range(len(st)-1)]
    rt = np.array(rt)
    rf34 += [rt.sum()/len(rt)]
rf34 = np.array(rf34)
rf34m = rf34.sum()/len(rf34)

#### Plot
plt.subplots(3,2,figsize=(12,9))
## s
bins=70
plt.subplot(3,2,1)
plt.hist(s3,bins=bins, density=1, alpha=0.5)
plt.plot(np.linspace(0,max(s3),1000),WD(np.linspace(0,max(s3),1000)), 'r', label='WD')
plt.plot(np.linspace(0,max(s3),1000),SP(np.linspace(0,max(s3),1000)), 'k', label='SP')
plt.plot(np.linspace(0,max(s3),1000),P(np.linspace(0,max(s3),1000)), 'g', label='P')
plt.xlim(0,4)
plt.legend(title='$H_3$, $L=%i$'%N)
plt.xlabel('spacing')
#
plt.subplot(3,2,2)
plt.hist(s34,bins=bins, density=1, alpha=0.5)
plt.plot(np.linspace(0,max(s34),1000),WD(np.linspace(0,max(s34),1000)), 'r', label='WD')
plt.plot(np.linspace(0,max(s34),1000),SP(np.linspace(0,max(s34),1000)), 'k', label='SP')
plt.plot(np.linspace(0,max(s34),1000),P(np.linspace(0,max(s34),1000)), 'g', label='P')
plt.xlim(0,4)
plt.legend(title='$H_{34}$, $L=%i$'%N)
plt.xlabel('spacing')
## r
bins=50
plt.subplot(3,2,3)
plt.hist(r3,bins=bins, density=1, alpha=0.5)
plt.plot(np.linspace(0,max(r3),1000),WDr(np.linspace(0,max(r3),1000)), 'r', label='WD')
plt.plot(np.linspace(0,max(r3),1000),SPr(np.linspace(0,max(r3),1000)), 'k', label='SP')
plt.plot(np.linspace(0,max(r3),1000),Pr(np.linspace(0,max(r3),1000)), 'g', label='P')
plt.xlim(0,1)
plt.legend(title='$H_3$, $L=%i$'%N)
#ax[0].text(0.3,1.8,'$\\langle r \\rangle = %f$' % r3m, size='x-large')
plt.xlabel('ratio')
#
plt.subplot(3,2,4)
plt.hist(r34,bins=bins, density=1, alpha=0.5)
plt.plot(np.linspace(0,max(r34),1000),WDr(np.linspace(0,max(r34),1000)), 'r', label='WD')
plt.plot(np.linspace(0,max(r34),1000),SPr(np.linspace(0,max(r34),1000)), 'k', label='SP')
plt.plot(np.linspace(0,max(r34),1000),Pr(np.linspace(0,max(r34),1000)), 'g', label='P')
plt.xlim(0,1)
plt.legend(title='$H_{34}$, $L=%i$'%N)
plt.xlabel('ratio')
## Steps
plt.subplot(3,2,5)
plt.plot(steps,rf3,'*')
plt.hlines(rf3m,step[0],step[-1],'b',linestyles='dashed', label='$\\langle r \\rangle = %.3f$' %rf3m)
plt.hlines(0.536,step[0],step[-1],'r',linestyles='dashed', label='WD')
plt.hlines(0.5,step[0],step[-1],'k',linestyles='dashed', label='SP')
plt.hlines(0.386,step[0],step[-1],'g',linestyles='dashed', label='P')
plt.legend(title='$H_3$, $L=%i$'%N)#, loc='upper right', bbox_to_anchor=(0.55, 0.42, 0.85, 0.61))
plt.ylabel('r')
plt.xlabel('epsilon')
plt.xlim(step[0],step[-1]+0.2)
#
plt.subplot(3,2,6)
plt.plot(steps,rf34,'*')
plt.hlines(rf34m,step[0],step[-1],'b',linestyles='dashed', label='$\\langle r \\rangle = %.3f$' % rf34m)
plt.hlines(0.536,step[0],step[-1],'r',linestyles='dashed', label='WD')
plt.hlines(0.5,step[0],step[-1],'k',linestyles='dashed', label='SP')
plt.hlines(0.386,step[0],step[-1],'g',linestyles='dashed', label='P')
plt.legend(title='$H_{34}$, $L=%i$'%N)#, loc='upper right', bbox_to_anchor=(0.55, 0.42, 0.85, 0.61))
plt.xlabel('epsilon')
plt.xlim(step[0],step[-1]+0.2)
#plt.savefig('../Figs/stat_r', bbox_inches='tight')
plt.show()

