from pathlib import Path
dirname = Path('.')     #os.path.dirname(os.path.abspath(__file__))
import QuSpin_functions as qs
import numpy as np
from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis.user import user_basis # Hilbert space user basis
import matplotlib.pyplot as plt
import functions as fs

sps, N, Np = fs.input_routine()
no_checks = dict(check_symm=False, check_pcon=False, check_herm=False)
states_text = "spin-" + fs.spin_text[sps-2] + \
    "_N=%i_(q,p)=(%i,%i).npy" % (N, Np[0], Np[1])
path_to_load = dirname / "Data" / "States" / states_text
states = np.array(np.load(path_to_load), dtype=np.uint32)
#### Function args and dictionaries
op_args = np.array([sps], dtype=np.uint32)
next_state_args = states
P_args = np.array([sps], dtype=np.uint32)
maps = dict(P_block=(qs.parity, 2, 0, P_args))
pcon_dict = dict(Np=Np, next_state=qs.next_state, next_state_args=next_state_args,
                 get_Ns_pcon=qs.get_Ns_pcon, get_s0_pcon=qs.get_s0_pcon)
op_dict = dict(op=qs.op, op_args=op_args)
allowed_ops = ["+-zi", "+-zPMZ"]
#### Basis
basis = user_basis(np.uint32, N, op_dict, allowed_ops=set(allowed_ops[sps-2]), 
                   sps=sps, pcon_dict=pcon_dict, parallel=True, **maps)
#### Lists
J_s = fs.interaction_list(N, fs.spin_ranges[sps-2][0])
J_l = fs.interaction_list(N, fs.spin_ranges[sps-2][1])
if sps == 2:
    static_s = [[op, J_s] for op in ["+--+", "-++-"]]
    static_l = [["+--+", J_s], ["-++-", J_s], ["+-i-+", J_l], ["-+i+-", J_l]]
elif sps == 3:
    static_s = [[op, J_s] for op in ["+M+", "-P-"]]
    static_l = [["+M+", J_s], ["-P-", J_s], ["+--+", J_l], ["-++-", J_l]]
#### Hamiltonian operators
no_checks = dict(check_symm=False, check_pcon=False, check_herm=False)
H_s = hamiltonian(static_s, [], basis=basis, dtype=np.float64, **no_checks)
H_l = hamiltonian(static_l, [], basis=basis, dtype=np.float64, **no_checks)
#####
##### Autocorrelator
#####
times = np.logspace(-2,2,100,base=10.0)
C_3 = qs.auto_correlator(times,H_s,N,basis)
C_34 = qs.auto_correlator(times,H_l,N,basis)
av3 = 0
av34 = 0
tott = 0
for t in range(len(times)):
    if times[t]>5:
        av3+=abs(C_3[t])
        av34+=abs(C_34[t])
        tott+=1
av3/=tott
av34/=tott
##### Plot
plt.figure(figsize=(6,3))
plt.plot(times,C_3,'y-',label='$H_3$')
plt.hlines(av3,0,100,'c',linestyles='dashed',label='Long-time average $H_3$')
plt.plot(times,C_34,'-r',label='$H_{34}$')
plt.hlines(av34,0,100,'b',linestyles='dashed',label='Long-time average $H_{34}$')
plt.legend()
plt.xscale('log')
plt.xlabel('time')
plt.ylabel('$\\langle C_0^z(t) \\rangle$')
#plt.savefig('Figs/autocorr_%i_(%i,%i)'%(N,Np[0],Np[1]))
plt.show()

