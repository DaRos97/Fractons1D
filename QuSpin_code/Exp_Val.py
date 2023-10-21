from pathlib import Path
dirname = Path('.')
import QuSpin_functions as qs
import functions as fs
from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis.user import user_basis # Hilbert space user basis
import numpy as np
import matplotlib.pyplot as plt


sps, N , Np = fs.input_routine()
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
# ## Import Ev and Evc
# text_s = "spin-" + sf.spin_text[sps-2] + "_short_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
# text_l = "spin-" + sf.spin_text[sps-2] + "_long_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
# path_to_load = dirname / "Data" / "Evs" / text_s
# ev_s = np.load(path_to_load)
# path_to_load = dirname / "Data" / "Evs" / text_l
# ev_l = np.load(path_to_load)
# path_to_load = dirname / "Data" / "Evcs" / text_s
# evc_s = np.load(path_to_load)
# path_to_load = dirname / "Data" / "Evcs" / text_l
# evc_l = np.load(path_to_load)
text_s = "spin-" + fs.spin_text[sps-2] + "_short_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
text_l = "spin-" + fs.spin_text[sps-2] + "_long_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
path_to_load = dirname / "Data" / "Evs" / text_s
ev_s = np.load(path_to_load)
path_to_load = dirname / "Data" / "Evcs" / text_s
evc_s = np.load(path_to_load)
path_to_load = dirname / "Data" / "Evs" / text_l
ev_l = np.load(path_to_load)
path_to_load = dirname / "Data" / "Evcs" / text_l
evc_l = np.load(path_to_load)
#ev_s, evc_s = H_s.eigh()
#ev_l, evc_l = H_l.eigh()
#%%
#### Sz2
Sz2 = hamiltonian([[ "z", [[1.0, N//2]] ]], [], dtype=np.float64, basis=basis, **no_checks)
ExpVal_s=[Sz2.expt_value(evecs) for evecs in evc_s.T]
ExpVal_l=[Sz2.expt_value(evecs) for evecs in evc_l.T]

#### Plot
plt.subplots(1,2,figsize=(12,3),sharey=True)
plt.subplot(1,2,1)
plt.plot(ev_s/N,ExpVal_s,'g.', label = 'N=%i'%N)
plt.legend()

plt.subplot(1,2,2)
plt.plot(ev_l/N,ExpVal_l,'g.', label = 'N=%i'%N)
plt.legend()

plt.show()
