import numpy as np
from quspin.basis.user import user_basis  # Hilbert space user basis
from quspin.operators import hamiltonian  # Hamiltonians and operators
from pathlib import Path
dirname = Path('.')
import functions as fs
import QuSpin_functions as qs

#### Input spin, N and Np
sps,N,Np = fs.input_routine()
# Load states
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
                   sps=sps, pcon_dict=pcon_dict, parallel=True)
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
#### Diagonalization
while True:
    also_evc = int(
        input("Compute all eigenstates(1) or only eigenvaues(0)? [1,0]: "))
    if also_evc != 0 and also_evc != 1:
        print("Error in input, retry")
        continue
    break
if also_evc:
    ev_s, evc_s = H_s.eigh()
    ev_l, evc_l = H_l.eigh()
else:
    ev_s = H_s.eigvalsh()
    ev_l = H_l.eigvalsh()
#### Save Data
text_s = "spin-" + fs.spin_text[sps-2] + "_short_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
text_l = "spin-" + fs.spin_text[sps-2] + "_long_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
path_to_save = dirname / "Data" / "Evs" / text_s
np.save(path_to_save, ev_s)
path_to_save = dirname / "Data" / "Evs" / text_l
np.save(path_to_save, ev_l)
if also_evc:
    path_to_save = dirname / "Data" / "Evcs" / text_s
    np.save(path_to_save, evc_s)
    path_to_save = dirname / "Data" / "Evcs" / text_l
    np.save(path_to_save, evc_l)
