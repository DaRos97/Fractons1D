import numpy as np
from quspin.operators import hamiltonian
# user basis data types signatures
from quspin.basis.user import next_state_sig_32, op_sig_32, map_sig_32
from numba import carray, cfunc  # numba helper functions
from numba import uint32, int32, float64  # numba data types
from pathlib import Path
dirname = Path('.')
import functions as fs


@cfunc("int32[:](uint32,int32,uint32)", locals=dict(res=int32[:]))
def int_to_list(s, N, sps):
    """
    Transforms an integer in binary or trinary representation depending on the spin value. We use the convention that in binary the single bit can be either +1 or -1, while in trinary can be also 0.

    Parameters
    ----------
    n : uint32
        Number to be transformed.
    N : int32
        Vector length.
    sps : uint32
        Namely 'states per site', is the usual sps = 2*s+1 where s is the spin, which can be either 1/2 or 1.

    Returns
    -------
    el : int32[:]
        Array containing the binary/trinary representation of n.

    """
    res=np.zeros(N,dtype=np.int32)
    for i in range(N):
        res[i]=-1
    if sps == 3:
        for i in range(N-1,-1,-1):
            if 3**i>s:
                res[N-i-1]=-1
                continue
            if 2*3**i>s:
                res[N-i-1]=0
            else:
                res[N-i-1]=1
            s-=(res[N-i-1]+1)*3**i
    elif sps == 2:
        for i in range(N-1,-1,-1):
            if 2**i>s:
                res[N-i-1] = -1
                continue
            else:
                res[N-i-1] = 1
            s -= 2**i
    return res


@cfunc("uint32(uint32,int32,uint32)", locals=dict(el=int32[:], el2=int32[:], res=uint32))
def P(n, N, sps):
    """
    Calculates the parity conjugate of a number n, passing from its binary/trinary representation. Essentially sends +1 to -1 and vice-versa.

    Parameters
    ----------
    n : uint32
        Number to be transformed.
    N : int32
        Vector length .
    sps : uint32
        Namely 'states per site', is the usual sps = 2*s+1 where s is the spin, which can be either 1/2 or 1.

    Returns
    -------
    res : uint32
        Parity transformed number.

    """
    el = int_to_list(n,N,sps)
    el2 = np.zeros(N, dtype=np.int32)
    for i in range(N):
        el2[i] = -el[i]
    res = 0
    if sps == 3:
        for i in range(N):
            res += (sps**(N-i-1))*(el2[i]+1)
    else:
        for i in range(N):
            res += (sps**(N-i-1))*((el2[i]+1)//2)
    return res


@cfunc(op_sig_32,
       locals=dict(n=float64, b=float64, sps=uint32), )
def op(op_struct_ptr, op_str, ind, N, args):
    op_struct = carray(op_struct_ptr, 1)[0]
    err = 0
    sps = args[0]
    n = int_to_list(op_struct.state, N, sps)[ind]
    b = sps**(N-ind-1)
    if sps == 3:
        if op_str == 43:  # "+" is integer value 43 = ord("+")
            if n == 1:
                op_struct.matrix_ele *= 0.
            else:
                op_struct.matrix_ele *= np.sqrt(2)
                op_struct.state += b
        elif op_str == 45:  # "-" is integer value 45 = ord("-")
            if n == -1:
                op_struct.matrix_ele *= 0.
            else:
                op_struct.matrix_ele *= np.sqrt(2)
                op_struct.state -= b
        elif op_str == 80:  # P
            if n == -1:
                op_struct.matrix_ele *= 2.
                op_struct.state += 2*b
            else:
                op_struct.matrix_ele *= 0.
        elif op_str == 77:  # M
            if n == 1:
                op_struct.matrix_ele *= 2.
                op_struct.state -= 2*b
            else:
                op_struct.matrix_ele *= 0.
        elif op_str == 122:  # "z" is integer value 122 ord("z")
            op_struct.matrix_ele *= n
        elif op_str == 90:  # Z (capital)
            op_struct.matrix_ele *= np.abs(n)
        else:
            op_struct.matrix_ele *= 0
            err = -1
    elif sps == 2:
        if op_str == 43:  # "+" is integer value 43 = ord("+")
            if n == 1:
                op_struct.matrix_ele *= 0.
            else:
                op_struct.state += b
        elif op_str == 45:  # "-" is integer value 45 = ord("-")
            if n == -1:
                op_struct.matrix_ele *= 0.
            else:
                op_struct.state -= b
        elif op_str == 105: # "i" is integer value 105 = ord("i")
            op_struct.matrix_ele *= 1
        elif op_str == 122:  # "z" is integer value 122 ord("z")
            op_struct.matrix_ele *= n
        else:
            op_struct.matrix_ele *= 0
            err = -1
    return err


@cfunc(next_state_sig_32,
       locals=dict(), )
def next_state(s, counter, N, args):
    return args[counter+1]


def get_s0_pcon(N, Np):
    dirname = Path('.')
    i_spin = fs.read_input()[1]["spin"]
    sps = int(2*i_spin+1)
    states_text = "spin-" + fs.spin_text[sps-2] + "_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
    path_to_load = dirname / "Data" / "States" / states_text
    states = np.load(path_to_load)
    return states[0]


def get_Ns_pcon(N, Np):
    dirname = Path('.')
    i_spin = fs.read_input()[1]["spin"]
    sps = int(2*i_spin+1)
    states_text = "spin-" + fs.spin_text[sps-2] + "_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
    path_to_load = dirname / "Data" / "States" / states_text
    states = np.load(path_to_load)
    return len(states)


@cfunc(map_sig_32,
       locals=dict())
def parity(x, N, sign_ptr, args):  # change + and -
    return P(x, N, args[0])


def corr(times, H, N, basis): #change to auto_correlator
    no_checks = dict(check_symm=False, check_pcon=False, check_herm=False)
    psi = np.random.rand(H.shape[0]) + 1j*np.random.rand(H.shape[0])
    psi /= np.linalg.norm(psi)
    psi_t = H.evolve(psi, times[0], times)
    Sz = hamiltonian([["z", [[1.0, N//2]]]], [],
                     dtype=np.float64, basis=basis, **no_checks)
    temp = Sz.dot(psi)
    temp_t = H.evolve(temp, times[0], times)
    C = [Sz.matrix_ele(psi_t.T[i], temp_t.T[i]) for i in range(len(times))]
    C = abs(np.array(C))
    return C

