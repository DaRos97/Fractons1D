import numpy as np
import sys
from colorama import Fore, Style  # ,Back
from tqdm import tqdm


#use file input = 1
#spin = 1
#length = 11
#use highest subsector = 1
#charge = 1
#dipole moment = -1


spin_text = ["0,5", "1"]
spin_ranges = [[4,5],[3,4]]


def read_input():
    error = 0
    f = open("input.txt",'r')
    inputs = dict()
    for line in f.readlines():
        sline = line.split(' = ')
        if sline[0]=="spin" and int(sline[1])==0:
            inputs[sline[0]] = 0.5
        else:
            inputs[sline[0]] = int(sline[1])
    # Define error checking on input
    #if len(inputs) != 5:
    #    error = 1
    f.close()
    return error, inputs

def check_in(sps, N, Np):
    bounds = [2, 3, 3, 20, -N, N, -N*(N+1)//2, N*(N+1)//2]      #check bounds for p value when computing from the center of the chain
    error = 0
    if sps < bounds[0] or sps > bounds[1]:
        print(Fore.RED + 'Error in spin: must be either 0 (spin-1/2) or 1 (spin-1)')
        error = 1
    if N < bounds[2] or N > bounds[3]:
        print(Fore.RED + 'Error in length: must be between ',bounds[2], ' and ',bounds[3])
        error = 1
    if Np[0] < bounds[4] or Np[0] > bounds[5]:
        print(Fore.RED + 'Error in charge (third input): for length = ',
              N, ' it must be between ', bounds[4], ' and ', bounds[5])
        error = 1
    if Np[1] < bounds[6] or Np[1] > bounds[7]:
        print(Fore.RED + 'Error in dipole (fourth input): for length = ',
              N, ' it must be between ', bounds[6], ' and ', bounds[7])
        error = 1
    if error:
        print(Style.RESET_ALL)
        exit()


def input_routine():
    err, inputs = read_input()
    if inputs["use file input"]:
        if err:
            print(Fore.RED + "Detected only ", len(inputs)-1,
                  " inputs from file, but 4 are required: ")
            print("- spin (0=spin-1/2, 1=spin-1)")
            print("- chain length (between 2 and 20")
            print("- charge")
            print("- dipole moment")
            print("Please retry.")
            print(Style.RESET_ALL)
            exit()
        sps = int(2*inputs["spin"] + 1)
        N = inputs["length"]
        if inputs["use highest subsector"]:
            Np = maxSector(sps,N)
            check_in(sps, N, Np)
            rewrite_input(Np)
        else:
            Np = (inputs["charge"], inputs["dipole moment"])
            check_in(sps, N, Np)
        print("Using input.txt file as input: N=",N,", sps=",sps,", q=",Np[0],", p=",Np[1])
    else:
        print("No suitable file input.txt detected, Retry.")
        exit()
    return sps, N, Np

def int_to_list(n,N,sps):
    el = np.array(list(np.base_repr(n, sps)), dtype=np.int32())
    el = np.append(np.zeros(N-len(el), dtype=np.int32()), el)
    if sps == 3:
        el -= 1
    elif sps == 2:
        el *= 2
        el -= 1
    else:
        print(Fore.RED + "Invalid sps value")
        print(Style.RESET_ALL)
        exit()
    return el

def totQ(n, N, sps):
    el = int_to_list(n, N, sps)
    q = np.sum(el)
    return int(q)


def totP(n, N, sps):
    el = int_to_list(n, N, sps)
    p = 0
    for i in range(N):
        p += (i-N//2)*(el[i])
    return int(p)

def states(sps,N,q,p):
    sts = np.ndarray(0,dtype=np.uint32)
    for i in tqdm(range(sps**N)):
        if totQ(i,N,sps) == q and totP(i,N,sps) == p:
            sts = np.append(sts,i)
    return sts

def maxSector(sps,N):
    maxQ = N
    if N%2:
        maxP = N//2*(N//2+1)
    else:
        maxP = N//2*N//2
    sectors = np.zeros((2*maxQ+1,2*maxP+1),dtype=np.uint32)
    for n in tqdm(range(sps**N)):
        sectors[maxQ+totQ(n,N,sps)][maxP+totP(n,N,sps)] += 1
    max_index = np.unravel_index(sectors.argmax(), sectors.shape)
    Np = (max_index[0]-maxQ,max_index[1]-maxP)
    return Np


def rewrite_input(Np):
    datas = (0, Np[0], Np[1])
    f = open("input.txt")
    string_list = f.readlines()
    f.close()
    for i in range(3,6):
        sline = string_list[i].split(' = ')
        sline[1] = datas[i-3]
        string_list[i] = sline[0] + ' = ' + str(sline[1]) + '\n'
    f = open("input.txt",'w')
    new_content = "".join(string_list)
    f.write(new_content)
    f.close()


def interaction_list(N,interaction_length):
    J = []
    for n in range(N-interaction_length+1):
        temp = [-1.]
        temp += [i for i in range(n,n+interaction_length)]
        J.append(temp)
    return J
