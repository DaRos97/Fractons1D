import numpy as np

#We start by computing, for given chain length N and spin value S, all the (2S+1)**N possible states and classify them by charge and momentum
#Charge is defined by teh total magnetization: Q=\sum_i Z_i
#Momentum is defined wrt the chain center -> easy for N odd, for N even take right one (index N//2 in both cases)
#The definition of momentum will be then P=\sum_i i*Z_i

#Higher momenta can be defined through P_n=\sum_i i**n*Z_i, but not for now

N = 10
S = 0.5
S_txt = {0.5:"1/2", 1:"1", 1.5: "3/2", 2:"2"}
S_i = int(2*S+1)    #number of possible configurations of each site
tot_states = S_i**N

print("Computing states for chain of length ",N," for spin ",S_txt[S])

print("Total number of states: ",tot_states)

#Index states -> 1-1 correspondence b/w index and states
def state(index_,S_i_,N_):
    state = np.zeros(N)
    state[:len(np.base_repr(index_,base=S_i_))] = list(np.base_repr(index_,base=S_i_))
    return state

for i in range(tot_states):
    if i % 123 == 0:
        print("State associated to ",i," is ",state(i,S_i,N))
        input()


