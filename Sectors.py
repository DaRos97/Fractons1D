from pathlib import Path
import numpy as np
import sys
dirname = Path('.')
import functions as fs
import matplotlib.pyplot as plt
from tqdm import tqdm

sps,N,Np = fs.input_routine()

maxQ = N
if N%2:
    maxP = N//2*(N//2+1)
else:
    maxP = N//2*N//2

sectors = np.zeros((2*maxQ+1,2*maxP+1),dtype=np.uint32)
for n in tqdm(range(sps**N)):
    sectors[maxQ+fs.totQ(n,N,sps)][maxP+fs.totP(n,N,sps)] += 1

max_index = np.unravel_index(sectors.argmax(), sectors.shape)
print("Index of the max value: ",max_index[0]-maxQ,max_index[1]-maxP)
print("Highest dimension: ",sectors[max_index[0]][max_index[1]])
#### Plot
q_list = np.linspace(-N,N,2*N+1,dtype=int)

for q in q_list:
    print(q)
    plt.plot(range(-maxP,maxP+1),sectors[maxQ+q],label='$q=$%i' % q)

plt.plot(max_index[1]-maxP,sectors[max_index[0]][max_index[1]],"k*")

plt.xlabel('P')
plt.ylabel('$D_{(q,p)}$')

plt.legend()
plt.show()
