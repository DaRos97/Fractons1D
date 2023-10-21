import sys
import numpy as np
from pathlib import Path
from colorama import Fore, Style # ,Back
                    #Fore: BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE, RESET.
                    #Back: BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE, RESET.
                    #Style: DIM, NORMAL, BRIGHT, RESET_ALL
import functions as fs

#### Input spin, N and Np
sps,N,Np = fs.input_routine()
#### Calculate states
states = fs.states(sps,N,Np[0],Np[1])
tot = len(states)
#### Save
dirname = Path('.')
states_text = "spin-" + fs.spin_text[sps-2] + "_N=%i_(q,p)=(%i,%i).npy" % (N,Np[0],Np[1])
path_to_save = dirname / "Data" / "States" / states_text
np.save(path_to_save,states)
print(Fore.GREEN)
print(tot, ' states saved in file ', path_to_save)
print(Style.RESET_ALL)
