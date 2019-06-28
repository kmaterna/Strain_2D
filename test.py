import numpy as np
import sys
import subprocess
import scipy.io.netcdf as netcdf
import matplotlib.pyplot as plt

print('Hello');
print("Exiting");
# sys.exit(0);
print("Haven't exited yet");
subprocess.call(['mkdir','-p','hooray'],shell=False)
print("Made directory!")
# subprocess.call(['pwd'], shell=False, cwd="Strain_Code")
subprocess.call(['pwd'], shell=False, cwd="Strain_Code")
print("printed directory")