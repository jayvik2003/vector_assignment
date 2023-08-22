import numpy as np
import matplotlib.pyplot as plt
#If using termux
import subprocess
import shlex
y = np.random.randint(-6,6, size=(1, 2))

A = np.random.randint(-6,6, size=(1, 2))
B = np.random.randint(-6,6, size=(1, 2))
C = np.random.randint(-6,6, size=(1, 2))

d = B- A
e = C - B
f = A - C

print("The direction vector of AB is ",d)
print("The direction vector of BC is ",e)
print("The direction vector of CA is ",f)


