import numpy as np
import matplotlib.pyplot as plt
import os

corrFilePath = "C:\\Users\\Junliu\\Desktop\\160Line\\Z-Z\\430L160"
corrFileNames = os.listdir(path)

for corrFileName in corrFileNames:
    rawCorrData = np.loadtxt(corrFileName)
    