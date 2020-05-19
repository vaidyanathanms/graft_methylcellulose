# To copy the resultant files for further analysis

import numpy
import os
import shutil
import subprocess
import sys
import glob

iarr = [1,2,3,4,5,6]
flag1 = 1
flag2 = 0
val1 = 1
val2 = 1
for ival in range(len(iarr)):
    print(iarr[ival])

    if flag1 == 0:
        continue
    else:
        val1 = val1 + 1

    if flag2 == 0:
        break
    else:
        val2 = val2 + 1
        
print(val1)
print(val2)
