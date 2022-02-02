# StackPy
*A collection of methods for data stacking*

![plot1](/figs/stackpy_logo.png)

## Introduction
A collection of methods for stacking of generic time series data.

## Available modules
This package is under active development. The currently available modules are listed here.

1. `utils`: This module contains frequently used utility functions.

2. `core`: This module contains core stacking functions.

## Installation
1. Install `stackpy` package functions using `pip`

```
$ pip install stackpy
```

This step will install the **StackPy** package. The modules would then be imported under any working directory.

2. Install with local copy:
```
$ conda create -n stackpy -c conda-forge jupyter numpy scipy pandas python
$ conda activate stackpy
```

`cd` to the directory you want to save the package files. Then run:

```
$ pip install .
```

1. Test the installation

Run the following commands to test your installation, under the root directory of StackPy.

```python
import os,pickle
import numpy as np
import matplotlib.pyplot as plt
from stackpy.core import stack
from scipy.signal import sosfiltfilt, butter

dataroot='./data'
dfile=dataroot+"/stackpy_testdataset.pk"
d=pickle.load(open(dfile,'rb'))

sos=butter(4,[0.05,0.5],fs=1/dt,btype="bandpass",output='sos')
plt.figure(figsize=(10,5))
scale=40
data,dt,lag,d_id=[d["data"],d["dt"],d['lag'],d['id']]
tx=np.arange(-lag,lag+0.5*dt,dt)
extent=[-lag,lag,data.shape[0],0]
dn=data.copy()

stack_method="robust"

for i in range(data.shape[0]):
    dn[i,:]=sosfiltfilt(sos,data[i,:]/np.max(np.abs(data[i,:])))

plt.imshow(dn,extent=extent,cmap="seismic",aspect="auto",alpha=0.7)

dstack=stack(dn,method=stack_method)
plt.plot(tx,scale*dstack+0.25*data.shape[0],'k',lw=2,label=stack_method)
plt.vlines(0,0,data.shape[0],'k')
plt.xlim([-200,200])
plt.ylim([0,data.shape[0]])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title(d_id)
plt.xlabel("time (s)",fontsize=14)
plt.ylabel("order",fontsize=14)
plt.legend(fontsize=12)
plt.show()

```

You should get the following plot:

![plot1](/figs/stack_example.png)


## Tutorials on key functionalities
See https://github.com/xtyangpsp/StackPy for tutorials and more detailed descriptions.


## Contribute
Any bugs and ideas are welcome. Please file an issue through GitHub https://github.com/xtyangpsp/StackPy.


## References
* To be added.
