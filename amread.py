import matplotlib.pyplot as plt
import numpy as np


def reader(fil,ylabel, ylabel2):
    f=open(str(fil),"r")
    lines=f.readlines()
    result0=[]
    result1=[]
    result2=[]

    for x in lines:
        result0.append(x.split(' ')[0])
        result1.append(x.split(' ')[1])
        result2.append(x.split(' ')[2])

    f.close()

    r0 = np.array([result0]).astype(np.float)
    r1 = np.array([result1]).astype(np.float)
    r2 = np.array([result2]).astype(np.float)

    r0 = r0[0]
    r1 = r1[0]
    r2 = r2[0]

#    plt.figure(1, figsize = (10,5))
#   plt.plot(r0,r1)
#   plt.xlabel('Frequency[GHz]')
#   plt.ylabel(str(ylabel))

#   plt.figure(2, figsize = (10,5))
#   plt.plot(r0,r2)

#   plt.xlabel('Frequency[GHz]')
#   plt.ylabel(str(ylabel2))

#plt.show()
    return r0, r1, r2


print '''Usage: Calling "python -i amread.py" in the command line: input a file name as a string with def reader(path/to/fil,ylabel, ylabel2) to get the transmission spectra and brightness temperature plots coming from an am atmospheric model. It returns the 3 column data that would make up the plots
    
    Example: 
    >>> a,b,c = reader('example2.1.out', 'Transmittance', '$T_b ~[K]$')
    '''

