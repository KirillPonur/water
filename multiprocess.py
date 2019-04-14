from multiprocessing import Pool
from pylab import *
import os
import multiprocessing as mp
from numpy import random
import ctypes
 


def water(x):
    model=0
    N=10
    for m in range(N):
        model+=cos(x)
    proc_name = current_process().name
    print('{0} by: {1}'.format(model, proc_name))
    return model

def water_worker(x,queue):
    nP = len(x)
    N=10
    D = np.zeros(len(x))

    while True:
        job = queue.get()
        if job == None:
            break

        start = int(job[0])
        stop = int(job[0] + job[1])

        D[start:stop] = np.sum([cos(x) for m in range(N)])
        queue.task_done()
    queue.task_done()

def mpCalc(x):
    # allocate shared array
    nP = len(x)
    # setup jobs
    nCPU = mp.cpu_count()
    # nCPU = 4
    nJobs = nCPU * 1
    q = nP / nJobs
    r = nP % nJobs
    jobs = []
    firstRow = 0
    for i in range(nJobs):
        rowsInJob = q
        if (r > 0):
            rowsInJob += 1
            r -= 1
        jobs.append((firstRow, rowsInJob))
        firstRow += rowsInJob

    queue = mp.JoinableQueue()
    
    for job in jobs:
        queue.put(job)
    for i in range(nCPU):
        queue.put(None)
    
    # run workers
    workers = []
    for i in range(nCPU):
        worker = mp.Process(target = water_worker,
                            args = (x,queue))
        workers.append(worker)
        worker.start()
        queue.join()
   
    # # make array from shared memory    
    # D = np.reshape(np.frombuffer(arrD), (nP, nQ))
    return D

def main():
    N=100
    psi=array([
    [ random.uniform(0,2*pi) for j in range(N)] 
      for i in range(N)              ])
    x=linspace(0,10,N)
    print(mpCalc(x))

if __name__ == '__main__':
    main()