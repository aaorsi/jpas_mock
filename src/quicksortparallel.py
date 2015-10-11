import random, time, sys
from multiprocessing import Process, Pipe

def quicksortParallel(lyst, conn, procNum):
    """
    Partition the list, then quicksort the left and right
    sides in parallel.
    """

    if procNum <= 0 or len(lyst) <= 1:
        #In the case of len(lyst) <= 1, quicksort will
        #immediately return anyway.
        conn.send(quicksort(lyst))
        conn.close()
        return

    #Create two independent lists (independent in that
    #elements will never need be compared between lists).
    pivot = lyst.pop(random.randint(0, len(lyst)-1))

    leftSide = [x for x in lyst if x < pivot]
    rightSide = [x for x in lyst if x >= pivot]

    #Creat a Pipe to communicate with the left subprocess
    pconnLeft, cconnLeft = Pipe()
    #Create a leftProc that executes quicksortParallel on
    #the left half-list.
    leftProc = Process(target=quicksortParallel, \
                       args=(leftSide, cconnLeft, procNum - 1))
    
    #Again, for the right.
    pconnRight, cconnRight = Pipe()
    rightProc = Process(target=quicksortParallel, \
                       args=(rightSide, cconnRight, procNum - 1))

    #Start the two subprocesses.
    leftProc.start()
    rightProc.start()

    #Our answer is the concatenation of the subprocesses' 
    #answers, with the pivot in between. 
    conn.send(pconnLeft.recv() + [pivot] + pconnRight.recv())
    conn.close()

    #Join our subprocesses.
    leftProc.join()
    rightProc.join()

if __name__ == '__main__':
	quicksortParallel(argv)	
