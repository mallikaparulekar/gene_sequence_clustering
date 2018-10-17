import numpy as np
import random

a = np.array([1, 2, 3])
a = np.array([1, 2, 3])
import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

# alpha distribution
numClusters = 3
x = 0.01
length = 100
numSamples = 100
refArray = ["A", "C", "T", "G"]
masterCluster = []

# num  = total number of sequences


# s = np.array([[0.250226366, 0.0002297601832, 0.165059640, 0.584416391],
# [0.0000194563325, 0.000000000245599067, 0.337103929, 0.662876614]])
# s = np.array([[2.50226366e-01, 2.97601832e-04, 1.65059640e-01, 5.84416391e-01],
# [1.94563325e-05, 2.45599067e-10, 3.37103929e-01, 6.62876614e-01]])


# printing a and b prob. of first cluster as one array

# setting a and b prob. of first cluster


# printing aI prob
#clearing the written file and the shuffle file
f = open('/Users/mallika/PycharmProjects/DirichletBio/venv/lib/writtenFile.txt', 'r+')
f.truncate(0)
f.close()

f = open('/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledwrittenFile.txt', 'r+')
f.truncate(0)
f.close()

f = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/writtenFile.txt", "w+")
for m in range(numClusters):
    s = np.random.dirichlet((x, x, x, x), length)
    #print("s:", s)
    allSequences = []
    for i in range(numSamples):
        currentSeq = []
        for j in range(length):
            # generating the random data
            # print("")
            randA = random.random()
            # print("apoint rand is: ", randA)
            # should it be less than equal to?
            # sometimes data is printed in decimal sometimes with e-why?
            pointProb = s[j]
            if randA < pointProb[0]:
                a = 0
            elif randA < (pointProb[0] + pointProb[1]):
                a = 1
            elif randA < (pointProb[0] + pointProb[1] + pointProb[2]):
                a = 2
            elif randA <= (pointProb[0] + pointProb[1] + pointProb[2] + pointProb[3]):
                a = 3
            currentSeq.append(refArray[a])
            print(refArray[a], end="")
            #writes each letter of the sequence
            f.write(refArray[a])

        print(" ",  m, "")
        #puts space cluster and goes to anew line in the file f
        print(" ",  m, "", file = f)
        allSequences.append(currentSeq)
    masterCluster.append(allSequences)
    print("al;lseq :", allSequences.__len__())



'''
#steps to refresh simulation
- run refined gen
- run shuffler
- run plotting 2.0
'''
