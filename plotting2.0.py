# revised for clusters of > 2 length
import numpy as np
import math
import matplotlib.pyplot as plt
from sequenceObj import Sequence
from sklearn.metrics.cluster import adjusted_rand_score
K = 5
# K = number of clusters
N = 100
# N = length of each sequence
S = 1000
# S = number of data points

# A = 0
# C = 1
# T = 2
# G = 3

# data = [[A,A],[C,C],[T,T],[G,G],[A,A],[C,C],[T,T]]
# print(data)
# data_full = [[G,A]]
# size = 5
# for i in range(size):
# data_full = np.concatenate((data_full, data))

# print(data_full)
#this code segments reads ATCG values from a flie, and converts them to arrays of 1,2,3,4
#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledExtendedseq.txt", "r")
#writtenText is unshuffled
yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/leeluMtest.txt", "r")
#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledwrittenFile.txt", "r")
finallist = []
seqlist = []
seqOBJArr = []
# seqlist is the list of all dna converted to a an array of 2 elements
for line in yep:
    values = line.split()
    seq = []
    #line = yep.readline()
    for j in range(N):
        if values[0][j] == "A":
            seq.append(0)
        elif values[0][j] == "C":
            seq.append(1)
        elif values[0][j] == "T":
            seq.append(2)
        else:
            seq.append(3)
    seqlist.append(seq)
    s = Sequence.make_sequence(int(values[1]), seq)
    seqOBJArr.append(s)
yep.close()

seqarray = np.asarray(seqlist)


cluster_id = []
# the following code creates the default cluster id array, with an equal number being in each cluster
for i in range (S):
    c = (int) (i/(S/K))
    cluster_id.append(c)
    #assign current cluster to the objects
    seqOBJArr[i].currentCluster = c


#matrix = [[]]
#matrix = [[0 for i in range(a)] for i in range(b)]
#creates a 2D array with b rows and a columns, set to O

counter = [[]]
counter = [[1 for i in range(N)] for i in range(4*K)]

probability = [[]]
probability =  [[0 for i in range(N)] for i in range(4*K)]

def counterCalc (seqOBJArr, counter):
#iterating through the sequence array and entering the number of ACTG for each cluster. For example, in a scenario with four clusters, the table would have 4 times 4 = 16 rows and 2 columns (N)
    #resets all the rows and columns of the counter to one
    for row in range(4*K):
        for column in range (N):
            counter[row][column]= 1
    #recalculates the counts of each ACTG of each cluster
    for i in range (S):
        k = seqOBJArr[i].currentCluster
        for j in range (N):
            m = seqOBJArr[i].value[j]
            counter [4*k + m][j] = counter [4*k + m][j] + 1

def probabilityCalc (probability, counter):
#goes 1 column, 1 cluster at a time and uses the counter array and a sum function calculate the probability of each ACTG
    for k in range (K):
        for n in range (N):
            summ = 0
            for b in range (4):
                summ = summ + counter [4*k + b][n]
            for b in range (4):
                probability[4*k+b][n] = counter [4*k + b][n]/summ

def ClusterCalc(seqOBJArr, probability, counter):
    changeCount = 10
    while (changeCount != 0):
        changeCount = 0
        #counts the number of times a sequence's cluster is reassigned
        counterCalc(seqOBJArr, counter)
        probabilityCalc (probability, counter)
        # print("cluster: ",cluster_id)
        # print("probability: ", probability)
        for s in range (S):
            prob = -1000000000000000
            optimclust = -1
            for k in range (K):
                probcurrent= 0
                for n in range (N):
                    seq = seqarray[s][n]
                    p = probability [4*k+seq][n]
                    probcurrent = probcurrent + math.log10(p)
                if (probcurrent > prob):
                    #print("probcurrent: ",probcurrent," prob: ", prob)
                    prob = probcurrent
                    optimclust = k

            old_cluster_id1 = seqOBJArr[s].currentCluster
            #print("old :", old_cluster_id1)
            seqOBJArr[s].currentCluster= optimclust
            #print("optim :", optimclust)
            #old_cluster_id = cluster_id[s]
            #print("old2 :", old_cluster_id1)
            #cluster_id[s]= optimclust
            #print("optim2 :", optimclust)
            #print(cluster_id[s])
            #print(optimclust, old_cluster_id1)
            if (optimclust !=  old_cluster_id1):
                changeCount = changeCount + 1
                #print("change:", changeCount)
        print("changeCount: ", changeCount)
        #checking if there were no swaps made:

ClusterCalc(seqOBJArr, probability, counter)


sortedseqOBJList = []
#sorted by new cluster ids
for k in range(K):
    for s in range(S):
        if(seqOBJArr[s].currentCluster== k):
            sortedseqOBJList.append(seqOBJArr[s])
print(sortedseqOBJList)

def checkClustering(seqOBJArr):
    totalCorrect = 0
    for i in range(seqOBJArr.__len__()):
        print("cc: ", seqOBJArr[i].currentCluster, "tc: ", seqOBJArr[i].trueCluster, totalCorrect)
        if(seqOBJArr[i].currentCluster == seqOBJArr[i].trueCluster):
            totalCorrect += 1
    return totalCorrect

print("clusterCheck: ",checkClustering(seqOBJArr))


'''
#plot of unsorted list
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(6,10))

ax2.imshow(seqarray, extent=[0,100,0,1], aspect='auto')
ax2.set_title('Auto-scaled Aspect')
plt.tight_layout()
plt.show()
'''
#plot sorted list

#checking the clustering using ARI
TrueClusters = []
for i in range (seqOBJArr.__len__()):
    TrueClusters.append(seqOBJArr[i].trueCluster)

AssignedClusters = []
for i in range (seqOBJArr.__len__()):
    AssignedClusters.append(seqOBJArr[i].currentCluster)

print("True Clusters: ", TrueClusters)
print("Assgined Clusters: ", AssignedClusters)

#check ARI
print(adjusted_rand_score(TrueClusters, AssignedClusters))

'''
alphaArr = [0.01, 0.1,0.2,0.5,0.75,1.0,2.0]
ariArr = [0.56206, 0.96054, 0.96053, 0.82823, 0.33733, 0.211045, 0.07455]
plt.plot(alphaArr, ariArr)
plt.show()
'''
#xhecking with leelu m actual
ans = []
q = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/actual.txt", "r")

for line in q:
    values=line.split()
    if (values[0]=="2"):
        ans.append(2)
    elif (values[0]=="1"):
        ans.append(1)
    elif (values[0]=="3"):
        ans.append(3)
    elif (values[0]=="4"):
        ans.append(4)
    elif (values[0]=="5"):
        ans.append(5)

print(ans)

print(adjusted_rand_score(ans, AssignedClusters))
