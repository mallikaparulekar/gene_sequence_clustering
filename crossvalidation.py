import numpy as np
import math
import matplotlib.pyplot as plt
from sequenceObj import Sequence
from sklearn.metrics.cluster import adjusted_rand_score
# K = number of clusters
N = 100
# N = length of each sequence
S = 300
# S = number of data points
T = 270

kArr = [1,2,3,4,5,6,7,8,9,10]
#kArr = array of different clusters to try
weightProbMaster = []
startFromLine = T
#where to start from for the last part (S+1)

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



for i in range (kArr.__len__()):
    #yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/mALphashuffle.txt", "r")
    yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledwrittenFile.txt", "r")
#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledwrittenFile.txt", "r")
    finallist = []
    seqlist = []
    seqOBJArr = []
# seqlist is the list of all dna converted to a an array of 2 elements
    linNum = 0
    for line in yep:
        if (linNum < T):
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
        linNum+= 1
    yep.close()
    seqarray= np.asarray(seqlist)
    K = kArr [i]
    cluster_id = []
    # the following code creates the default cluster id array, with an equal number being in each cluster
    for i in range (T):
        c = (int) (i/(T/K))
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
        for i in range (T):
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
            for s in range (T):
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
    print(probability)
    yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/mALphashuffle.txt", "r")#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledwrittenFile.txt", "r")
    smallSeqOBJ= []
    linNum = 0
    for line in yep:
        #extracting the last 30 sequences
        if (linNum >= startFromLine):
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
            s = Sequence.make_sequence(int(values[1]), seq)
            smallSeqOBJ.append(s)
        linNum += 1
    yep.close()
    #empty array of length K (cluster number)
    trialFreq = [1]*K
    #the following loop calculates the trial frequencies of each cluser
    for i in range (seqOBJArr.__len__()):
        trialFreq[seqOBJArr[i].currentCluster]+= 1
    print(trialFreq)

    trialFreqsum = 0
    for i in range(trialFreq.__len__()):
        trialFreqsum += trialFreq[i]
    #find sum of trial frequencies

    #make trial Freq actuall a frequency by dividing it by sum
    for i in range(trialFreq.__len__()):
        trialFreq[i]= trialFreq[i]/trialFreqsum
    print("trialFreq", trialFreq)


    #calculate the probability of belonging to a cluster

    #create an array of K columns and length S-T
     #matrix = [[0 for i in range(a)] for i in range(b)]
    #creates a 2D array with b rows and a columns, set to O
    trialProb=[[0 for i in range(K)]for i in range(S-T)]
    for s in range (S-T):
        prob = -1000000000000000
        optimclust = -1
        for k in range (K):
            probcurrent= 0
            for n in range (N):
                 #seq = seqarray[s][n]
                 seq= seqOBJArr[s].value[n]
                 p = probability [4*k+seq][n]
                 probcurrent = probcurrent + math.log10(p)
            trialProb[s][k]= probcurrent
    print(trialProb)

    weightedProbability = [0]*(S-T)
    for l in range (trialProb.__len__()):
        wProb = 0
        for k in range (K):
            wProb += (10**(trialProb [l][k])) * trialFreq[k]
            #print("prob= ", 10**(trialProb [l][k]))
            #10 raised to to convert trial freq to probability (it is currently logs)
        weightedProbability[l]= wProb
    weightProbMaster.append(weightedProbability)


print("wPM", weightProbMaster)

sumArr = [0]*(weightProbMaster.__len__())

for r in range(weightProbMaster.__len__()):
    sum = 0
    for c in range (S-T):
        sum = sum + weightProbMaster[r][c]
    sumArr[r]= sum
    


print("sumArr", sumArr)

'''
#following code only works for when you are comparing 2 arrays
betterArr = [0]*(weightProbMaster.__len__())
#checks for each sequence which clustering was better
for c in range(S-T):
    weights = []
    for i in range(weightProbMaster.__len__()):
        weights.append(weightProbMaster[i][c])

    if (weights[0]>weights[1]):
        betterArr[0]+=1
    if (weights[1]>weights[0]):
        betterArr[1]+=1


print(betterArr)
'''
#sorts the array copy- using bubble sort ayyy
#in descending order--best to worst
sumArrCopy = sumArr.copy()
for j in range(sumArrCopy.__len__()):
    for i in range(sumArrCopy.__len__()-j-1):
        if (sumArrCopy[i]<sumArrCopy[i+1]):
            temp = sumArrCopy[i]
            sumArrCopy[i]= sumArrCopy[i+1]
            sumArrCopy[i+1]= temp

print("copy", sumArrCopy)

arrayRanking = [1]*(sumArrCopy.__len__())
#set to one from the start to since we are testing cluster numbers from 1-10, not 0-9
for i in range(sumArrCopy.__len__()):
    for j in range(sumArr.__len__()):
        if (sumArrCopy[i]==sumArr[j]):
            arrayRanking[i] += j

print(arrayRanking)
print("sumArr", sumArr)

