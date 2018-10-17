import numpy as np
import random
import matplotlib.pyplot as plt

numClusters = 3
aLarge = 5
aSmall = 0.2
length = 100
numSamples = 100 #per cluster
refArray = ["A", "C", "T", "G"]
masterCluster = []
randNum = 50

#assign a low alpha to everything--for now we are using the same small alpha???? later random
for m in range(numClusters):
    s = np.random.dirichlet((aLarge,aLarge, aLarge, aLarge ), length)
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
            currentSeq.append(a)
            #can do refarray for previous (to converst to letters)
        #puts space cluster and goes to anew line in the file f
        masterCluster.append(currentSeq)

def MasterClusterPrint(cluster):
    for i in range(cluster.__len__()):
        for j in range(cluster[0].__len__()):
            print(cluster[i][j])
    print("")


#plotting:
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(6,10))

ax2.imshow(masterCluster, extent=[0,100,0,1], aspect='auto')
ax2.set_title('Auto-scaled Aspect')

plt.tight_layout()
plt.show()


#changing the alpha of each cluster to a higher number
print((random.randint(1,4))/10)



for m in range(numClusters):
    #for now aSmall is a fixed 0.1, but it can also be a random number between 0.1 and 0.4 like: ((random.randint(1,4))/10
    s = np.random.dirichlet((aSmall, aSmall, aSmall, aSmall), length)
    #the random sumple function takes in 2 args, the first one is a range (first num is included, second is not), second argument is the number of random nums you want from that range
    randomArr = (random.sample(range(0,100), randNum))
    print(randomArr)
    #print("s:", s)
    for i in range(randNum):
        for j in range(numSamples):
            rowNum = j + (100*m)
            pos = randomArr[i]
            randA = random.random()
            pointProb = s[pos]
            if randA < pointProb[0]:
                a = 0
            elif randA < (pointProb[0] + pointProb[1]):
                a = 1
            elif randA < (pointProb[0] + pointProb[1] + pointProb[2]):
                a = 2
            elif randA <= (pointProb[0] + pointProb[1] + pointProb[2] + pointProb[3]):
                a = 3
            masterCluster[rowNum][pos]= a

            #can do refarray for previous (to converst to letters)
        #puts space cluster and goes to anew line in the file f


#plots master array with new low alphas
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(6,10))

ax2.imshow(masterCluster, extent=[0,100,0,1], aspect='auto')
ax2.set_title('Auto-scaled Aspect')

plt.tight_layout()
plt.show()





#print masterCluster to a file
f= open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/mAlphaunshuffle.txt", "w+")
for m in range(numClusters):
    for n in range(numSamples):
        for l in range(length):
            print(refArray[masterCluster[m*100 + n][l]], end="")
            f.write(refArray[masterCluster[m*100 + n][l]])
        print(" ",  m, "")
        print(" ",  m, "", file = f)
        print("")

