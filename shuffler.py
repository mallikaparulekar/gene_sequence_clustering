import random
#this prints out the shuffled sequence
#yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/mAlphaunshuffle.txt", "r")
#f = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/mAlphashuffle.txt", "w+")
yep = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/writtenFile.txt", "r")
f = open("/Users/mallika/PycharmProjects/DirichletBio/venv/lib/shuffledWrittenFile.txt", "w+")
S = 300
#S = num of clusters
seqList = []


for i in range(S):
    line= yep.readline()
    seqList.append(line)
yep.close()


random.shuffle(seqList)
for i in range(S):
    print(seqList[i], end="")
    print(seqList[i], end="",  file = f)
