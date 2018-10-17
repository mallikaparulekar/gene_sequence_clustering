# importing the required module
import matplotlib.pyplot as plt

#correctness graph- how many times for each alpha (and actual k = 3) did algo work- check google doc
# x axis values
x = [0.01,0.1,0.2, 0.3, 0.4]
# corresponding y axis values
y = [0.7,0.9, 0.9,0.1, 0.1]

# plotting the points
plt.plot(x, y)

# naming the x axis
plt.xlabel('alpha')
# naming the y axis
plt.ylabel('# of correct classification')

# giving a title to my graph
plt.title('Frequency of accuracy of clustering (K= 3)')

# function to show the plot
plt.show()


#score graph- if 3 is the correct cluster: 100, if it is second: 90, third: 80............last : 10 points--plot average points
# x axis values
x = [0.01,0.1,0.2, 0.3, 0.4]
# corresponding y axis values
y = [80,99,99,86,72]

# plotting the points
plt.plot(x, y)

# naming the x axis
plt.xlabel('alpha')
# naming the y axis
plt.ylabel('classification score')

# giving a title to my graph
plt.title('Average clustering score when (K= 3)')

# function to show the plot
plt.show()
