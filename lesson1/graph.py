import sys
import matplotlib.pyplot as plt

x=[]
y=[]
with open(sys.argv[1]) as f:
    for line in f:
      try:
	x.append( float(line.split()[0]))
	y.append( float(line.split()[1]))
      except Exception:
	pass
      
plt.plot(x, y)
plt.show()
