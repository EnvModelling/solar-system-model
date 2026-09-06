import matplotlib.pyplot as plt
import numpy as np

if __name__=="__main__":
	t=np.linspace(0,1000,1000)
	sin1=np.sin(2*np.pi*t/200.)
	cos2=0.1*np.cos(2*np.pi*t/100.)
	sin2=0.1*np.sin(2*np.pi*t/(200/3))
	cos3=0.01*np.cos(2*np.pi*t/(200/4))
	plt.ion()
	plt.plot(t,sin1)
	plt.plot(t,cos2)
	plt.plot(t,sin2)
	plt.plot(t,cos3)
	plt.plot(t,sin1+cos2+sin2+cos3,lw=3)
	plt.show()
	