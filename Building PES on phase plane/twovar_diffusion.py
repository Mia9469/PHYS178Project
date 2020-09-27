import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

a = 270
b = 108
d = 0.154
gamma = 0.641
taus = 100*10**-3
tauampa = 2*10**-3
Jn11, Jn12 = 0.2609, 0.0497
Jext = 5.2*(10**-4)
I0 = 0.3255
mu0 = 30
dt = 10**-4
cprime = 0
stdnoise = 0.02
D=stdnoise**2/tauampa*2

starttime = -0.5
endtime = 2.5
steps = int(abs(starttime - endtime)/dt)
time = np.linspace(starttime, endtime, steps)

def H(xi):
	return (a*xi-b)/(1-np.exp(-d*(a*xi-b)))

def dH(xi):
	return (a*(1-np.exp(-d*(a*xi-b)))-(a*xi-b)*np.exp(-d*(a*xi-b))*d*a)/(1-np.exp(-d*(a*xi-b)))**2

def experiment():
	H1, H2 = np.zeros(steps+1), np.zeros(steps+1)
	S1, S2 = np.random.rand(steps+1), np.random.rand(steps+1)
	#S1, S2 = np.zeros(steps+1), np.zeros(steps+1)
	sig = np.zeros((2,2))
	for index, t in enumerate(time):
		if t > 0:
			x1 = Jn11*S1[index] - Jn12*S2[index] + I0 + Jext*mu0*(1+cprime/100)
			x2 = Jn11*S2[index] - Jn12*S1[index] + I0 + Jext*mu0*(1-cprime/100)
		else:
			x1 = Jn11*S1[index] - Jn12*S2[index] + I0 
			x2 = Jn11*S2[index] - Jn12*S1[index] + I0 
		H1[index+1] = H(x1)
		H2[index+1] = H(x2)
		S1[index+1] = S1[index] + dt*(-S1[index]/taus + (1 - S1[index])*gamma*H1[index])
		S2[index+1] = S2[index] + dt*(-S2[index]/taus + (1 - S2[index])*gamma*H2[index])
		f11 = -1/taus - gamma*H1[index+1] + (1 - S1[index+1])*gamma*Jn11*dH(x1)
		f22 = -1/taus - gamma*H2[index+1] + (1 - S2[index+1])*gamma*Jn11*dH(x2)
		f12 = (S1[index+1] - 1)*gamma*Jn12*dH(x1)
		f21 = (S2[index+1] - 1)*gamma*Jn12*dH(x2)
		A = [[f11,f12],[f21,f22]]
		sig = sig+dt*(sig@np.transpose(A)+A@sig+[D,D])
	return H1[1:], H2[1:], S1[1:], S2[1:],sig

def slided(data):
	timestep = 5*10**-3
	slided_data = []
	for index, value in enumerate(data):
		if index % int(timestep/dt) == 0:
			slided_data.append(value)
	return slided_data


def smoothing(data):
	length = len(data)
	smoothed_data = np.zeros(length)
	width = int(10*10**-3/dt)
	for i in range(length):
		if length - (i+1) < width:
			smoothed_data[i] = np.average(data[i:])
		else: smoothed_data[i] = np.average(data[i: i+width])
	return smoothed_data


def nullclines(cprime,mu):

	def gate(S):
		# Input Currents 
		# Zero on the end is for easy moving to zero input
		# I[2]
		I = [Jext*mu*(1+(cprime/100)), 
			 Jext*mu*(1-(cprime/100)),0]
		# x variable
		x = [Jn11*S[0] - Jn12*S[1] + I0 + I[0],
			 Jn11*S[1] - Jn12*S[0] + I0 + I[1]]

		# H values	 
		H = [(a*x[i] - b)/(1 - np.exp(-d*(a*x[i] - b))) for i in range(2)] 

		# dS1/dt and dS2/dt
		expr = [-(S[i]/taus) + (1-S[i])*gamma*H[i] for i in range(2)]
		return expr

	# Brute Force approach

	iterr = np.arange(0,1,0.001)
	ds1 = []
	ds2 = []
	dr1 = []
	dr2 = []
	fix_points = []
	fix_points_r = []

	for i in iterr:
	
		for j in iterr:
			y,x = gate([i,j])
			I = [Jext*mu*(1+(cprime/100)), 
			 	Jext*mu*(1-(cprime/100)),0]
			xx = [Jn11*i - Jn12*j + I0 + I[0], Jn11*j - Jn12*i + I0 + I[1]]	 
				
			if 0.01 > x > -0.01:
				ds1 += [[i,j]]
				dr1 += [[(a*xx[k] - b)/(1 - np.exp(-d*(a*xx[k] - b))) for k in range(2)]]
			if 0.01 > y > -0.01:
				ds2 += [[i,j]]
				dr2 += [[(a*xx[k] - b)/(1 - np.exp(-d*(a*xx[k] - b))) for k in range(2)]]
			if 0.01 > x > -0.01 and 0.01 > y > -0.01:
				fix_points += [[i,j]]
				fix_points_r += [[(a*xx[k] - b)/(1 - np.exp(-d*(a*xx[k] - b))) for k in range(2)]]

	ds1=np.vstack(ds1)
	ds2=np.vstack(ds2)
	fix_points=np.vstack(fix_points)
	dr1=np.vstack(dr1)
	dr2=np.vstack(dr2)
	fix_points_r=np.vstack(fix_points_r)
	return ds1,ds2,fix_points,dr1,dr2,fix_points_r


plt.figure
for i in range(1):
	result = experiment()
	hue = ['red','blue']
	plt.plot(time*1000, smoothing(result[0]), color = hue[0])
	plt.plot(time*1000, smoothing(result[1]), color = hue[1])
plt.plot(time*1000, 15*np.ones(steps),color = 'black')
plt.xlabel('Time(ms)')
plt.ylabel('Firing rate/r(Hz)')
plt.ylim(top = 20)
plt.text(-125, 16, 'Threshold')
plt.show()

for i in range(8):
	cprime=6.4*(i+1)
	n_s=nullclines(cprime,mu0)
	sio.savemat('fixpoints_cp'+str(i+1)+'.mat', mdict={'fix_'+str(i+1):n_s[2]})

mu0 = 30
n_s=nullclines(cprime,mu0)
result = experiment()
# plot S
plt.figure
plt.scatter(n_s[0][:,0],n_s[0][:,1],s=1)
plt.scatter(n_s[1][:,0],n_s[1][:,1],s=1)
plt.scatter(n_s[2][:,0],n_s[2][:,1])

for i in range(1):
	plt.scatter(result[2],result[3],s=1)

plt.xlabel('S2')
plt.ylabel('S1')
plt.legend(['dS2/dt=0','dS1/dt=0','Fixed Points'])
plt.title('Nullclines cprime: {}'.format(cprime))
plt.show()
result


mu0=30

cprime=6.4*8
result = experiment()
result[4]


