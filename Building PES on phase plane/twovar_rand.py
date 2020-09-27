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
stdnoise = 0.02
mu0 = 30
dt = 0.1*10**-3
cprime = 0

def H(xi):
	return (a*xi-b)/(1-np.exp(-d*(a*xi-b)))
starttime = -0.5
endtime = 1
steps = int(abs(starttime - endtime)/dt)
time = np.linspace(starttime, endtime, steps)

def experiment():
	H1, H2, S1, S2 = np.zeros(steps+1), np.zeros(steps+1), np.random.rand(steps+1),np.random.rand(steps+1)
	Inoise1, Inoise2 = np.zeros(steps+1), np.zeros(steps+1)
	for index, t in enumerate(time):
		Inoise1[index+1] = Inoise1[index] + dt*(-Inoise1[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		Inoise2[index+1] = Inoise2[index] + dt*(-Inoise2[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		if t > 0:
			x1 = Jn11*S1[index] - Jn12*S2[index] + I0 + Jext*mu0*(1+cprime/100) + Inoise1[index]
			x2 = Jn11*S2[index] - Jn12*S1[index] + I0 + Jext*mu0*(1-cprime/100) + Inoise2[index]
		else:
			x1 = Jn11*S1[index] - Jn12*S2[index] + I0 + Inoise1[index]
			x2 = Jn11*S2[index] - Jn12*S1[index] + I0 + Inoise2[index]
		H1[index+1] = H(x1)
		H2[index+1] = H(x2)
		S1[index+1] = S1[index] + dt*(-S1[index]/taus + (1 - S1[index])*gamma*H1[index])
		S2[index+1] = S2[index] + dt*(-S2[index]/taus + (1 - S2[index])*gamma*H2[index])
	return H1[1:], H2[1:],S1[1:], S2[1:]

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

def quasi_dist(times,T):
	s1 = []
	s2 = []
	r1 = []
	r2 = []
	for _ in range(times):
		result = experiment()
		s1.append(result[2][T-1:-1])
		s2.append(result[3][T-1:-1])
		r1.append(result[0][T-1:-1])
		r2.append(result[1][T-1:-1])
	s1=np.hstack(s1)
	s2=np.hstack(s2)
	r1=np.hstack(r1)
	r2=np.hstack(r2)

	return s1,s2,r1,r2


'''plt.figure
cprime = 20
starttime = -0.1
endtime = 0.1
time = np.linspace(starttime, endtime, steps)
result = experiment()
plt.plot(time*1000, smoothing(result[0]), color = 'red')
plt.plot(time*1000, smoothing(result[1]), '--', color = 'blue')
plt.plot(time*1000, 15*np.ones(steps))
plt.xlabel('Time(ms)')
plt.ylabel('Firing rate(Hz)')
plt.ylim(top = 20)
plt.text(-75, 16, 'Threshold')
plt.show()'''

plt.figure
for i in range(1):
	result = experiment()
	if cprime == 0: hue = ['red','blue']
	else: hue = ['blue','green']
	plt.plot(time*1000, smoothing(result[0]), color = hue[0])
	plt.plot(time*1000, smoothing(result[1]), color = hue[1])
plt.plot(time*1000, 15*np.ones(steps),color = 'black')
plt.xlabel('Time(ms)')
plt.ylabel('Firing rate(Hz)')
plt.ylim(top = 20)
plt.text(-125, 16, 'Threshold')
plt.show()

plt.figure
for i in range(10):
	result = experiment()
	hue = ['red','blue']
	plt.plot(time*1000, smoothing(result[2]), color = hue[0])
	plt.plot(time*1000, smoothing(result[3]), color = hue[1])
plt.xlabel('Time(ms)')
plt.ylabel('Gate variable/S(%)')
plt.show()


n_s=nullclines(0,30)

# plot S
plt.figure
plt.scatter(n_s[0][:,0],n_s[0][:,1],s=1)
plt.scatter(n_s[1][:,0],n_s[1][:,1],s=1)
plt.scatter(n_s[2][:,0],n_s[2][:,1])

for i in range(10):
	result = experiment()
	plt.scatter(result[2],result[3],s=1)

plt.xlabel('S2')
plt.ylabel('S1')
plt.legend(['dS2/dt=0','dS1/dt=0','Fixed Points'])
plt.title('Nullclines cprime: {}'.format(cprime))
plt.show()

s1,s2,r1,r2 = quasi_dist(500,10000)

plt.figure
h=plt.hist2d(x=s1, y=s2, bins=100,cmin=0)
plt.xlabel('S2')
plt.ylabel('S1')
plt.show()

plt.figure
h_r=plt.hist2d(x=r1, y=r2, bins=100,cmin=0)
plt.xlabel('r2')
plt.ylabel('r1')
plt.show()

#plt.contourf(np.array(h[0]))
#plt.pcolormesh(h[0])
#range=[[0,0.75],[0,0.75]]

sio.savemat('hist_30_ran.mat', mdict={'hist':h[0],'hist_r':h_r[0],'edges':[h[1],h[2]],'edges_r':[h_r[1],h_r[2]], 'S': [s1,s2],'r':[r1,r2]})

