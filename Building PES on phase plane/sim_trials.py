import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

## parameters set
a = 270
b = 108
d = 0.154
gamma = 0.641
taus = 100*10**-3
tauampa = 2*10**-3
Jn11, Jn12 = 0.2609, 0.0497
Jext = 5.2*(10**-4)
I0 = 0.3255
stdnoise = 0.5
mu0 = 30
dt = 10**-4
cprime = 0
trial = 1

starttime = -0.5
endtime = 2
steps = int(abs(starttime - endtime)/dt)
time = np.linspace(starttime, endtime, steps)

## Functions 
def experiment():
	H1, H2, S1, S2 = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
	#H1, H2 = np.zeros(steps+1), np.zeros(steps+1)
	#S1, S2 = np.random.rand(steps+1), np.random.rand(steps+1)
	Inoise1, Inoise2 = np.zeros(steps+1), np.zeros(steps+1)
	for index, t in enumerate(time):
		Inoise1[index+1] = Inoise1[index] + dt*(-Inoise1[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		Inoise2[index+1] = Inoise2[index] + dt*(-Inoise2[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		if t > 0:
			x1 = Jn11*S1[index] - Jn12*S2[index] + I0 + Jext*mu0*(1-cprime/100) + Inoise1[index]
			x2 = Jn11*S2[index] - Jn12*S1[index] + I0 + Jext*mu0*(1+cprime/100) + Inoise2[index]
		else:
			x1 = Jn11*S1[index] - Jn12*S2[index] + I0 + Inoise1[index]
			x2 = Jn11*S2[index] - Jn12*S1[index] + I0 + Inoise2[index]
		H1[index+1] = (a*x1-b)/(1-np.exp(-d*(a*x1-b)))
		H2[index+1] = (a*x2-b)/(1-np.exp(-d*(a*x2-b)))
		S1[index+1] = S1[index] + dt*(-S1[index]/taus + (1 - S1[index])*gamma*H1[index])
		S2[index+1] = S2[index] + dt*(-S2[index]/taus + (1 - S2[index])*gamma*H2[index])
	return H1[1:], H2[1:], S1[1:], S2[1:]

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


# Not an efficient method will try solving simultaneous equations
# and using sympy to solve
def nullclines(cprime,mu0):

	def gate(S):
		# Input Currents 
		# Zero on the end is for easy moving to zero input
		# I[2]
		I = [Jext*mu0*(1+(cprime/100)), 
			 Jext*mu0*(1-(cprime/100)),0]
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
			I = [Jext*mu0*(1+(cprime/100)), 
			 	Jext*mu0*(1-(cprime/100)),0]
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

def devi(times):
	s1=0
	for _ in range(times):
		result = experiment()
		if result[2][-1]>result[3][-1]: s1+=1
	return s1/times

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

## Firing rate r
plt.figure
for i in range(trial):
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

## Gate variable S
plt.figure
for i in range(trial):
	result = experiment()
	hue = ['blue','red']
	plt.plot(time*1000, smoothing(result[2]), color = hue[0])
	plt.plot(time*1000, smoothing(result[3]), color = hue[1])
plt.xlabel('Time(ms)')
plt.ylabel('Gate variable/S(%)')
plt.show()

## Nullclines	
n_s=nullclines(cprime)

# plot S
plt.figure
plt.scatter(n_s[0][:,0],n_s[0][:,1],s=1)
plt.scatter(n_s[1][:,0],n_s[1][:,1],s=1)
plt.scatter(n_s[2][:,0],n_s[2][:,1])

for i in range(trial):
	result = experiment()
	plt.scatter(result[2],result[3],s=1)

plt.xlabel('S2')
plt.ylabel('S1')
plt.legend(['dS2/dt=0','dS1/dt=0','Fixed Points'])
plt.title('Nullclines cprime: {}'.format(cprime))
plt.show()

# plot r
plt.figure
plt.scatter(n_s[3][:,0],n_s[3][:,1],s=1)
plt.scatter(n_s[4][:,0],n_s[4][:,1],s=1)
plt.scatter(n_s[5][:,0],n_s[5][:,1])

for i in range(trial):
	result = experiment()		
	plt.scatter(smoothing(result[0]),smoothing(result[1]),s=1)

plt.xlabel('r2')
plt.ylabel('r1')
plt.legend(['dS2/dt=0','dS1/dt=0','Fixed Points'])
plt.title('Nullclines cprime: {}'.format(cprime))
plt.show()

## Quasi-distribution

s1,s2,r1,r2 = quasi_dist(5000,2000)

plt.figure
h=plt.hist2d(x=s1, y=s2, bins=100,cmin=100)
plt.xlabel('S2')
plt.ylabel('S1')
plt.show()

plt.figure
h_r=plt.hist2d(x=r1, y=r2, bins=100,cmin=100)
plt.xlabel('r2')
plt.ylabel('r1')
plt.show()

#plt.contourf(np.array(h[0]))
#plt.pcolormesh(h[0])
#range=[[0,0.75],[0,0.75]]

sio.savemat('hist_s30_0.mat', mdict={'hist':h[0],'hist_r':h_r[0],'edges':[h[1],h[2]],'edges_r':[h_r[1],h_r[2]], 'S': [s1,s2],'r':[r1,r2]})

'''plt.figure
plt.hist2d(x=s2, y=s1, bins=100,range=[[0,0.2],[0.55,0.75]])
plt.xlabel('S2')
plt.ylabel('S1')
plt.show()

plt.figure
plt.hist2d(x=s2, y=s1, bins=100,range=[[0.55,0.75],[0,0.2]])
plt.xlabel('S2')
plt.ylabel('S1')
plt.show()'''

cprime=70
p=devi(10000)
p