'''
	Is it necessary that neurons subserving integration during simulation also
	show persistent activity during working memory, Article asks this question explicitly 
	on page 2

	It is necessary that the neurons continue to fire in the attractor sites
	during working memory

	This version of phase plane is without AMPA.
'''
import matplotlib.pyplot as plt
import numpy as np
from sympy.solvers import solve
from scipy.optimize import fsolve 

#=================================
#			Constants
#=================================
a = 270
b = 108
d = 0.154
gamma = 0.641
taus = 100*10**-3
Jn11, Jn12 = 0.2609, 0.0497
Jext = 5.2*(10**-4)
I0 = 0.3255


#=================================
# 			Phase Plane
#=================================
'''
H1,H2,S1, S2,sig = experiment()
plt.plot(S2,S1)
plt.xlabel('S2')
plt.ylabel('S1')
plt.title('Phase Plane')
plt.show()

'''
# Plotting Nuclines for phase plane analysis
# dS1/dt == 0
# dS2/dt == 0
# Since S1 and S2 cant be solved directly need numerical solutions

# Not an efficient method will try solving simultaneous equations
# and using sympy to solve
def nullclines(cprime,mu0):

	def gate(S):
		# Input Currents 
		# Zero on the end is for easy moving to zero input
		# I[2]
		I = [Jext*mu0*(1-(cprime/100)), 
			 Jext*mu0*(1+(cprime/100)),0]
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
			I = [Jext*mu0*(1-(cprime/100)), 
			 	Jext*mu0*(1+(cprime/100)),0]
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
	
	# Stack the results for easy plotting	
		
	ds1=np.vstack(ds1)
	ds2=np.vstack(ds2)
	fix_points=np.vstack(fix_points)

	dr1=np.vstack(dr1)
	dr2=np.vstack(dr2)
	fix_points_r=np.vstack(fix_points_r)

	plt.scatter(ds1[:,0],ds1[:,1],s=1)
	plt.scatter(ds2[:,0],ds2[:,1],s=1)
	plt.scatter(fix_points[:,0],fix_points[:,1])
	plt.xlabel('S2')
	plt.ylabel('S1')
	plt.legend(['dS2/dt=0','dS1/dt=0','Fixed Points'])
	plt.title('Nullclines cprime: {}'.format(cprime))
	plt.show()

	plt.scatter(dr1[:,0],dr1[:,1],s=1)
	plt.scatter(dr2[:,0],dr2[:,1],s=1)
	plt.scatter(fix_points_r[:,0],fix_points_r[:,1])
	plt.xlabel('r2')
	plt.ylabel('r1')
	plt.legend(['dS2/dt=0','dS1/dt=0','Fixed Points'])
	plt.title('Nullclines cprime: {}'.format(cprime))
	plt.show()
	return ds1,ds2,fix_points

cprime=6.4
mu=30
n_s=nullclines(cprime,mu)


plt.figure('Figure Object 1',       # 图形对象名称  窗口左上角显示
           figsize = (6, 6),        # 窗口大小
           dpi = 120,               # 分辨率
           facecolor = 'white', # 背景色
           )
plt.scatter(n_s[0][:,0],n_s[0][:,1],s=1)
plt.scatter(n_s[1][:,0],n_s[1][:,1],s=1)
plt.scatter(n_s[2][:,0],n_s[2][:,1])

for i in range(5):
	result = experiment()
	plt.scatter(result[2],result[3],s=1)

plt.xlabel('S2')
plt.ylabel('S1')
plt.legend(['dS2/dt=0','dS1/dt=0','Fixed Points'])
plt.title('Nullclines cprime: {}'.format(cprime))
plt.show()
