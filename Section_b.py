import numpy as np
import matplotlib.pyplot as plt




#SECTION b)


#Initial set up

V = 1e13




#We define a function that calculate the Helmholtz free energy divided by T

def F(n, fill):

	#First we calculate the total number of rods

	N = n * V


	#After that, we create a vector of values for both Nx and Ny. They both can take values between zero and the total number of rods N

	Nx = np.linspace(0, N, 1000)

	Ny = np.linspace(0, N, 1000)


	#Now we transform those vectors into matrices in order to be able to represent the contourmap

	NX, NY = np.meshgrid(Nx, Ny)


	#We also create a matrix to store the values of F/T. Instead of filling it with zeros, we fill it with a value specified at the beggining of the program in order to adjust the scale of the elements of the matrix that will not be filled (those that correspond to impossible combinations of Nx and Ny)

	f = np.ones([1000,1000]) * fill


	#We calculate F/T for all the possible combinations of Nx and Ny

	for i in range (0, 1000):

		for j in range (0, int(1000-i)):

			#The several if statements are used for not calculate values of log that are not valid

			if NX[i,j] == 0 and NY[i,j] != 0 and (N - NY[i,j]) != 0:

				f[i,j] = NY[i,j]*np.log(NY[i,j]) + (N - NY[i,j])*np.log(N - NY[i,j]) - N*np.log(V) + 10*n*NY[i,j] - 10*(1/V)*(NY[i,j]**2)


			elif NX[i,j] != 0 and NY[i,j] == 0 and (N - NX[i,j]) != 0:

				f[i,j] = NX[i,j]*np.log(NX[i,j]) + (N - NX[i,j])*np.log(N - NX[i,j]) - N*np.log(V) + 10*n*(NX[i,j]) - 10*(1/V)*(NX[i,j]**2)


			elif NX[i,j] != 0 and NY[i,j] != 0 and (N - NX[i,j] - NY[i,j]) == 0:

				f[i,j] = NX[i,j]*np.log(NX[i,j]) + NY[i,j]*np.log(NY[i,j]) - N*np.log(V) + 10*n*(NX[i,j]+NY[i,j]) - 10*(1/V)*( (NX[i,j]**2) + (NY[i,j]**2) + (NX[i,j]*NY[i,j]) )


			elif NX[i,j] == 0 and NY[i,j] != 0 and (N - NX[i,j] - NY[i,j]) == 0:

				f[i,j] = NY[i,j]*np.log(NY[i,j]) - N*np.log(V) + 10*n*(NY[i,j]) - 10*(1/V)*(NY[i,j]**2)


			elif NX[i,j] != 0 and NY[i,j] == 0 and (N - NX[i,j] - NY[i,j]) == 0:

				f[i,j] = NX[i,j]*np.log(NX[i,j]) - N*np.log(V) + 10*n*(NX[i,j]) - 10*(1/V)*(NX[i,j]**2)


			elif NX[i,j] == 0 and NY[i,j] == 0:

				f[i,j] = N*np.log(N) - N*np.log(V)

			else:

				f[i,j] = NX[i,j]*np.log(NX[i,j]) + NY[i,j]*np.log(NY[i,j]) + (N - NX[i,j] - NY[i,j])*np.log(N - NX[i,j] - NY[i,j]) - N*np.log(V) + 10*n*(NX[i,j]+NY[i,j]) - 10*(1/V)*( (NX[i,j]**2) + (NY[i,j]**2) + (NX[i,j]*NY[i,j]) )


	return(f, NX, NY)






#Let's plot F/T for several values of n

plt.contourf(F(1e-3, -7.46e10)[1], F(1e-3, -7.46e10)[2], F(1e-3, -7.46e10)[0], 50 , cmap = 'RdGy')
plt.colorbar()
plt.xlabel(r'$N_x$')
plt.ylabel(r'$N_y$')
plt.title(r'F/T vs $N_x$ and $N_y$ $(n = 1 \cdot 10^{-3})$')
plt.savefig("F_n=1e-3.pdf")


plt.show()
plt.contourf(F(0.27, -3.8e12)[1], F(0.27, -3.8e12)[2], F(0.27, -3.8e12)[0], 50 , cmap = 'RdGy')
plt.colorbar()
plt.xlabel(r'$N_x$')
plt.ylabel(r'$N_y$')
plt.title(r'F/T vs $N_x$ and $N_y$ (n = 0.27)')
plt.savefig("F_n=0.27.pdf")


plt.show()
plt.contourf(F(0.28, -3.785e12)[1], F(0.28, -3.785e12)[2], F(0.28, -3.785e12)[0], 50 , cmap = 'RdGy')
plt.colorbar()
plt.xlabel(r'$N_x$')
plt.ylabel(r'$N_y$')
plt.title(r'F/T vs $N_x$ and $N_y$ (n = 0.28)')
plt.savefig("F_n=0.28.pdf")


plt.show()
plt.contourf(F(0.4, -3.14e12)[1], F(0.4, -3.14e12)[2], F(0.4, -3.14e12)[0], 50 , cmap = 'RdGy')
plt.colorbar()
plt.xlabel(r'$N_x$')
plt.ylabel(r'$N_y$')
plt.title(r'F/T vs $N_x$ and $N_y$ (n = 0.4)')
plt.savefig("F_n=0.4.pdf")


plt.show()
plt.contourf(F(1, 1.15e13)[1], F(1, 1.15e13)[2], F(1, 1.15e13)[0], 50 , cmap = 'RdGy')
plt.colorbar()
plt.xlabel(r'$N_x$')
plt.ylabel(r'$N_y$')
plt.title(r'F/T vs $N_x$ and $N_y$ (n = 1)')
plt.savefig("F_n=1.pdf")






#We can see that there, from one n on, we observe three minima instead of one that are clearly on one of the Ni being equal to N and the other ones zero

#But, for the cases where there is just one minima, that will appears to be in the point where Nx = Ny = Nz = N/3. Let's check that drawing lines at those quantities 

N = (1e-3) * V

Nx = Ny = N/3


plt.show()
plt.contourf( F(1e-3, -7.46e10)[1], F(1e-3, -7.46e10)[2], F(1e-3, -7.46e10)[0], 50 , cmap = 'RdGy' )
plt.vlines( Nx, 0, Ny, linewidth = 2, color = 'black' )
plt.hlines( Ny, 0, Nx, linewidth=2, color ='black' )
plt.colorbar()
plt.xlabel(r'$N_x$')
plt.ylabel(r'$N_y$')
plt.title(r'F/T vs $N_x$ and $N_y$ $(n = 1 \cdot 10^{-3})$')
plt.savefig("F_n=1e-3_lines.pdf")



#Now, we are going to plot F/T with respect to n

n = 0.28

n1 = np.linspace(0, n, 1000)

FT1 = V*(n1 * np.log(n1/3) + (10/3) * n1**2)


n2 = np.linspace(n, 1, 1000)

FT2 = V*n2*np.log(n2)


ymin = V*(n * np.log(n/3) + (10/3) * n**2)

ymax = V*n*np.log(n)


plt.show()
plt.plot(n1, FT1, color = 'hotpink')
plt.plot(n2, FT2, color = 'hotpink')
plt.vlines( n, ymin, ymax, ls = 'dashed', color = 'hotpink')
plt.xlabel(r'n/(rods/m^3)')
plt.ylabel(r'(F/T)/rods')
plt.title(r'F/T  vs n')
plt.grid(True)
plt.savefig("F_n.pdf")




#Finally, we are going to plot P/T vs n

PT1 = n1 + (10/3) * n1**2

PT2 = n2

ymin = n + (10/3) * n**2
ymax = n


plt.show()
plt.plot(n1, PT1, color = 'hotpink')
plt.plot(n2, PT2, color = 'hotpink')
plt.vlines( n, ymin, ymax, ls = 'dashed', color = 'hotpink')
plt.xlabel(r'n/(rods/m^3)')
plt.ylabel(r'(P/T)/(rods / JPa)')
plt.title(r'P/T  vs n')
plt.grid(True)
plt.savefig("P_n.pdf")
