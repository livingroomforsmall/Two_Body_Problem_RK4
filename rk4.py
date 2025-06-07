import numpy as np
import matplotlib.pyplot as plt

G = 100
M = 20
u0 = np.array([10.0, 5.00, 10.0, 10.0])
E0 = 1/2*(u0[2]**2 + u0[3]**2) - G*M/(np.sqrt(u0[0]**2+u0[1]**2))

def du(t, u):
    r = np.sqrt(u[0]**2 + u[1]**2)
    return np.array([u[2], u[3], -G*M*u[0]/r**3, -G*M*u[1]/r**3])

def rk4singlestep(fun, dt, t0, y0):
    f1 = fun(t0, y0)
    f2 = fun(t0 + dt/2, y0 + (dt/2)*f1)
    f3 = fun(t0 + dt/2, y0 + (dt/2)*f2)
    f4 = fun(t0 + dt/2, y0 + dt*f3)
    yout = y0 + (dt/6) * (f1 + 2*f2 + 2*f3 + f4)
    return yout

dt = 0.01 #LRLRLRLRLRLRLRLRLLRLRLRR
T = 10
num_time_pts = int(T/dt)
t = np.linspace(0, T, num_time_pts) #von 0 bis T-0.01
E = np.zeros(num_time_pts)
E[0] = E0
U = np.zeros((4, num_time_pts))
U[:,0] = u0
uin = u0
for _ in range(num_time_pts - 1):
    uout = rk4singlestep(du, dt, t[_], uin)
    U[:, _+1] = uout
    E[_+1] = 1/2*(uout[2]**2 + uout[3]**2) - G*M/(np.sqrt(uout[0]**2+uout[1]**2))
    uin = uout
#Matrix U hat rows x,y,vx,vy, now plot 1st and 2nd row to get trajectory
plt.plot(U[0,:], U[1,:], 'b')
plt.title('Trajectory 2D orbit system')
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.show()

plt.plot(t, E, 'r')
plt.title('Total energy')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.show()

"""
Erwartung ist ne geschlossene Kurve und const E, allerdings schon gut
im Vergleich zu Euler weil Euler Restglied O(delta2) RK4 O(delta5)
also für kleine Zeitverrückung ist Fehler deutlich smaller
Erst nach 1500 steps divergiert die Kurve... Man erkennt schon
wo die Zentralmasse ist.
Wahrscheinlich ist 0.01 nicht klein genug habe ich noch kleiner angepasst
und geschlossene Kurven bekommen
"""