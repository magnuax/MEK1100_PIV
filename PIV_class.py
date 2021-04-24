"""
Class for analyzing measurements from PIV (Particle Imaging Velocimetry),
Dependencies: scipy.io, numpy, matplotlib.pyplot
@author Magnus Axelsen
@version 1.0, 21/04/21
"""
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch

class PIV:
    def __init__(self, datafile):
        data = sio.loadmat(datafile)
        self.x = data.get("x")
        self.y = data.get("y")
        self.u = data.get("u")
        self.v = data.get("v")
        self.xit = data.get("xit")
        self.yit = data.get("yit")

        self.x_start = self.x[0,0]; self.x_stop = self.x[0,-1]
        self.y_start = self.y[0,0]; self.y_stop = self.y[-1,0]
        self.x_len = self.x.shape[1]; self.dx = self.x[0,1] - self.x[0,0]
        self.y_len = self.y.shape[0]; self.dy = self.y[1,0] - self.y[0,0]

        self.speed = np.sqrt(self.u**2 + self.v**2)/1000 #[m/s]


    def __call__(self):
        print(f"Matrice shapes, formatted (ny, nx):\n \
            x:{self.x.shape}    y:{self.y.shape}\n \
            u:{self.u.shape}    v:{self.v.shape}\n \
            xit:{self.xit.shape}    yit:{self.yit.shape}")
        print(f"\n x = \n{self.x}")
        print(f"\n y = \n{self.y}")


    def test_spacing(self, expected):
        msg = "Grid unevenly spaced! calculations may be incorrect!"
        for i in range(self.y_len-1):
            for j in range(self.x_len-1):
                assert (self.x[i,j+1]-self.x[i,j] == expected), msg
                assert (self.y[i+1,j]-self.y[i,j] == expected), msg


    def plot_wrap(self, splitcolor="crimson"):
        plt.plot(self.xit[0], self.yit[0], splitcolor, linewidth=3)
        plt.axis([self.x_start,self.x_stop,self.y_start,self.y_stop])
        plt.xlabel("x [mm]")
        plt.ylabel("y [mm]")
        plt.tight_layout()


    def split(self, vector):
        vec_gass = np.full_like(vector, np.nan)
        vec_fluid = np.full_like(vector, np.nan)
        for i in range(200):
            for j in range(193):
                if self.y[i,j]>self.yit[0,j]:
                    vec_gass[i,j] = vector[i,j]
                else:
                    vec_fluid[i,j] = vector[i,j]
        return vec_gass, vec_fluid


    def plot_speed(self, levels):
        speed_gass, speed_fluid = self.split(self.speed)
        plt.figure(figsize=(8,5))

        plt.contourf(self.x, self.y, speed_gass, levels,  cmap="viridis")
        plt.colorbar(label="\nGass speed, [m/s]\n")
        plt.contourf(self.x, self.y, speed_fluid, levels, cmap="magma")
        plt.colorbar(label="\nFluid speed, [m/s] ")
        self.plot_wrap()


    def plot_velocity(self, skipsize):
        n = skipsize
        L = self.speed[::n,::n]
        U, V  = self.u[::n,::n]/L, self.v[::n,::n]/L
        X, Y  = self.x[::n,::n], self.y[::n,::n]

        plt.figure(figsize=(8,5))
        plt.quiver(X, Y, U, V, L, cmap="plasma")
        plt.colorbar(label=rf"||$\vec v$||,  [m/s]")
        self.plot_wrap()


    def plot_div(self, levels):
        dudx = np.gradient(self.u, self.dx, axis=1)
        dvdy = np.gradient(self.v, self.dy, axis=0)
        div  = dudx + dvdy

        plt.figure(figsize=(8,5))
        plt.contourf(self.x,self.y, div, levels)
        plt.colorbar(label=r"$\nabla\cdot\vec v$,    $[s^{-1}]$ ")
        self.plot_wrap()


    def plot_curl(self, levels):
        dudy = np.gradient(self.u, self.dy, axis=0)
        dvdx = np.gradient(self.v, self.dx, axis=1)
        curl  = dvdx - dudy

        plt.figure(figsize=(8,5))
        plt.contourf(self.x, self.y, curl, levels)
        plt.colorbar(label=r"$(\nabla\times\vec v) \cdot \hat k$,  $[s^{-1}]$ ")
        self.plot_wrap()


    def streamplot(self, density):
        plt.figure(figsize=(8,5))
        plt.streamplot(self.x, self.y, self.u, self.v, density=density, color=self.speed)
        plt.colorbar(label="\nSpeed, [m/s] ")
        self.plot_wrap()
