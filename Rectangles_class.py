import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt

class Rectangles:
    def __init__(self, PIV, rectangles):
        # PIV is instance of PIV class
        # rectangles a thrice nested list, rectangles=[[corners0], [corners1]]
        # for rectangle ABCD: corners=[[Ax,Ay], [Cx,Cy]], contains indices
        self.PIV = PIV
        self.rectangles=rectangles


    def plot(self):
        for corners in self.rectangles:
            Ax, Ay = corners[0]
            Cx, Cy = corners[1]
            x = self.PIV.x
            y = self.PIV.y
            AB = [ [x[Ay,Ax], x[Cy,Cx]], [y[Ay,Ax], y[Ay,Ax]] ]
            BC = [ [x[Cy,Cx], x[Cy,Cx]], [y[Ay,Ax], y[Cy,Cx]] ]
            CD = [ [x[Cy,Cx], x[Cy,Ax]], [y[Cy,Cx], y[Cy,Cx]] ]
            DA = [ [x[Cy,Ax], x[Cy,Ax]], [y[Cy,Cx], y[Ay,Ax]] ]

            plt.plot(AB[0], AB[1], "tab:red", linewidth=2)
            plt.plot(BC[0], BC[1], "tab:green", linewidth=2)
            plt.plot(CD[0], CD[1], "tab:blue", linewidth=2)
            plt.plot(DA[0], DA[1], "k", linewidth=2)


    def circulation(self):
        u = self.PIV.u
        v = self.PIV.v
        dx = self.PIV.dx
        dy = self.PIV.dy
        circ_total = []
        circ_sides = []

        for corners in self.rectangles:
            Ax, Ay = corners[0]
            Cx, Cy = corners[1]

            intAB = np.sum( u[Ay, Ax:Cx+1] ) * dx
            intBC = np.sum( v[Ay:Cy+1, Cx] ) * dy
            intCD = np.sum( u[Cy, Ax:Cx+1] ) * (-dx)
            intDA = np.sum( v[Ay:Cy+1, Ax] ) * (-dy)

            circ_total.append(intAB + intBC + intCD + intDA)
            circ_sides.append([intAB, intBC, intCD, intDA])

        return np.array(circ_total), np.array(circ_sides)


    def circulation_stokes(self):
        u = self.PIV.u
        v = self.PIV.v
        dx = self.PIV.dx
        dy = self.PIV.dy
        dudy = np.gradient(u, dx, axis=0)
        dvdx = np.gradient(v, dy, axis=1)
        curl  = dvdx - dudy
        circulation = []

        for corners in self.rectangles:
            Ax, Ay = corners[0]
            Cx, Cy = corners[1]
            int = 0 # [m^2 /s ]
            for i in range(Ay,Cy+1):
                for j in range(Ax,Cx+1):
                    int += curl[i,j]*dx*dy
            circulation.append(int)

        return np.array(circulation)


    def flux(self):
        u = self.PIV.u
        v = self.PIV.v
        dx = self.PIV.dx
        dy = self.PIV.dy
        flux_total = []
        flux_sides = []

        for corners in self.rectangles:
            Ax, Ay = corners[0]
            Cx, Cy = corners[1]

            intAB = np.sum(-v[Ay, Ax:Cx+1] )*dx
            intBC = np.sum( u[Ay:Cy+1, Cx] )*dy
            intCD = np.sum( v[Cy, Ax:Cx+1] )*(-dx)
            intDA = np.sum(-u[Ay:Cy+1, Ax] )*(-dy)

            flux_total.append(intAB + intBC + intCD + intDA)
            flux_sides.append([intAB, intBC, intCD, intDA])

        return np.array(flux_total), np.array(flux_sides)
