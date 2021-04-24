from PIV_class import PIV
from Rectangles_class import Rectangles
import matplotlib.pyplot as plt
plt.style.use("ggplot")
import matplotlib.patches as mpatch

# Setup
indices = [[[35,160],[70,170]], [[35,85],[70,100]], [[35,50],[70,60]]]
levels=200
tube = PIV("data.mat")
rectangles = Rectangles(tube, indices)

# Oppgave a)
tube()
tube.test_spacing(0.5)

# Oppgave b)
tube.plot_speed(levels)

# Oppgave c)
tube.plot_velocity(skipsize=9)
rectangles.plot()

# Oppgave d)
tube.plot_div(levels)
rectangles.plot()

# Oppgave e)
tube.plot_curl(levels)
rectangles.plot()
tube.streamplot(density=0.8)

circ_direct, circ_sides = rectangles.circulation()
circ_stokes = rectangles.circulation_stokes()

for i in range(0,3):
    print(f"\nRectangle {i+1} circulation, [mm^2/s] :")
    print(f"   Direct: {circ_direct[i]:.3f},\t (AB:{circ_sides[i,0]:.3f}, BC:{circ_sides[i,1]:.3f}, CD:{circ_sides[i,2]:.3f}, DA:{circ_sides[i,3]:.3f} ")
    print(f"   Stokes: {circ_stokes[i]:.3f}\n")

# Oppgave g)
flux_total, flux_sides = rectangles.flux()
for i in range(0,3):
    print(f"\nRectangle {i+1} flux, [mm^3/s]:")
    print(f"   Direct: {flux_total[i]:.3f}\t (AB={flux_sides[i,0]:.3f}, BC={flux_sides[i,1]:.3f}, CD={flux_sides[i,2]:.3f}, DA={flux_sides[i,3]:.3f})")
plt.show()
