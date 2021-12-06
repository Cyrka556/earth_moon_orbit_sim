# Imports
import numpy as np
import matplotlib.pyplot as plt


# Constants
pi = np.pi
G = 6.67408 * (10 ** (-11))  # m^3/kg s
Me = 5.972 * (10 ** 24)  # kg
Ms = 2755 # kg
Mm = 7.34767309 * (10 ** 22)
Re = 6371000 # m
Rm = 1737100
xm = -384400000


# Initial Conditions x0 and v0 defined in the function EarthMoon
y0 = 0  # m
vx0 = 0  # m/s
t0 = 0  # s




def acc(x, y):
    return [-G * Me * x / ((x ** 2 + y ** 2) ** (3 / 2)), -G * Me * y / ((x ** 2 + y ** 2) ** (3 / 2))]



def acc13(x, y):
    return [-((G * Me * x / ((x**2 + y**2) ** (3 / 2))) + (G * Mm * (x - xm) / ((x - xm)**2 + y**2)**(3/2))), -((G * Me * y / ((x**2 + y**2) ** (3/2))) + (G * Mm * y / ((x - xm)**2 + y**2)**(3/2)))]



# Defines r and theta from Earth
def polarE(x, y):
    rvector = np.sqrt(x ** 2 + y ** 2)


    if x == 0:
        if y > 0:
            angle = pi/2

        else:
            angle = 3*pi/2

    else:
        if x > 0:
            if y > 0:
                angle = np.arctan(np.abs(y/x))

            else:
                angle = 2*pi - np.arctan(np.abs(y/x))


        else:
            if y > 0:
                angle = pi - np.arctan(np.abs(y/x))

            else:
                angle = pi + np.arctan(np.abs(y/x))


    return [rvector, angle]



# Defines r and theta from the Moon
def polarM(x, y):
    x -= xm
    rvector = np.sqrt(x ** 2 + y ** 2)


    if x == 0:
        if y > 0:
            angle = pi/2

        else:
            angle = 3*pi/2

    else:
        if x > 0:
            if y > 0:
                angle = np.arctan(np.abs(y/x))

            else:
                angle = 2*pi - np.arctan(np.abs(y/x))


        else:
            if y > 0:
                angle = pi - np.arctan(np.abs(y/x))

            else:
                angle = pi + np.arctan(np.abs(y/x))


    return [rvector, angle]




def potential(x, y):
    return -G * Me * Ms / np.sqrt(x**2 + y**2)



def Mpotential(x, y):
    return -G * Ms * ((Me / np.sqrt(x**2 + y**2)) + (Mm / np.sqrt((x - xm)**2 + y**2)))




def Earth(x0, vy0):
    """4th order Runge-Kutta method for solving orbits, accelaration calculated is for just the earth, x0 and vy0 (initial
    conditions) defined on input so this can be generalised for circular and elliptical orbits"""

    a = 0.03
    x = [x0]
    y = [y0]
    vx = [vx0]
    vy = [vy0]
    t = [t0]



    Etheta = 0
    while Etheta < 2*pi - a:
        # Adjusts h based on distance, as calculation was taking too long when far away
        h = 8 * 10 ** (-7) * polarE(x[-1], y[-1])[0] + 0.1



        # k arrays
        kcoef = np.array([1, 2, 2, 1])
        kx = np.zeros(4)
        ky = np.zeros(4)
        kvx = np.zeros(4)
        kvy = np.zeros(4)



        # k1
        kx[0] = vx[-1]
        ky[0] = vy[-1]
        kvx[0] = acc(x[-1], y[-1])[0]
        kvy[0] = acc(x[-1], y[-1])[1]



        # k2
        kx[1] = vx[-1] + h * kvx[0] / 2
        ky[1] = vy[-1] + h * kvy[0] / 2
        kvx[1] = acc(x[-1] + h * kx[0] / 2, y[-1] + h * ky[0] / 2)[0]
        kvy[1] = acc(x[-1] + h * kx[0] / 2, y[-1] + h * ky[0] / 2)[1]

        # print(kx, ky, kvx, kvy, '\n')

        # k3
        kx[2] = vx[-1] + h * kvx[1] / 2
        ky[2] = vy[-1] + h * kvy[1] / 2
        kvx[2] = acc(x[-1] + h * kx[1] / 2, y[-1] + h * ky[1] / 2)[0]
        kvy[2] = acc(x[-1] + h * kx[1] / 2, y[-1] + h * ky[1] / 2)[1]



        # k4
        kx[3] = vx[-1] + h * kvx[2]
        ky[3] = vy[-1] + h * kvy[2]
        kvx[3] = acc(x[-1] + h * kx[2], y[-1] + h * ky[2])[0]
        kvy[3] = acc(x[-1] + h * kx[2], y[-1] + h * ky[2])[1]



        # Appending New x and y
        xi1 = x[-1] + h * (np.sum(kcoef * kx)) / 6
        x.append(xi1)
        yi1 = y[-1] + h * (np.sum(kcoef * ky)) / 6
        y.append(yi1)

        # Appending New Velocity Values
        vxi1 = vx[-1] + h * (np.sum(kcoef * kvx)) / 6
        vx.append(vxi1)
        vyi1 = vy[-1] + h * (np.sum(kcoef * kvy)) / 6
        vy.append(vyi1)

        # Find Polar Coordinate of Rocket Relative to Earth
        Edist = polarE(xi1, yi1)[0]
        Etheta = polarE(xi1, yi1)[1]


        # Appending Time
        t.append(t[-1] + h)



        if Edist < Re:
            print('The Rocket Crashed in to the Earth\n')
            break




    return [x, y, vx, vy, t]






def EarthMoon(x0, vy0):
    """4th order Runge-Kutta method for Earth-Moon orbit, accelaration calculated is for the earth and moon, x0 and vy0 (initial
        conditions) defined on input to established appropriate conditions for the orbit"""

    a = 0.03
    x = [x0]
    y = [y0]
    vx = [vx0]
    vy = [vy0]
    t = [t0]



    Etheta = 0
    while Etheta < 2*pi - a:
        # Adjusts h based on distance, as calculation was taking too long when far away
        if x[-1] > xm/2:
            h = 2 * 10 ** (-6) * polarE(x[-1], y[-1])[0] + 0.1

        else:
            h = 2 * 10 ** (-6) * polarM(x[-1], y[-1])[0] + 0.1


        # k arrays
        kcoef = np.array([1, 2, 2, 1])
        kx = np.zeros(4)
        ky = np.zeros(4)
        kvx = np.zeros(4)
        kvy = np.zeros(4)



        # k1
        kx[0] = vx[-1]
        ky[0] = vy[-1]
        kvx[0] = acc13(x[-1], y[-1])[0]
        kvy[0] = acc13(x[-1], y[-1])[1]



        # k2
        kx[1] = vx[-1] + h * kvx[0] / 2
        ky[1] = vy[-1] + h * kvy[0] / 2
        kvx[1] = acc13(x[-1] + h * kx[0] / 2, y[-1] + h * ky[0] / 2)[0]
        kvy[1] = acc13(x[-1] + h * kx[0] / 2, y[-1] + h * ky[0] / 2)[1]

        # print(kx, ky, kvx, kvy, '\n')

        # k3
        kx[2] = vx[-1] + h * kvx[1] / 2
        ky[2] = vy[-1] + h * kvy[1] / 2
        kvx[2] = acc13(x[-1] + h * kx[1] / 2, y[-1] + h * ky[1] / 2)[0]
        kvy[2] = acc13(x[-1] + h * kx[1] / 2, y[-1] + h * ky[1] / 2)[1]



        # k4
        kx[3] = vx[-1] + h * kvx[2]
        ky[3] = vy[-1] + h * kvy[2]
        kvx[3] = acc13(x[-1] + h * kx[2], y[-1] + h * ky[2])[0]
        kvy[3] = acc13(x[-1] + h * kx[2], y[-1] + h * ky[2])[1]



        # Appending New x and y
        xi1 = x[-1] + h * (np.sum(kcoef * kx)) / 6
        x.append(xi1)
        yi1 = y[-1] + h * (np.sum(kcoef * ky)) / 6
        y.append(yi1)

        # Appending New Velocity Values
        vxi1 = vx[-1] + h * (np.sum(kcoef * kvx)) / 6
        vx.append(vxi1)
        vyi1 = vy[-1] + h * (np.sum(kcoef * kvy)) / 6
        vy.append(vyi1)

        # Find Polar Coordinate of Rocket Relative to Earth
        Edist = polarE(xi1, yi1)[0]
        Etheta = polarE(xi1, yi1)[1]

        # Find Polar Coordinate of Rocket Relative to Moon
        Mdist = polarM(xi1, yi1)[0]




        # Appending Time
        t.append(t[-1] + h)



        if Edist < Re:
            print('The Rocket Crashed in to the Earth\n')
            break

        if Mdist < Rm:
            print('The Rocket Crashed in to the Moon\n')
            break


    return [x, y, vx, vy, t]





print('Welcome to my orbit simulation program\n')


ProgramLoop = 5
while ProgramLoop < 6:
    """Menu System"""
    print('Main Menu: Please select which model you would like\nRocket initial conditions set for each option to provide an example')
    print('a - Circular Earth Orbit\nb - Eccentric Earth Orbit\nc - Apollo 13 style Earth-Moon Orbit\nd - Earth Crash\ne - Moon Crash\nq - Quit')
    OrbitOption = input()
    print()

    if OrbitOption == 'a':
        """Circular Earth Orbit"""

        [x, y, vx, vy, t] = Earth(8000000, 7060)

        # Earth Figure
        earth = plt.Circle((0, 0), Re, fill=True, color='darkslateblue')


        KE = []
        PE = []
        Etot = []
        for i in range(len(t)):
            Kinetici = (1 / 2) * Ms * (vx[i] ** 2 + vy[i] ** 2)
            KE.append(Kinetici)
            Potentiali = potential(x[i], y[i])
            PE.append(Potentiali)
            Etot.append(Kinetici + Potentiali)

        print('Maximum Kinetic Energy is', max(KE))
        print('Minimum Kinetic Energy is', min(KE))
        print('The difference is', max(KE)-min(KE))
        print()
        print('Maximum Kinetic Energy is', max(Etot))
        print('Minimum Kinetic Energy is', min(Etot))
        print('The difference is', max(Etot) - min(Etot))
        print()
        print('Maximum time is', t[-1])
        print()


        fig1, orbit = plt.subplots(1)
        plt.plot(x, y, linestyle='--', linewidth='2', color='red')
        plt.axis('equal')
        plt.title('Circular orbit for body around the Earth')
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
        fig1 = plt.gcf()
        orbit = fig1.gca()
        orbit.add_patch(earth)


        fig2, earthenergy = plt.subplots(1)
        earthenergy.plot(t, KE, color='darkred')
        earthenergy.plot(t, PE, color='darkblue')
        earthenergy.plot(t, Etot, color='darkgreen')
        plt.ylabel('Energy (J)')
        plt.xlabel('Time (s)')
        earthenergy.legend(['Kinetic Energy', 'Potential Energy', 'Total Energy'], loc='best')

        fig3, earthkinetic = plt.subplots(1)
        earthkinetic.plot(t, KE)
        plt.ylabel('Kinetic Energy (J)')
        plt.xlabel('Time (s)')

        plt.show()




    elif OrbitOption == 'b':
        """Eccentric Earth Orbit"""

        [x, y, vx, vy, t] = Earth(8000000, 9600)

        # Earth Figure
        earth = plt.Circle((0, 0), Re, fill=True, color='darkslateblue')

        radius = []
        KE = []
        PE = []
        Etot = []
        for i in range(len(t)):
            Kinetici = (1 / 2) * Ms * (vx[i] ** 2 + vy[i] ** 2)
            KE.append(Kinetici)
            Potentiali = potential(x[i], y[i])
            PE.append(Potentiali)
            Etot.append(Kinetici + Potentiali)
            radius.append(np.sqrt((x[i])**2+(y[i])**2))

        print('Maximum radius is', max(radius))
        print('Minimum radius is', min(radius))
        print()
        print('The Maximum Energy is', max(Etot))
        print('The Minimum Energy is', min(Etot))
        print('The energy difference is', (max(Etot) - min(Etot)))
        print()

        fig1, orbit = plt.subplots(1)
        plt.plot(x, y, linestyle='--', linewidth='2', color='red')
        plt.axis('equal')
        plt.title('Eccentric Earth orbit')
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
        fig1 = plt.gcf()
        orbit = fig1.gca()
        orbit.add_patch(earth)

        fig2, earthenergy = plt.subplots(1)
        earthenergy.plot(t, KE, color='darkred')
        earthenergy.plot(t, PE, color='darkblue')
        earthenergy.plot(t, Etot, color='darkgreen')
        plt.ylabel('Energy (J)')
        plt.xlabel('Time (s)')
        earthenergy.legend(['Kinetic Energy', 'Potential Energy', 'Total Energy'], loc='best')

        plt.show()


    elif OrbitOption == 'c':
        """Apollo 13 style Earth-Moon Orbit"""

        [x, y, vx, vy, t] = EarthMoon(7000000, 10560.6)

        KE = []
        PE = []
        Etot = []
        for i in range(len(t)):
            Kinetici = (1 / 2) * Ms * (vx[i] ** 2 + vy[i] ** 2)
            KE.append(Kinetici)
            Potentiali = Mpotential(x[i], y[i])
            PE.append(Potentiali)
            Etot.append(Kinetici + Potentiali)

        print('Maximum Total Energy is', max(Etot))
        print('Minimum Total Energy is', min(Etot))
        print('The variation in total energy is', max(Etot) - min(Etot))
        print('Total time is', t[-1])
        print()


        # Earth and Moon Figures
        earth = plt.Circle((0, 0), Re, fill=True, color='darkslateblue')
        moon = plt.Circle((xm, 0), Rm, fill=True, color='grey')

        # Dotted line along y=0
        xaxisline = np.arange(xm, 0, 10 ** 2)
        yaxisline = np.zeros(len(xaxisline))

        fig1, orbit = plt.subplots(1)
        plt.plot(x, y, linestyle='--', linewidth='2', color='red')
        plt.plot(xaxisline, yaxisline, linestyle='--', color='black', linewidth='0.2')
        plt.axis('equal')
        plt.title('Orbit for body around the earth and moon')
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
        fig1 = plt.gcf()
        orbit = fig1.gca()
        orbit.add_patch(earth)
        orbit.add_patch(moon)

        fig2, earthenergy = plt.subplots(1)
        earthenergy.plot(t, KE, color='darkred')
        earthenergy.plot(t, PE, color='darkblue')
        earthenergy.plot(t, Etot, color='darkgreen')
        plt.ylabel('Energy (J)')
        plt.xlabel('Time (s)')
        earthenergy.legend(['Kinetic Energy', 'Potential Energy', 'Total Energy'], loc='best')

        plt.show()


    elif OrbitOption == 'd':
        """Earth Crash"""

        [x, y, vx, vy, t] = Earth(8000000, 6600)

        # Earth Figure
        earth = plt.Circle((0, 0), Re, fill=True, color='darkslateblue')

        KE = []
        PE = []
        Etot = []
        for i in range(len(t)):
            Kinetici = (1 / 2) * Ms * (vx[i] ** 2 + vy[i] ** 2)
            KE.append(Kinetici)
            Potentiali = potential(x[i], y[i])
            PE.append(Potentiali)
            Etot.append(Kinetici + Potentiali)

        fig1, orbit = plt.subplots(1)
        plt.plot(x, y, linestyle='--', linewidth='2', color='red')
        plt.axis('equal')
        plt.title('Rocket crash in to the Earth')
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
        fig1 = plt.gcf()
        orbit = fig1.gca()
        orbit.add_patch(earth)

        fig2, earthenergy = plt.subplots(1)
        earthenergy.plot(t, KE, color='darkred')
        earthenergy.plot(t, PE, color='darkblue')
        earthenergy.plot(t, Etot, color='darkgreen')
        plt.ylabel('Energy (J)')
        plt.xlabel('Time (s)')
        earthenergy.legend(['Kinetic Energy', 'Potential Energy', 'Total Energy'], loc='best')

        plt.show()



    elif  OrbitOption == 'e':
        """Moon Crash"""

        [x, y, vx, vy, t] = EarthMoon(7000000, 10567)

        KE = []
        PE = []
        Etot = []
        for i in range(len(t)):
            Kinetici = (1 / 2) * Ms * (vx[i] ** 2 + vy[i] ** 2)
            KE.append(Kinetici)
            Potentiali = Mpotential(x[i], y[i])
            PE.append(Potentiali)
            Etot.append(Kinetici + Potentiali)

        # Earth and Moon Figures
        earth = plt.Circle((0, 0), Re, fill=True, color='darkslateblue')
        moon = plt.Circle((xm, 0), Rm, fill=True, color='grey')

        # Dotted line along y=0
        xaxisline = np.arange(xm, 0, 10 ** 2)
        yaxisline = np.zeros(len(xaxisline))

        fig1, orbit = plt.subplots(1)
        plt.plot(x, y, linestyle='--', linewidth='2', color='red')
        plt.plot(xaxisline, yaxisline, linestyle='--', color='black', linewidth='0.2')
        plt.axis('equal')
        plt.title('Rocket crash in to the Moon')
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
        fig1 = plt.gcf()
        orbit = fig1.gca()
        orbit.add_patch(earth)
        orbit.add_patch(moon)

        fig2, earthenergy = plt.subplots(1)
        earthenergy.plot(t, KE, color='darkred')
        earthenergy.plot(t, PE, color='darkblue')
        earthenergy.plot(t, Etot, color='darkgreen')
        plt.ylabel('Energy (J)')
        plt.xlabel('Time (s)')
        earthenergy.legend(['Kinetic Energy', 'Potential Energy', 'Total Energy'], loc='best')

        plt.show()





    elif OrbitOption == 'q':
        print('Exiting Program')
        break

    else:
        print('Invalid input, returning to main menu\n')