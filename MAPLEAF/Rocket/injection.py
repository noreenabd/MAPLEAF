import math, numpy as np, scipy.optimize as so
from CoolProp.CoolProp import PropsSI

from MAPLEAF.Motion import Vector

__all__ = [ "Injection" , "computeForces"]

class Injection(): 
    def __init__(self, componentDictReader, rocket, stage):
        self.injectant = componentDictReader.getString("Injection.injectant")
        self.diameter = componentDictReader.getFloat("Injection.diameter")
        self.mach = componentDictReader.getFloat("Injection.mach") # Mach of the injection
        self.stagPress = componentDictReader.getFloat("Injection.stagPress")
        self.Cd = componentDictReader.getFloat("Injection.Cd")
        self.wswm = componentDictReader.getFloat("Injection.wswm")
        self.location = componentDictReader.getFloat("Injection.location")
        self._precomputeGeometry()
    
    def _precomputeGeometry(self):
        self.area = self.diameter**2*math.pi/4 #m^2, area of secondary injection

def computeForces(nozzle, exhaust, injection, rocket):
    environment = rocket.environment.getAirProperties(rocket.rigidBody.state.position, 0)

    def GetExhaustPropx(Nozzle, exhaust, x):
        # Obtains properties at a point in the nozzle using stagnation properties based on isentropic flow conditions
        # x is the distance in meters from the nozzle throat
        D = lambda x: (Nozzle.exitDia - Nozzle.throatDia)*x/Nozzle.length + Nozzle.throatDia # 

        A1 = math.pi/4*Nozzle.throatDia**2
        M1 = exhaust.throat.mach
        
        class Prop: 
            def __init__(self): 
                self.T = None
                self.P = None 
                self.vel = None 
                self.M = None 
                self.c = None 
                self.rho = None
        
        prop = Prop()

        Mach = lambda M2: A2*M2/(1+(exhaust.throat.gamma-1)/2*M2**2)**((exhaust.throat.gamma+1)/2/(exhaust.throat.gamma-1)) - A1*M1/(1+(exhaust.throat.gamma-1)/2*M1**2)**((exhaust.throat.gamma+1)/2/(exhaust.throat.gamma-1))
        # Initialize arrays 
        prop.M = np.array([])
        prop.T = np.array([])
        prop.P = np.array([])
        prop.rho = np.array([])
        prop.c = np.array([])
        prop.vel = np.array([])

        for i in range(len(x)):
            Dia = D(x[i])
            A2 = math.pi/4*Dia**2
            prop.M = np.append(prop.M, so.brentq(Mach, 1, 30))
            prop.T = np.append(prop.T, exhaust.T0/(1+(exhaust.throat.gamma-1)/2*prop.M[i]**2))
            prop.P = np.append(prop.P, exhaust.P0/(exhaust.T0/prop.T[i])**(exhaust.throat.gamma/(exhaust.throat.gamma-1)))
            prop.rho = np.append(prop.rho, exhaust.rho0/(exhaust.P0/prop.P[i])**(1/exhaust.throat.gamma))
            prop.c = np.append(prop.c, math.sqrt(exhaust.throat.gamma*exhaust.R*prop.T[i]))
            prop.vel = np.append(prop.vel, prop.c[i]*prop.M[i])
        
        return prop
    # Computes side forces, efficiency, and other important parameters 
    # Calculate the specific heat ratio (gamma) at the injection pressure and ~room temp
    injection.cP = PropsSI('CPMASS', 'P', injection.stagPress, 'T', 300, injection.injectant) 
    injection.cV = PropsSI('CVMASS', 'P', injection.stagPress, 'T', 300, injection.injectant) 
    injection.gamma = injection.cP/injection.cV

    # Assume that the injection Mach is 1
    injection.mach = 1 # Mach of the injection

    # Calculate molecular weight (MW)
    injection.MW = 1000*PropsSI('M', injection.injectant) # kg/kmol Molar Mass

    # Calculate specific gas constant (R)
    injection.R = 8314/injection.MW # J/kgK

    # Using M, gamma, Pt, find the injection pressure 
    injection.press = injection.stagPress/(1+(injection.gamma-1)/2*injection.mach**2)**(injection.gamma/(injection.gamma-1))
    # print('pressure of air injection', injection.press)

    # # Define the secondary to primary mass flow rate ratio(s)
    # StoP_ratio = np.array([0.032, 0.039])#0.10 # ratio of secondary to primary mass flow rate

    # Initialize variables of interest
    efficiency = [] # degree/percent injection
    side_force = np.array([]) # newtons
    impulse_side_specific =np.array([]) # seconds, side-specific impulse
    F_axial = np.array([])
    force_ratio = np.array([])
    F = np.array([])
    force_ratio_Broadwell = np.array([])

    Cpmax = lambda M, gamma: 1/(1/2*gamma*M**2)*(((gamma+1)/2*M**2)**(gamma/(gamma-1))/(2*gamma/(gamma+1)*M**2-(gamma-1)/(gamma+1))**(1/(gamma-1))-1)

    m_secondary = injection.wswm*exhaust.m_primary
    print('m_secondary', m_secondary)
    injection.density = m_secondary**2/(injection.gamma*injection.press*injection.area**2*injection.mach**2) #kg/m^3, density of the injection
    injection.temp = injection.press/injection.density/injection.R #  Kelvin, temperature of the injection
    print('temperature of injection', injection.temp)
    injection.c = math.sqrt(injection.gamma*injection.R*injection.temp) # speed of sound of the injectant at the injeciton port
    print('injection sonic vel', injection.c)
    
    # Computing the axial force
    exit = GetExhaustPropx(nozzle, exhaust, [nozzle.length]) # get isentropic properties at nozzle exit
    # print('pressure at the exit', press_exit)
    # THIS IS WRONG FIX THIS. NEED TO REMEMBER THAT THE FLOW GOES OUT AT AN ANGLE!!! 
    F_axial = np.append(F_axial, exhaust.m_primary*exit.vel) #+(press_exit-injection.P_amb)*nozzle.exit_dia**2*math.pi/4) # newtons, compute axial force FIX THIS!!!
    # print('primary mass flow rate', exhaust.m_primary)
    print('Axial Force', F_axial)


    # COMPUTATION OF THE SEPARATION POINT 
    x_inj = injection.location*nozzle.length # meters
    x_sep_old = 0.9*x_inj # meters, initial guess
    error = 5 #initialize error
    while error>0.0001:
        # Compute the exhaust flow properties at the separation location 
        separation = GetExhaustPropx(nozzle, exhaust, [x_sep_old])
        # print('separation pressure', P_separation)
        # print('separation mach', M_separation)
        # print('gamma', exhaust.throat.gamma)

        # Compute the plateau pressure using the Schilling Criterion 
        P_plateau = separation.P/0.582/(1+(exhaust.throat.gamma-1)/2*separation.M**2)**(-0.1197*exhaust.throat.gamma/(exhaust.throat.gamma-1)) # Pa
        # print('P_plateau', P_plateau)
        # Compute the shock wave angle (beta) in degrees using Oblique Shock relations (John & Keith p. 612-613)
        beta = math.asin(math.sqrt((P_plateau/separation.P*(exhaust.throat.gamma+1)+(exhaust.throat.gamma-1))/(2*exhaust.throat.gamma*separation.M**2)))
        # print('beta (shock wave angle)', (180/math.pi)*beta, 'degrees')

        # Compute the deflection angle (delta) in degrees
        delta = 1/math.tan(beta)*((separation.M*math.sin(beta))**2-1)/((exhaust.throat.gamma+1)/2*separation.M**2-((separation.M*math.sin(beta))**2-1))
        delta = (180/math.pi)*delta
        # print('delta', delta, 'degrees')

        # Compute Cpmax
        Cpmax_value = Cpmax(separation.M, exhaust.throat.gamma) #unitless

        # Compute dynamic pressure at the separation location 
        q0 = 1/2*separation.rho*separation.vel**2 # Pa, dynamic pressure before the shock

        # Compute the penetration height 
        # denominator = P_plateau - injection.P_amb + 1/2*q0*Cpmax_value
        denominator = separation.P - exit.P + 1/2*q0*Cpmax_value
        # print('denom', denominator)
        numerator = injection.Cd*injection.gamma*injection.stagPress*(2/(injection.gamma+1))**(1/(injection.gamma-1))*(1/(injection.gamma**2-1)*(1-(environment.Pressure/injection.stagPress)**((injection.gamma-1)/injection.gamma)))**0.5
        # print('numerator', numerator)
        h = injection.diameter*np.sqrt(numerator/denominator)

        # Compute the new separation location 
        x_sep_new = x_inj + h - h/math.sin(delta*math.pi/180)
        error = abs(x_sep_new - x_sep_old)
        x_sep_old = x_sep_new
    print('deflection angle', delta)
    # # SEPARATION POINT COMPUTATION 
    # x_inj = location*nozzle.length # meters

    # x_sep_old = 0.9*x_inj # meters, initial guess
    # error = 5 #initialize error
    # while error>0.0001:
    #     M_separation, rho_separation, vel_separation, P_separation = GetExhaustPropx(nozzle, exhaust, x_sep_old)
    #     print('separation pressure',P_separation)
    #     output_vars, output_vals = calc_oblique_with_gam_const("M1", M_separation, "M2", M_separation*0.78, exhaust.throat.gamma, ["THETA", "BETA"], quiet=True)
    #     delta = output_vals[0]
    #     theta = output_vals[1]
    #     P_plateau = exhaust.P0/0.582/(1+(exhaust.throat.gamma-1)/2*M_separation**2)**(-0.1197*exhaust.throat.gamma/(exhaust.throat.gamma-1)) # Pa
    #     print('plateau pressure', P_plateau)
    #     Cpmax_value = Cpmax(M_separation, exhaust.throat.gamma) #unitless
    #     q0 = 1/2*rho_separation*vel_separation**2 # Pa, dynamic pressure before the shock

    #     # denominator = P_plateau - injection.P_amb + 1/2*q0*Cpmax_value
    #     denominator = P_plateau - press_exit + 1/2*q0*Cpmax_value
    #     # print('denom', denominator)
    #     numerator = injection.Cd*injection.gamma*injection.stagPress*(2/(injection.gamma+1))**(1/(injection.gamma-1))*(1/(injection.gamma**2-1)*(1-(injection.P_amb/injection.stagPress)**((injection.gamma-1)/injection.gamma)))**0.5
    #     # print('numerator', numerator)
    #     h = injection.diameter*np.sqrt(numerator/denominator)
    #     x_sep_new = x_inj - h*(1/math.tan(delta)-1)
    #     error = abs(x_sep_new - x_sep_old)
    #     x_sep_old = x_sep_new
    #     # print('x_sep', x_sep_new)
    
    print('h',h)
    print('x_sep', x_sep_new)
    # print('nozzle length', nozzle.length)
    # efficiency.append(math.degrees(delta)/(injection.wswm[i]*100)) # degree/percent injection
    # side_force = np.append(side_force, math.pi*h**2*((P_plateau-injection.P_amb)/2+0.25*q0*(Cpmax_value))) # newtons

    # side_force = np.append(side_force, math.pi*h**2*((separation.P-exit.P)/2+0.25*q0*(Cpmax_value))) # newtons
    # ^ this is completely wrong lol. This is the axial force on the control volume. not even the lateral force at all!! 

    # 
    Ls = x_inj - x_sep_new
    print('Ls', Ls)
    # Compute radius of curvature, Rc, for the hyperbola at y = 0 (z = 0)
    Rc = h*1.143*math.exp(0.54/(separation.M-1)**1.2)
    # Compute Delta, ▲, the distance from the shock to the stopping point
    Delta = h*0.143*math.exp(3.24/separation.M**2)
    # Create linearly spaced array of x positions from the separation point to the injection location. Each element is xk
    x = np.linspace(x_sep_new, x_inj, num = 50)
    # print('x', x)
    
    # Compute Rs, the separation distance from the point c
    Rs = Ls + h*(1.143*np.exp(0.54/(separation.M-1)**1.2)-0.143*np.exp(3.24/separation.M**2))
    print('Rs', Rs)
    x_c = x_sep_new + Rs
    # plt.plot(x_c,0,'o')

    # Compute the Mach angle in radians (characterizing the asymptote of the hyperbola after in the injection location)
    mu = math.asin(1/separation.M)

    # Compute and plot the separation line 
    b = -Rs/np.tan(mu)
    a = Rs/(np.tan(mu))**2
    k = x_c - Rs - a

    z = -b*np.sqrt((x-k)**2/a**2-1)
    # plt.plot(x,z, label='separation line 1')


    k = (x_sep_new*math.sqrt(1+(math.tan(mu))**2)-x_c)/(math.sqrt(1+(math.tan(mu))**2)-1)
    a = x_sep_new-k
    print('mu', mu)
    print('a', a)
    b = a*np.tan(mu)
    print('b', b)
    z = b*np.sqrt((x-k)**2/a**2-1)
    # plt.plot(x,z, label='separation line 2')


    # Compute and plot the asymptote
    print('k/nozzle length', k/nozzle.length)
    x_asymp = np.linspace(k,x_inj,50)
    y=b/a*(x_asymp-k)
    # plt.plot(x_asymp, y, '.')
    

    # Compute the radius at various x positions along the nozzle diverging section
    x_nozzle = np.linspace(0,nozzle.length, num = 50)
    r_nozzle = nozzle.throatDia/2+x_nozzle*math.tan(nozzle.angle/2*math.pi/180)
    # plt.plot(x_nozzle, r_nozzle, label = 'nozzle radius')
    r = nozzle.throatDia/2+x*math.tan(nozzle.angle/2*math.pi/180)

    # Compute the maximum phi angle
    phi_max = np.arcsin(z/r)
    phi_max[np.isnan(phi_max)] = np.pi/2

    phi_max = z/r
    # print('phi_max', phi_max*180/np.pi)
    # a = (x - x_sep_new)*np.tan(mu)
    # phi_max = a/r

    # print('phi_max', phi_max*180/np.pi)
    # print('y/r',y/r)

    # Compute and plot the detached shock 
    k_detach = ((x_inj-Delta)*np.sqrt(1+(math.tan(mu))**2)-x_c)/(np.sqrt(1+(math.tan(mu))**2)-1)
    x_detach = np.linspace(x_inj-Delta,x_inj+h, 50)
    a_detach = x_inj-Delta-k_detach 
    b_detach = a_detach*np.tan(mu)
    y_detach = b_detach*np.sqrt((x_detach-k_detach)**2/a_detach**2 - 1)
    # plt.plot(x_detach, y_detach, label='detached shock')
    
    # Compute and plot the asymptotes 
    x_asymp_detach = np.linspace(k_detach,x_inj,50)
    y_detach=b_detach/a_detach*(x_asymp_detach-k_detach)
    # plt.plot(x_asymp_detach, y_detach, '.')

    # Compute and plot the injection 
    x_circ = np.linspace(x_inj,x_inj+h, 50)
    y_circ = np.sqrt(h**2-(x_circ-(x_inj+h))**2)
    # fig = plt.plot(x_circ,y_circ, label='injection')

    # Configuring plot
    # plt.xlabel('x [mm]')
    # plt.ylabel('z [mm]')
    # ax = plt.gca()
    # ax.set_aspect('equal', 'box')

    # plt.legend()
    def xpos2ratio(x):
        return x/nozzle.length*100


    def ratio2xpos(x):
        return x*nozzle.length/100

    # secax = ax.secondary_xaxis('top', functions=(xpos2ratio, ratio2xpos))
    # secax.set_xlabel('Percent of Nozzle length [%]')
    # plt.show()

    # NEED TO COMPUTE PRESSURE DISTRIBUTION
    # phi_max = 70*math.pi/180 # temporary
    
    properties = GetExhaustPropx(nozzle, exhaust, x)
    # print('r', r)


    # F_yb is the normal force in the area in the hyperbola due to the force on the nozzle wall. The pressure is interpolated between the P_plateau and the isentropic freestream pressure. 
    # Compute the pressure at each index i,j 
    # using elliptical interpolation between the plateau pressure and the freestream pressure 
    # at each location x (index j) and angular position phi (index j)
    
    # Initialize empty 2D array to hold pressures p_b and elemental force F_yb
    phi = np.linspace(0, phi_max, num = 50) 
    p_b = np.zeros((len(x), len(phi)))
    delF_yb = np.zeros((len(x), len(phi)))
    delF_xb = np.zeros((len(x), len(phi)))
    
    for k in range(len(x)-1):
        phi = np.linspace(0, phi_max[k], num = 50) 
        Lrk = phi_max[k]*r[k]
        if k == 0: 
            continue
        for j in range(len(phi)-1): 
            if j == 0: 
                continue
            p_b[k][j] = (P_plateau-properties.P[k])*math.sqrt(1-(phi[j]*r[k]/Lrk)**2)
            delF_yb[k][j] = p_b[k][j]*(x[k+1]-x[k-1])/2*r[k]*(math.sin(phi[j+1])-math.sin(phi[j-1]))/2#*np.cos(phi[j])
            # if phi[j] > np.pi/2 and phi[j] < 3*np.pi/2:
            #     delF_yb[k][j] = -delF_yb[k][j]
            delF_xb[k][j] = p_b[k][j]*(r[k+1]-r[k-1])/2*(phi[j+1]-phi[j-1])*r[k]
    # np.savetxt('delF_yb.csv', delF_yb, delimiter=',')
    # print('delF_yb size', np.shape(delF_yb))
    # print('delF_yb', delF_yb)
    F_yb = 2*np.sum(delF_yb)
    F_xb = 2*np.sum(delF_xb)
    print('F_yb', F_yb)

    # F_ytuy is the normal force on the nozzle wall not disturbed by the injection (i.e. before separation). The normal force cancels due to symmetry: 
    F_ytuy = 0

    # F_ys is the normal force in the region of separation where the plateau pressure is constant for an annular injection with an angle of phi. 
    # For an orifice injection, this area is zero. As a result, F_ys goes to zero
    F_ys = 0

    # F_yj is the normal force due to the momentum of the injected jet 
    inj_plane_prop = GetExhaustPropx(nozzle, exhaust, [injection.location]) # get isentropic properties at injection location
    F_yj = m_secondary*injection.c + injection.area*(injection.press-inj_plane_prop.P)# N, lateral force due to the momentum of the injection
    print('F_yj', F_yj)

    # F_yd is the normal force that arises due to overexpanded nozzles
    F_yd = 0

    # F_yr is the normal force that arises due to underexpanded nozzles (without free separation). i.e. there is reattachement of the flow. 
    # The effects of reattachment are neglected in the case of circular injection orifice. 
    F_yr = 0

    F_y = F_ytuy + F_ys + F_yb + F_yj + F_yd + F_yr
    print('F_y', F_y)


    side_force = F_y
    impulse_side_specific = np.append(impulse_side_specific, side_force/m_secondary/9.81) # seconds, side-specific impulse
    # print('efficiency', efficiency,'°/%')
    # print('h/(dc^.5)', h/injection_1.diameter/math.sqrt(injection_1.Cd))
    print('Side Force', side_force)
    # print('Side-specific impulse', impulse_side_specific, 's')


    F_throat = np.pi*(nozzle.throatDia/2)**2*exhaust.throat.pressure*(1+exhaust.throat.gamma*exhaust.throat.mach**2)
    F_x = F_throat + F_xb 
    print('F_x', F_x)

    angle = math.atan(side_force/F_axial)*180/math.pi
    print('deflection angle', angle, 'degrees')
    force_ratio = np.append(force_ratio, side_force/F_axial) # ratio of side to axial force
    print('Fs/Fa', force_ratio*100 ,'%')
    # print('exit velocity', vel_exit)
    # print('exit Mach', M_exit)

    # BROADWELL BLASTWAVE THEORY 
    inj_plane_prop = GetExhaustPropx(nozzle, exhaust, [injection.location]) # get isentropic properties at injection location
    F_interaction = 0.1*inj_plane_prop.M*inj_plane_prop.vel*m_secondary
    # print('injectant mass flow rate', m_secondary)
    # print('Interaction Force', F_interaction)
    injection.velocity = m_secondary/injection.density/(math.pi*injection.diameter**2/4)
    F_jet_reaction = m_secondary*injection.velocity
    # print('Jet reaction force', F_jet_reaction)
    F = np.append(F, F_interaction + F_jet_reaction)
    force_ratio_Broadwell = np.append(force_ratio_Broadwell, F/F_axial) # ratio of side to axial force
    return Vector(0, F_y, F_axial) 
