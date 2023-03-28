import subprocess, math

__all__ = [ "PrimaryFlow", "RunCEA" ]

def RunCEA(fileName, folderPath, nozzle):
    ''' 
        Combustion Analysis
        This function computes chemical equilibrium compositions for rockets using NASA CEA 
        It runs the executable file from the command prompt. 
    '''
    # Path to the executable file
    exe_path = "NASA_CEA\\FCEA2.exe" # C:\\Users\\Noreen\\OneDrive - University of Calgary\\Documents\\GitHub\\Research\\New folder\

    # Running the executable file and passing the argument as the user input
    p = subprocess.Popen(exe_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    input = fileName + "\n"
    stdout, stderr = p.communicate(input=input.encode())

    # Print the output from the executable file
    # print(stdout.decode(), stderr.decode())

    # Read output files using Python
    parameter_names = ['p', 't', 'rho', 'cp', 'gam', 'son', 'pip', 'mach', 'aeat', 'ivac', 'isp', 'mw']

    def parseFile(filePath):
        with open(filePath, 'r') as f:
            lines = f.readlines()
        
        header = lines[0].strip().split()
        parameter_indices = {parameter_name: header.index(parameter_name) for parameter_name in parameter_names}
        
        chamber_values = [float(x) for x in lines[1].strip().split()]
        chamber_parameter_values = [chamber_values[index-1] for index in parameter_indices.values()]
        
        throat_values = [float(x) for x in lines[2].strip().split()]
        throat_parameter_values = [throat_values[index-1] for index in parameter_indices.values()]
        
        exit_values = [float(x) for x in lines[3].strip().split()]
        exit_parameter_values = [exit_values[index-1] for index in parameter_indices.values()]
        return chamber_parameter_values, throat_parameter_values, exit_parameter_values

    filePath = (folderPath + "\\" + fileName + ".plt")
    chamber_parameter_values, throat_parameter_values, exit_parameter_values = parseFile(filePath)

    chamber = Parameter(*chamber_parameter_values)
    throat = Parameter(*throat_parameter_values)
    exit = Parameter(*exit_parameter_values)

    exhaust = Exhaust(chamber, throat, exit)

    exhaust.chamber.pressure = exhaust.chamber.pressure*1E5 # changing from bar to Pa
    exhaust.throat.pressure = exhaust.throat.pressure*1E5 # changing pressure units from bar to Pa

    exhaust.R = exhaust.throat.specific_heat_at_constant_pressure*(exhaust.throat.gamma-1)/exhaust.throat.gamma*1000  # J/kgK - specific gas constant of the exhaust
    exhaust.V_throat = exhaust.throat.sonicVelocity*exhaust.throat.mach # m/s - velocity of the exhaust flow at the nozzle throat
    exhaust.m_primary = exhaust.throat.density*exhaust.V_throat*nozzle.throatArea # kg/s - primary (exhaust) mass flow rate 

    # Exhaust Stagnation Properties
    exhaust.T0 = exhaust.throat.temperature*(1+(exhaust.throat.gamma-1)/2) # K - stagnation temperature of the exhaust
    exhaust.P0 = exhaust.throat.pressure*(exhaust.T0/exhaust.throat.temperature)**(exhaust.throat.gamma/(exhaust.throat.gamma-1)) # Pa - stagnation pressure of the exhaust
    exhaust.rho0 = exhaust.throat.density*(exhaust.P0/exhaust.throat.pressure)**(1/exhaust.throat.gamma) # kg/m3 - stagnation density of the exhaust
    return exhaust

class PrimaryFlow(): 
    ### Init Functions ### 
    def __init__(self, componentDictReader, rocket, stage):
        self.folderPath = componentDictReader.tryGetString("PrimaryFlow.folderPath", defaultValue=None)
        self.fileName = componentDictReader.tryGetString("PrimaryFlow.fileName", defaultValue=None)
        self.oxToFuelRatio = componentDictReader.tryGetFloat("PrimaryFlow.oxToFuelRatio", defaultValue=None)

class Parameter:
    def __init__(self, pressure, temperature, density, specific_heat_at_constant_pressure, gamma, sonicVelocity, pressure_ratio, mach, ratio_of_exit_area_to_throat_area, isp_vac, specific_impulse, molecular_weight):
        self.pressure = pressure # bar
        self.temperature = temperature # Kelvin
        self.density = density # kg/m^3
        self.specific_heat_at_constant_pressure = specific_heat_at_constant_pressure # kJ/kgK
        self.gamma = gamma # unitless
        self.sonicVelocity = sonicVelocity # m/s
        self.pressure_ratio = pressure_ratio # unitless
        self.mach = mach # unitless 
        self.ratio_of_exit_area_to_throat_area = ratio_of_exit_area_to_throat_area # unitless 
        self.isp_vac = isp_vac  # m/s 
        self.specific_impulse = specific_impulse # m/s
        self.Q = None # volumetric flow rate m3/s

class Exhaust:
    def __init__(self, chamber, throat, exit):
        self.chamber = chamber
        self.throat = throat
        self.exit = exit

class Nozzle():
    '''Defines the diverging section of the nozzle geometry'''

    #### Init Functions ####
    def __init__(self, componentDictReader, rocket, stage):
        self.angle = componentDictReader.tryGetFloat("Nozzle.angle", None) # degrees - full angle of the diverging nozzle section
        self.areaRatio = componentDictReader.tryGetFloat("Nozzle.areaRatio", None) # (non-dimensional) - area ratio of the nozzle, Aexit/Athroat
        self.throatDia = componentDictReader.tryGetFloat("Nozzle.throatDia", None) # m - nozzle throat diameter

        self._precomputeGeometry()

    def _precomputeGeometry(self):
        # Calculate dependent geometry
        self.exitDia = self.throatDia*math.sqrt(self.areaRatio) # m - nozzle exit diameter
        self.length = (self.exitDia - self.throatDia)/2/math.tan(math.radians(self.angle)/2) # m - nozzle diverging section length
        self.throatArea = self.throatDia**2*math.pi/4 # m2 - nozzle throat area
        self.exitArea = self.areaRatio*self.throatArea
