import configparser



def create_summary(filebasename):

    summary = []


    config = configparser.ConfigParser()
    with open(filebasename+".ini", "r") as f:
        data = f.readlines()[1:] # remove the first line as fields without a section are not allowed
        config.read_string("".join(data))

    sim = config["Simulation"]
    summary.append(r"Number of simulation steps: \num{"+sim['numTimeSteps']+"}")
    summary.append(r"Length of a simulation step: $\Delta t = \SI{"+sim["dt"]+r"}{\second}$")
    summary.append(r"Number of real particles per simulated particle: $F_N = \num{"+sim["F_N"]+"}$")
    summary.append(r"Size of the simulation area: \SI{"+sim["V_c"]+r"}{\metre^2}")

    num_particles = 0
    with open(filebasename+".sim", "r") as f:
        num_rows, num_particle_species = [int(x) for x in f.readline().split(" ")]
        sim_steps, particles = [], []
        for l in f.readlines():
            d = l.split(" ")
            sim_steps.append(int(d[0]))
            particles.append([int(x) for x in d[1:1+num_particle_species]])
    for i in range(num_rows):
        if i < num_rows-1:
            sim_steps_in_section = sim_steps[i+1] - sim_steps[i]
        else:
            sim_steps_in_section = int(sim["numTimeSteps"])+1-sim_steps[i]
        for j in range(num_particle_species):
            num_particles += sim_steps_in_section * particles[i][j]

    with open(filebasename+".particles", "r") as f:
        n = int(f.readline())

    summary.append(r"Total number of simulated particles: $\num{"+str(num_particles+n)+"}$")



    qt = config["Quadtree"]
    r = int(qt["numStatisticsCellRows"])
    c = int(qt["numStatisticsCellColumns"])
    summary.append(r"Number of sampling cells: \num{"+str(r*c)+"}"+f" (${r}$ rows x ${c}$ columns)")
    summary.append(r"Quadtree parameters: $T_\text{upper} = "+qt["splitThreshold"]+r", T_\text{lower} = "+qt["mergeThreshold"]+"$")

    MM = config["MolecularModel"]
    n = int(MM["numParticleSpecies"])
    model = "VHS" if MM["collisionModel"] == "2" else "HS"
    summary.append(f"Molecular model: {model}" + (r" with $\nu = \num{"+MM["nu"]+"}") if model=="VHS" else "")
    summary.append(f"Number of particle species: {n}")
    m = MM["m"].split(",")
    d_ref = MM["d_ref"].split(",")
    T_ref = MM["T_ref"].split(",")
    for i in range(n):
        if model == "VHS":
            p = r"$m = \SI{"+m[i]+r"}{\kg}, d_\text{ref} = \SI{"+d_ref[i]+r"}{\metre}, T_\text{ref} = \SI{"+T_ref[i]+r"}{\kelvin}$"
        elif model == "HS":
            p = r"$m = \SI{" + m[i] + r"}{\kg}, d = \SI{" + d_ref[i] + r"}{\kg}$"
        summary.append(f"Particle species {i+1}: {p}")

    s = config["Surface"]
    surface_reflection_model = "diffuse" if int(s["surfaceCollisionModel"]) == 2 else "specular"
    surface = f"Surface reflection model: {surface_reflection_model}"
    if surface_reflection_model == "diffuse":
        surface += r", surface temperature: "+r"\SI{"+s['surfaceTemperature']+r"}{\kelvin}"
    summary.append(surface)

    return "\\\\\n".join(summary)


if __name__ == "__main__":
    filebasename = "../data/apollo_cm/cm7"
    filebasename = "../data/cylinder/c3"
    # filebasename = "../data/shear_flow/s1"
    # filebasename = "../data/temperature_exchange/T2"

    print(create_summary(filebasename))


