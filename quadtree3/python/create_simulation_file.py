import numpy as np
from scipy.constants import Boltzmann as k_B



def write_simulation_from_sequence_dict(filename, sim_seq, sep=" "):

    num_steps = len(sim_seq)
    num_particle_classes = len(sim_seq[list(sim_seq.keys())[0]])

    file_content = []
    file_content.append(f"{num_steps}{sep}{num_particle_classes}")

    for step in sorted(sim_seq.keys()):
        if len(sim_seq[step]) != num_particle_classes:
            raise ValueError("Invalid format of the simulation sequence dictionary. "
                             "The number of entries at each step needs to be the same.")

        num_particles = []
        spawn_areas = []
        velocity_distributions = []
        for c in range(num_particle_classes):
            conf = sim_seq[step][c]
            if conf.get("particle_number",0) == 0:
                num_particles.append(0)
                spawn_areas.extend([0.]*4)
                velocity_distributions.extend([0.]*6)
            else:
                num_particles.append(conf["particle_number"])
                area = conf.get("spawn_area")
                if area is None:
                    raise ValueError("Invalid format of the simulation sequence dictionary. "
                                     f"No spawn area is given at step {step}, particle class {c+1}.")
                spawn_areas.extend(area)
                vel_dist = conf.get("velocity_distribution")
                if vel_dist is None:
                    raise ValueError("Invalid format of the simulation sequence dictionary. "
                                     f"No velocity distribution is given at step {step}, particle class {c+1}.")
                velocity_distributions.extend(vel_dist)

        line_list = [step]+num_particles+spawn_areas+velocity_distributions
        file_content.append(sep.join(str(x) for x in line_list))

    with open(filename, "w") as f:
        f.write("\n".join(file_content))

    print(f"Successfully written the simulation sequence to '{filename}'.")


def calculate_influx(number_density, influx_surface_area, velocity_distribution, dt):

    if np.all(np.array(influx_surface_area[2:]) == 0):
        raise ValueError("influx_surface_area can't be zero")
    elif influx_surface_area[2] == 0:
        A = influx_surface_area[3]
    elif influx_surface_area[3] == 0:
        A = influx_surface_area[2]
    else:
        raise ValueError("influx_surface_area needs to be either in x or in y direction")

    v = np.linalg.norm(velocity_distribution[:3])

    # number_density -> n = N * v*dt * A
    N = number_density/(v*dt*A)
    return {
        "particle_number": int(N),
        "spawn_area": influx_surface_area,
        "velocity_distribution": velocity_distribution
    }


if __name__ == "__main__":
    filename = "../data/vertical_line/vl2.sim"

    particle_masses = [
        4.6518e-27,  # mass of N2
        10,
    ]
    stream_temperature = 300 # Kelvin
    sigma = lambda T, m: np.sqrt(k_B * T/m)

    simulation_sequence = {
        1: [
            calculate_influx(number_density=1e1,
                             influx_surface_area=[0,0,0,15],
                             velocity_distribution=[1000,0,0]+[sigma(stream_temperature, particle_masses[0])]*3,
                             dt=5e-6),
            {
                "particle_number": 10,
                "spawn_area": [0,0,0,15],
                "velocity_distribution": [1000,0,0]+[sigma(stream_temperature, particle_masses[0])]*3,
            },
            # {
            #     "particle_number": 0,
            #     "spawn_area": [0,0,100, 100],
            #     "velocity_distribution": [100,100,100]+[sigma(stream_temperature, particle_masses[1])]*3,
            # }
        ],
        1000: [
            {
                "particle_number": 0,
            },
            {
                "particle_number": 0,
            }
        ]
    }

    write_simulation_from_sequence_dict(filename, simulation_sequence)