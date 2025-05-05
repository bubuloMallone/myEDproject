import math

def generate_time_evolution_file(U_max, U_min, total_time, Dt, N_timeSteps):
    """
    Generate a time evolution protocol file.

    Parameters:
    U_max (float): Maximum value of U at the start of the ramp.
    U_min (float): Minimum value of U at the end of the ramp.
    total_time (float): Total duration of the ramp in the same time units as dt.
    dt (float): Time step for each entry in the protocol.
    """
    num_steps = round(total_time / Dt)
    dt = Dt
    time = Dt
    DU = (U_max - U_min)/ num_steps  
    U_values = [U_max - ((x+1) * DU) for x in range(num_steps)]

    # Calculate the number of decimal places to represent Dt without scientific notation
    if Dt != 0:
        decimal_places = max(-int(math.floor(math.log10(abs(Dt)))), 0)
    else:
        decimal_places = 0
    # Format Dt with dynamic decimal places
    Dt_format = f'{{:.{decimal_places}f}}'
    Dt_str = Dt_format.format(Dt)

    # Filename includes total_time and dt
    filename = f'TimeEvoFile_Umax.{U_max}.t.{total_time}.Dt.{Dt_str}.txt'
    
    with open(filename, 'w') as file:
        for i, U in enumerate(U_values):
            
            # Write time entry
            file.write(f"time \n")
            # num 
            file.write(f"{i} \n")
            # time value
            file.write(f"{time} \n")

            # Write dt entry
            file.write(f"dt\n")
            # num
            file.write(f"{i} \n")
            # dt value
            file.write(f"{dt} \n")

            # Write U factor entry
            file.write(f"H1a \n")
            # num
            file.write(f"{i} \n")
            # U value
            file.write(f"{U} \n")

    print(f"File '{filename}' has been generated with {num_steps} entries.")

# Example usage
U_max = 1.0
U_min = 0.0
total_time_list = [0.001]  # Total time of the ramp 
Dt = 0.000001        # Time step
N_timeSteps = 50
for total_time in total_time_list:
    generate_time_evolution_file(U_max, U_min, total_time, Dt, N_timeSteps)
