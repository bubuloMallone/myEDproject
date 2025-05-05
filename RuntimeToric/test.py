# Python script to compute column sums from a file

# Specify the file name
file_name = "test.txt"

# Initialize sums
sum_first_column = 0.0
sum_second_column = 0.0

try:
    # Open the file for reading
    with open(file_name, 'r') as file:
        # Process each line in the file
        for line in file:
            # Split the line into two columns
            values = line.strip().split()
            
            # Ensure the line contains two columns
            if len(values) == 2:
                # Convert columns to floats and add to the sums
                sum_first_column += float(values[0])
                sum_second_column += float(values[1])
except FileNotFoundError:
    print(f"File '{file_name}' not found.")
except ValueError:
    print("Encountered a value that could not be converted to float.")
else:
    # Print the results
    print(f"Sum of the first column: {sum_first_column}")
    print(f"Sum of the second column: {sum_second_column}")
