import csv
import itertools


def process_text_file(input_file):
    data = []

    with open(input_file, 'r') as file:
        lines = file.readlines()

        # Flag to determine if we are inside the data section
        inside_data = False

        for line in lines:
            if inside_data:
                # Check if the line contains data (not dashes or empty)
                if line.strip() and not line.startswith('-'):
                    # Split the line into columns
                    columns = line.split()

                    # Extract relevant information
                    alg, C, Cn, Crob, CrobScaled = columns[0], columns[1], columns[2], columns[3], columns[4]

                    # Append the data as a tuple
                    data.append((alg, C, Cn, Crob, CrobScaled))

            # Check if the line starts the data section
            if line.startswith('------------'):
                inside_data = not inside_data

    return data


def write_to_csv(data, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write header
        csv_writer.writerow(['alg', 'C', 'Cn', 'Crob', 'CrobScaled'])

        # Write data
        csv_writer.writerows(data)


if __name__ == "__main__":
    input_file = "Bridgeness/robustness.txt"
    # Process the text file and get data
    data = process_text_file(input_file)

    # Group data by 'alg' and write to separate CSV files
    for alg, group in itertools.groupby(data, key=lambda x: x[0]):
        output_file = f"Bridgeness/{alg}Robustness.csv"
        write_to_csv(list(group), output_file)

    print("CSV files generated successfully.")


