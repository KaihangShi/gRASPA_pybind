# Path to the text file
file_path = 'output.txt'

# Function to extract sections from the text
def extract_section_values(file_path):
    sections_data = {
        "FINAL STAGE": {},
        "ENERGY DRIFT": {}
    }

    with open(file_path, 'r') as file:
        lines = file.readlines()

    current_section = None
    section_started = False

    for line in lines:
        # Check for section headers
        if "*** FINAL STAGE ***" in line:
            current_section = "FINAL STAGE"
            section_started = True
            continue
        elif "*** ENERGY DRIFT " in line:
            current_section = "ENERGY DRIFT"
            section_started = True
            continue

        # If we're in a relevant section, extract the values
        if section_started and current_section:
            if "VDW [Host-Host]" in line or \
               "VDW [Host-Guest]" in line or \
               "VDW [Guest-Guest]" in line or \
               "Real Coulomb [Host-Host]" in line or \
               "Real Coulomb [Host-Guest]" in line or \
               "Real Coulomb [Guest-Guest]" in line or \
               "Ewald [Host-Host]" in line or \
               "Ewald [Host-Guest]" in line or \
               "Ewald [Guest-Guest]" in line or \
               "DNN Energy" in line or \
               "Tail Correction Energy" in line or \
               "Total Energy" in line:
                
                # Extract key and value
                key_value = line.split(":")
                key = key_value[0].strip()
                value = float(key_value[1].split()[0].strip())
                
                # Store the value in the current section
                sections_data[current_section][key] = value
            
            # Stop if we reach the next section or an empty line
            if line.strip() == "":
                section_started = False

    return sections_data

# Call the function to extract values
final_data = extract_section_values(file_path)

# Display the extracted values
for section, values in final_data.items():
    print(f"\n{section} Section:")
    for key, value in values.items():
        print(f"{key}: {value}")

