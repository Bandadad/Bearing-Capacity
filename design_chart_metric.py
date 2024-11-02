import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fox_embedment as fox
import influence_factor
import weighted_modulus as _mod
import math
import yaml

# Functions

def create_dataframe(B_values, qu_values, q_values, three_quarter_q_values, half_q_values, quarter_q_values, modulus_values, I_s_values, I_f_values):
    data = {
        'B (m)': B_values,
        'q_allow (kPa)': qu_values,
        '25 mm (kPa)': q_values,
        '19 mm (kPa)': three_quarter_q_values,
        '12 mm (kPa)': half_q_values,
        '6 mm (kPa)': quarter_q_values,
        'Modulus (kPa)': modulus_values,
        'I_s': I_s_values,
        'I_f': I_f_values
    }
    df = pd.DataFrame(data)
    return df

def plot_and_save_results(B_values, qu_values, q_values, three_quarter_q_values, half_q_values, quarter_q_values, modulus_values, I_s_values, I_f_values, D, input_info, Shape, settlement_location):
    plt.figure(figsize=(10, 6))
    plt.plot(B_values, qu_values, color='k', linewidth=2, label=r'Strength Limit ($\phi$=0.45)')
    plt.plot(B_values, q_values, color='k', linestyle='dashed', label='Service Limit = 25 mm')
    plt.plot(B_values, three_quarter_q_values, color='k', linestyle='dashdot', label='Service Limit = 19 mm')
    plt.plot(B_values, half_q_values, color='k', linestyle='dotted', label='Service Limit = 12 mm')
    plt.plot(B_values, quarter_q_values, color='k', linestyle=(0, (3, 5, 1, 5, 1, 5)), label='Service Limit = 6 mm')
    plt.xlabel('Footing Width (m)', fontsize=14, weight='bold')
    plt.ylabel('Factored Net Bearing Resistance (kPa)', fontsize=14, weight='bold')
    plt.xlim(0, 10)
    plt.ylim(0, 600)
    details_text = f"Footing Depth: {D} m\nL/D Ratio: {Shape}\nCalculation Location: {settlement_location.capitalize()}"
    plt.annotate(details_text, xy=(0.5, 1.02), xycoords='axes fraction', ha='center', fontsize=10, bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))
    plt.legend(loc='upper right', fontsize=10, framealpha=1)
    plt.grid(which='major', linestyle='dashed', alpha=0.75)
    plt.xticks(np.arange(0, 12, 2))
    plt.yticks(np.arange(0, 650, 50))
    plt.tick_params(axis='both', labelsize=14)
    plt.minorticks_on()
    plt.tick_params(axis='y', which='minor', length=4)
    plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(10))
    plt.savefig("figure.png", dpi=1200, format='png', bbox_inches='tight')
    plt.show()

    df = create_dataframe(B_values, qu_values, q_values, three_quarter_q_values, half_q_values, quarter_q_values, modulus_values, I_s_values, I_f_values)

    # Save the results dataframe as an Excel file
    output_excel = 'output_results.xlsx'
    df.to_excel(output_excel, index=False)

    # Save the input information as a separate Excel file
    input_excel = 'input_information.xlsx'
    input_info_table = pd.DataFrame(input_info.items(), columns=['Variable', 'Value'])
    input_info_table.to_excel(input_excel, index=False)

    # Print the results dataframe and input information
    print("Output Results:")
    print(df)
    print("\nInput Information:")
    print(input_info_table)

    return output_excel, input_excel

def calculate_factors(phi):
    rad_phi = math.radians(phi)
    Nq = math.exp(math.pi * math.tan(rad_phi)) * (math.tan(math.pi / 4 + rad_phi / 2)) ** 2
    Nc = (Nq - 1) / math.tan(rad_phi)
    Ng = 2 * (Nq + 1) * math.tan(rad_phi)
    return Nc, Nq, Ng

def calculate_shape_factors(B, L, phi, Nq, Nc):
    rad_phi = math.radians(phi)
    Sc = 1 + (Nq / Nc) * (B / L)
    Sq = 1 + (B / L) * math.tan(rad_phi)
    Sg = max(0.6, 1 - 0.4 * (B / L))
    return Sc, Sq, Sg

def calculate_depth_factors(D, B, phi):
    rad_phi = math.radians(phi)
    ratio = D / B
    if ratio <= 1:
        k = np.minimum(np.arctan(ratio), D / B)
    else:
        k = np.arctan(ratio)
    Dc = 1 + 0.4 * k
    Dq = 1 + 2 * math.tan(rad_phi) * np.power(1 - math.sin(rad_phi), 2) * k
    Dg = 1
    return Dc, Dq, Dg

def calculate_bearing_capacity(B, D, L, gamma_soil, gamma_backfill, phi, c):
    Nc, Nq, Ng = calculate_factors(phi)
    Sc, Sq, Sg = calculate_shape_factors(B, L, phi, Nq, Nc)
    Dc, Dq, Dg = calculate_depth_factors(D, B, phi)
    qu = c * Nc * Sc * Dc + gamma_backfill * D * Nq * Sq * Dq + 0.5 * gamma_soil * B * Ng * Sg * Dg
    return qu * 0.45

def design_chart(Z, D, Shape, settlement_location, m, gamma_soil, gamma_backfill, phi, c, nu, modulus_file=False, file=None, Es=None):
    # Initialize arrays
    B_values = np.arange(0.25, 12, 0.25)  # foundation width in meters
    q_values = []
    qu_values = []
    three_quarter_q_values = []
    half_q_values = []
    quarter_q_values = []
    modulus_values = []
    I_s_values = []
    I_f_values = []

    for B in B_values:
        L = Shape * B  # This can be changed to a constant if L is a set dimension
        if m == 4:
            B_prime = B / 2
            L_prime = L / 2
        elif m == 2:
            B_prime = B / 2
            L_prime = L
        else:
            B_prime = B 
            L_prime = L
        H = min(Z, 5 * B)
        I_s = influence_factor.calculate_Is(B_prime, L_prime, H, nu)
        I_f = fox.fox_factor(B, L, D, nu)
        # Append values to respective lists
        I_s_values.append(I_s)
        I_f_values.append(I_f)

        qu = calculate_bearing_capacity(B, D, L, gamma_soil, gamma_backfill, phi, c)
        modulus = _mod.compute_weighted_average(file, H) if modulus_file else Es
        modulus_values.append(modulus)
        if B == 0:
            q_values.append(0)  # To avoid division by zero
        else:
            q = (0.0254) / (B_prime * ((1 - nu ** 2) / modulus) * m * I_s * I_f)
            q_values.append(q)
            three_quarter_q_values.append(q * 0.75)
            half_q_values.append(q * 0.5)
            quarter_q_values.append(q * 0.25)
        qu_values.append(qu)

    # Prepare input information
    input_info = {
        'Z': Z,
        'D': D,
        'Shape': Shape,
        'Settlement_location': settlement_location,
        'gamma_soil': gamma_soil,
        'gamma_backfill': gamma_backfill,
        'phi': phi,
        'c': c,
        'nu': nu,
        'modulus_file': modulus_file,
        'file': file,
        'Es': Es
    }

    # Call the plot_and_save_results function
    output_excel, input_excel = plot_and_save_results(B_values, qu_values, q_values, three_quarter_q_values, half_q_values, quarter_q_values, modulus_values, I_s_values, I_f_values, D, input_info, Shape, settlement_location)

    return output_excel, input_excel

def main():
    # Load configuration from YAML file with UTF-8 encoding
    with open('settle_config_metric.yaml', 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)

    # Extract parameters from the configuration
    D = config['settlement_parameters']['D']
    Shape = config['settlement_parameters']['shape']
    settlement_location = config['settlement_parameters'].get('settlement_location', 'center')

    gamma_backfill = config['soil_properties']['gamma_backfill']
    gamma_soil = config['soil_properties']['gamma_soil']
    Z = config['soil_properties']['Z']
    phi = config['soil_properties']['phi']
    c = config['soil_properties']['c']
    nu = config['soil_properties']['nu']

    modulus_file = config['modulus_options']['modulus_file']
    file = config['modulus_options']['file']
    Es = config['modulus_options'].get('Es')

    if settlement_location == 'center':
        m = 4
    elif settlement_location == 'edge':
        m = 2
    elif settlement_location == 'corner':
        m = 1
    else:
        raise ValueError("Invalid settlement location specified in the configuration file.")    

    # Run the settlement calculation function with loaded parameters
    output_excel, input_excel = design_chart(Z, D, Shape, settlement_location, m, gamma_soil, gamma_backfill, phi, c, nu, modulus_file, file, Es)

    print("\nOutput results saved as:", output_excel)
    print("Input information saved as:", input_excel)


if __name__ == "__main__":
    main()