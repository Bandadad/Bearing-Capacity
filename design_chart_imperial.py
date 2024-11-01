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
        'B (ft)': B_values,
        'q_allow (ksf)': qu_values,
        '1 inch (ksf)': q_values,
        '3/4 inch (ksf)': three_quarter_q_values,
        '1/2 inch (ksf)': half_q_values,
        '1/4 inch (ksf)': quarter_q_values,
        'Modulus (ksf)': modulus_values,
        'I_s': I_s_values,
        'I_f': I_f_values
    }
    df = pd.DataFrame(data)
    return df


def plot_and_save_results(B_values, qu_values, q_values, three_quarter_q_values, half_q_values, quarter_q_values, modulus_values, I_s_values, I_f_values, D, input_info):
    plt.figure(figsize=(10, 6))
    plt.plot(B_values, qu_values, color='k', linewidth=2, label=r'Strength Limit ($\phi$=0.45)')
    plt.plot(B_values, q_values, color='k', linestyle='dashed', label='Service Limit = 1 inch')
    plt.plot(B_values, three_quarter_q_values, color='k', linestyle='dashdot', label='Service Limit = 3/4 inch')
    plt.plot(B_values, half_q_values, color='k', linestyle='dotted', label='Service Limit = 1/2 inch')
    plt.plot(B_values, quarter_q_values, color='k', linestyle=(0, (3, 5, 1, 5, 1, 5)), label='Service Limit = 1/4 inch')
    plt.xlabel('Footing Width (ft)', fontsize=14)
    plt.ylabel('Factored Net Bearing Resistance (ksf)', fontsize=14)
    plt.xlim(0, 20)
    plt.ylim(0, 12)
    plt.title(f'Design Chart for Footing Depth={D} ft', fontsize=16)
    plt.legend(loc='upper right', fontsize=10, framealpha=1)
    plt.grid(which='both', linestyle='dashed', alpha=0.75)
    plt.xticks(np.arange(0, 22, 2))
    plt.yticks(np.arange(0, 14, 2))
    plt.tick_params(axis='both', labelsize=14)
    plt.savefig("figure.svg", dpi=1200, format='svg', bbox_inches='tight')
    plt.show()

    df = create_dataframe(B_values, qu_values, q_values, three_quarter_q_values, half_q_values,
                          quarter_q_values, modulus_values, I_s_values, I_f_values)

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
    else: k = np.arctan(ratio)
    Dc = 1 + 0.4 * k
    Dq = 1 + 2 * math.tan(rad_phi) * np.power(1-math.sin(rad_phi), 2) * k
    Dg = 1
    return Dc, Dq, Dg


def calculate_bearing_capacity(B, D, L, gamma_soil, gamma_backfill, phi, c):
    Nc, Nq, Ng = calculate_factors(phi)
    Sc, Sq, Sg = calculate_shape_factors(B, L, phi, Nq, Nc)
    Dc, Dq, Dg = calculate_depth_factors(D, B, phi)
    qu = c * Nc * Sc * Dc + gamma_backfill * D * Nq * Sq * Dq + 0.5 * gamma_soil * B * Ng * Sg * Dg
    return qu * 0.45 / 1000


def design_chart(Z, D, Shape, gamma_soil, gamma_backfill, phi, c, nu, modulus_file=False, file=None, Es=None):
    # Initialize arrays
    B_values = np.arange(0.25, 41, 0.25)  # foundation width in feet
    q_values = []
    qu_values = []
    three_quarter_q_values = []
    half_q_values = []
    quarter_q_values = []
    modulus_values = []
    I_s_values = []
    I_f_values = []

    for B in B_values:
        L = Shape * B # This can be changed to a constant if L is a set dimension
        H = min(Z, 5 * B)
        I_s = influence_factor.calculate_Is(B, L, H, nu)
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
            q = (1/12) / ((B / 2) * ((1 - nu ** 2) / modulus) * 4 * I_s * I_f)
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
    output_excel, input_excel = plot_and_save_results(B_values, qu_values, q_values, three_quarter_q_values, half_q_values,
                                                      quarter_q_values, modulus_values, I_s_values, I_f_values, D,
                                                      input_info)

    return output_excel, input_excel

def main():
    # Load configuration from YAML file with UTF-8 encoding
    with open('settle_config.yaml', 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)

    # Extract parameters from the configuration
    D = config['settlement_parameters']['D']
    Shape = config['settlement_parameters']['shape']

    gamma_backfill = config['soil_properties']['gamma_backfill']
    gamma_soil = config['soil_properties']['gamma_soil']
    Z = config['soil_properties']['Z']
    phi = config['soil_properties']['phi']
    c = config['soil_properties']['c']
    nu = config['soil_properties']['nu']

    modulus_file = config['modulus_options']['modulus_file']
    file = config['modulus_options']['file']
    Es = config['modulus_options'].get('Es')

    # Run the settlement calculation function with loaded parameters
    output_excel, input_excel = design_chart(Z, D, Shape, gamma_soil, gamma_backfill, phi, c, nu, modulus_file, file, Es)
    
    print("\nOutput results saved as:", output_excel)
    print("Input information saved as:", input_excel)

if __name__ == "__main__":
    main()
