import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fox_embedment as fox
import influence_factor
import weighted_modulus as _mod 
import math
import tkinter as tk
from tkinter import ttk, messagebox


# Functions
def resource_path(relative_path):
    """ Get the absolute path to a resource"""
    try:
        # PyInstaller creates a temporary folder and stores path in _MEIPASS.
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

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

def plot_and_save_results(B_values, qu_values, q_values, three_quarter_q_values, half_q_values, quarter_q_values, modulus_values, I_s_values, 
                          I_f_values, D, input_info, Shape, settlement_location, max_footing_width, max_bearing_pressure, x_ticks, y_ticks, r_factor):
    plt.figure(figsize=(10, 6))
    plt.plot(B_values, qu_values, color='k', linewidth=2, label=rf'Strength Limit ($\phi$={r_factor})')
    plt.plot(B_values, q_values, color='k', linestyle='dashed', label='Service Limit = 1 inch')
    plt.plot(B_values, three_quarter_q_values, color='k', linestyle='dashdot', label='Service Limit = 3/4 inch')
    plt.plot(B_values, half_q_values, color='k', linestyle='dotted', label='Service Limit = 1/2 inch')
    plt.plot(B_values, quarter_q_values, color='k', linestyle=(0, (3, 5, 1, 5, 1, 5)), label='Service Limit = 1/4 inch')
    plt.xlabel('Footing Width (ft)', fontsize=14, weight='bold')
    plt.ylabel('Factored Net Bearing Resistance (ksf)', fontsize=14, weight='bold')
    plt.xlim(0, max_footing_width)
    plt.ylim(0, max_bearing_pressure)
    details_text = f"Footing Depth: {D} ft\nL/D Ratio: {Shape}\nCalculation Location: {settlement_location.capitalize()}"
    plt.annotate(details_text, xy=(0.5, 1.02), xycoords='axes fraction', ha='center', fontsize=10, bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))
    plt.legend(loc='upper right', fontsize=10, framealpha=1)
    plt.grid(which='major', linestyle='dashed', alpha=0.75)
    plt.xticks(np.arange(0, max_footing_width + x_ticks, x_ticks))
    plt.yticks(np.arange(0, max_bearing_pressure + y_ticks, y_ticks))
    plt.tick_params(axis='both', labelsize=14, which='both', length=4)
    plt.minorticks_on()
    #plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.2))
    plt.savefig("figure.png", dpi=1200, format='png', bbox_inches='tight')
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

def calculate_bearing_capacity(B, D, L, gamma_soil, gamma_backfill, phi, c, r_factor):
    Nc, Nq, Ng = calculate_factors(phi)
    Sc, Sq, Sg = calculate_shape_factors(B, L, phi, Nq, Nc)
    Dc, Dq, Dg = calculate_depth_factors(D, B, phi)
    qu = c * Nc * Sc * Dc + gamma_backfill * D * Nq * Sq * Dq + 0.5 * gamma_soil * B * Ng * Sg * Dg
    return qu * r_factor / 1000

def design_chart(Z, D, Shape, settlement_location, m, gamma_soil, gamma_backfill, phi, c, nu, max_footing_width, max_bearing_pressure, 
                 x_ticks, y_ticks, r_factor, modulus_file=False, file=None, Es=None):
    # Initialize arrays
    B_values = np.arange(0.25, max_footing_width + 0.25, 0.25)  # foundation width in feet
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

        qu = calculate_bearing_capacity(B, D, L, gamma_soil, gamma_backfill, phi, c, r_factor)
        modulus = _mod.compute_weighted_average(file, H) if modulus_file else Es
        modulus_values.append(modulus)
        if B == 0:
            q_values.append(0)  # To avoid division by zero
        else:
            q = (1/12) / (B_prime * ((1 - nu ** 2) / modulus) * m * I_s * I_f)
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
    output_excel, input_excel = plot_and_save_results(
        B_values, qu_values, q_values, three_quarter_q_values, half_q_values,
        quarter_q_values, modulus_values, I_s_values, I_f_values,
        D, input_info, Shape, settlement_location, max_footing_width, max_bearing_pressure,
        x_ticks, y_ticks, r_factor
    )
    
    print("\nOutput results saved as:", output_excel)
    print("Input information saved as:", input_excel)
    return output_excel, input_excel

# ======================
# TKINTER GUI CODE
# ======================

class DesignChartApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Spread Footing Design Chart - US Customary Units")
        self.root.iconbitmap(resource_path("footing.ico"))       
        self.inputs = {}
        self.create_widgets()

    def create_widgets(self):
        sections = {
            "Settlement Parameters": [
                ("Foundation Depth (ft):", "D"),
                ("Shape (L/B ratio):", "Shape"),
                ("Settlement Location (center/edge/corner):", "settlement_location"),
                ("Resistance Factor (For ASD, use 1 / FoS):", "r_factor")
            ],
            "Soil Properties": [
                ("Unit Weight of Backfill (pcf):", "gamma_backfill"),
                ("Unit Weight of Foundation Soil (pcf):", "gamma_soil"),
                ("Depth to Hard Layer (ft):", "Z"),
                ("Soil Friction Angle (degrees):", "phi"),
                ("Soil Cohesion (psf):", "c"),
                ("Poisson's Ratio:", "nu")
            ],
            "Modulus Options": [
                ("Use Modulus File (True/False):", "modulus_file"),
                ("Modulus File Name (if used):", "file"),
                ("Settlement Modulus (ksf, None if no file):", "Es")
            ],
            "Plot Options": [
                ("Max X-Axis Footing Width (ft):", "max_footing_width"),
                ("Max Y-Axis Bearing Pressure (ksf):", "max_bearing_pressure"),
                ("X-axis Tick Interval:", "x_ticks"),
                ("Y-axis Tick Interval:", "y_ticks")
            ]
        }

        row = 0
        for section, fields in sections.items():
            frame = ttk.LabelFrame(self.root, text=section, padding=10)
            frame.grid(row=row, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
            for i, (text, var) in enumerate(fields):
                label = tk.Label(frame, text=text, anchor="e", width=35)
                label.grid(row=i, column=0, padx=5, pady=2, sticky="e")
                entry = tk.Entry(frame, width=20)
                entry.grid(row=i, column=1, padx=5, pady=2, sticky="w")
                self.inputs[var] = entry
            row += 1

        submit_btn = tk.Button(self.root, text="Generate Chart", command=self.process_inputs)
        submit_btn.grid(row=row, column=0, columnspan=2, pady=10)

    def process_inputs(self):
        try:
            # Gather and convert inputs from the GUI entries
            D = float(self.inputs["D"].get())
            Shape = float(self.inputs["Shape"].get())
            settlement_location = self.inputs["settlement_location"].get().lower()
            r_factor = float(self.inputs["r_factor"].get())
            gamma_backfill = float(self.inputs["gamma_backfill"].get())
            gamma_soil = float(self.inputs["gamma_soil"].get())
            Z = float(self.inputs["Z"].get())
            phi = float(self.inputs["phi"].get())
            c = float(self.inputs["c"].get())
            nu = float(self.inputs["nu"].get())
            max_footing_width = float(self.inputs["max_footing_width"].get())
            max_bearing_pressure = float(self.inputs["max_bearing_pressure"].get())
            x_ticks = float(self.inputs["x_ticks"].get())
            y_ticks = float(self.inputs["y_ticks"].get())
            modulus_file = self.inputs["modulus_file"].get().strip().lower() == "true"
            file = self.inputs["file"].get() if modulus_file else None
            Es = float(self.inputs["Es"].get()) if not modulus_file else None
            
            # Map settlement location to m value
            if settlement_location == 'center':
                m = 4
            elif settlement_location == 'edge':
                m = 2
            elif settlement_location == 'corner':
                m = 1
            else:
                raise ValueError("Invalid settlement location. Choose 'center', 'edge', or 'corner'.")
            
            # Call the original design_chart function
            design_chart(Z, D, Shape, settlement_location, m, gamma_soil, gamma_backfill, phi, c, nu, 
                         max_footing_width, max_bearing_pressure, x_ticks, y_ticks, r_factor, modulus_file, file, Es)
            messagebox.showinfo("Success", "Design Chart Generated Successfully")
        except Exception as e:
            messagebox.showerror("Input Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = DesignChartApp(root)
    root.mainloop()

