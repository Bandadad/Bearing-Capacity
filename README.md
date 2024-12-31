# Bearing Capacity and Settlement of Shallow Foundations

This project produces a design chart for spread footing using elastic theory, based on Bowles' *Foundation Analysis and Design*, 5th Edition.

### Settlement Calculation

Settlement is estimated from the equation:

ΔH = (q * B' * (1 - ν²) / Es) * m * Is * If

### Project Structure

The repository includes the following Python scripts:

- `design_chart_imperial.py`: Generates design charts in imperial units (feet, pounds).
- `design_chart_metric.py`: Generates design charts in metric units (meters, kilograms).
- `fox_embedment.py`: Calculates embedment factors using Fox's method.
- `influence_factor.py`: Computes influence factors for settlement calculations.
- `weighted_modulus.py`: Calculates weighted modulus based on soil layers.

### Configuration Files

The YAML configuration file allows customization of the parameters used in settlement and bearing capacity calculations. The parameters are grouped as follows:

1. **Settlement Parameters**:

   - `D`: Embedment depth of the foundation.
   - `shape`: Aspect ratio of the footing, given as L / B ratio
   - `settlement_location`: Location for calculating settlement (`center`, `edge`, or `corner`)
   - `resistance_factor`: Resistance Factor (or Inverse of Factor of Safety if ASD is used).

2. **Soil Properties**:

   - `gamma_backfill`: Unit weight of structural backfill (imperial: pcf, metric: kN/m³).
   - `gamma_soil`: Unit weight of foundation soil (imperial: pcf, metric: kN/m³).
   - `Z`: Depth to the hard layer (imperial: ft, metric: m).
   - `phi`: Friction angle of the foundation soil (degrees).
   - `c`: Cohesion of the foundation soil (imperial: psf, metric: kPa).
   - `nu`: Poisson's ratio of the foundation soil.

3. **Modulus Options**:

   - `modulus_file`: Set to `True` if modulus values vary by layer and are provided in a separate file.
   - `file`: Path to a CSV file with layer information, required if `modulus_file` is `True`. This file should have columns:
     - `layer_number`: Numeric identifier for the layer.
     - `top`: Depth to the top of the layer (ft or m).
     - `bottom`: Depth to the bottom of the layer (ft or m).
     - `modulus`: Modulus value for the layer (ksf or kPa).
   - `Es`: Optional constant modulus value if `modulus_file` is `False`.

4. **Plot Options**:

   - `max_footing_width`: Set to the maximum footing width to be plotted on the design chart (imperial: ft, metric: m).
   - `max_bearing_pressure`: Set to the maximum bearing pressure to be plotted on the design chart (imperial: ksf, metric: kPa).
   - `X_ticks`: Interval for major ticks on the x-axis.
   - `Y_ticks`: Interval for major ticks on the y-axis.

### Dependencies

The project requires the following Python libraries:

- `numpy`
- `pandas`
- `matplotlib`
- `pyyaml`

You can install them using:

```bash
pip install -r requirements.txt
```

### Running the Code

1. Set the desired parameters in the YAML configuration file.
2. Run the script with the name of the configuration file (.yaml extension optional):

   **Imperial Units**:
   ```bash
   python design_chart_imperial.py <config_file>
   ```
   **Metric Units**:
   ```bash
   python design_chart_metric.py <config_file>
   ```
   Replace `<config_file>` with the name of your YAML configuration file, e.g., `settle_config_imperial` or `settle_config_metric`.

### Outputs

Each script generates the following output files:

- `figure.png`: A plot of the bearing resistance vs. footing width.
- `output_results.xlsx`: An Excel file containing calculated results, including footing width, bearing resistance, modulus, and influence factors.
- `input_information.xlsx`: An Excel file with input values used in the calculations.

### Additional Notes

- Ensure the CSV file for modulus values follows the specified structure if `modulus_file` is set to `True`.
- The `fox_embedment.py` and `influence_factor.py` modules provide essential calculations for settlement and embedment factors, respectively. Do not modify them unless necessary.
- Separate configuration files are provided for imperial (`settle_config_imperial.yaml`) and metric (`settle_config_metric.yaml`) units. Ensure you use the appropriate script and configuration file based on your system of measurement.
