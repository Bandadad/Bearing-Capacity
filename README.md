# Bearing Capacity and Settlement of Shallow Foundations

This project produces a design chart for spread footing using elastic theory, based on Bowles' *Foundation Analysis and Design*, 5th Edition.

### Settlement Calculation

Settlement is estimated from the equation:

\[
\Delta H = \frac{q \cdot B' \cdot (1 - \nu^2)}{Es} \cdot m \cdot I_s \cdot I_f
\]

### Project Structure

The repository includes two main Python scripts for generating design charts:

- `design_chart_imperial.py`: Uses imperial units (feet, pounds).
- `design_chart_metric.py`: Uses metric units (meters, kilograms).

Each script loads configuration parameters from a corresponding YAML file (`settle_config_imperial.yaml` or `settle_config_metric.yaml`).

### Configuration Files

The YAML configuration file allows customization of the parameters used in settlement and bearing capacity calculations. The parameters are grouped as follows:

1. **Settlement Parameters**:
   - `D`: Embedment depth of the foundation.
   - `shape`: Aspect ratio of the footing, given as \( B/L \).
   - `settlement_location`: Location for calculating settlement (`center`, `edge`, or `corner`).

2. **Soil Properties**:
   - `gamma_backfill`: Unit weight of structural backfill (imperial: pcf, metric: kN/m³).
   - `gamma_soil`: Unit weight of foundation soil.
   - `Z`: Depth to the hard layer.
   - `phi`: Friction angle of the foundation soil (degrees).
   - `c`: Cohesion of the foundation soil (imperial: psf, metric: kPa).
   - `nu`: Poisson's ratio of the foundation soil.

3. **Modulus Options**:
   - `modulus_file`: Set to `True` if modulus values vary by layer and are provided in a separate file.
   - `file`: Path to a CSV file with layer information, required if `modulus_file` is `True`. This file should have columns `layer_number`, `top`, `bottom`, and `modulus`.
   - `Es`: Optional constant modulus value if `modulus_file` is `False`.

### Running the Code

1. Set the desired parameters in the YAML configuration file.
2. Run the desired script:
   - For imperial units: `python design_chart_imperial.py`
   - For metric units: `python design_chart_metric.py`

### Outputs

Each script generates the following output files:

- `figure.png`: A plot of the bearing resistance vs. footing width.
- `output_results.xlsx`: An Excel file containing calculated results, including footing width, bearing resistance, modulus, and influence factors.
- `input_information.xlsx`: An Excel file with input values used in the calculations.
