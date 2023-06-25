import pandas as pd

def compute_weighted_average(file, Z):
    # Read the csv file
    df = pd.read_csv(file)
    # Compute the thickness and weight of each layer
    df['thickness'] = df['bottom'] - df['top']
    df['weight'] = df['thickness'] * df['modulus']

    # Filter the layers that are fully within the depth Z
    complete_layers = df[df['bottom'] <= Z].copy()

    # Check if Z falls within a layer
    within_layer = df[(df['top'] < Z) & (df['bottom'] > Z)]

    if not within_layer.empty:
        # Compute the partial thickness and weight of the layer where Z falls within
        partial_thickness = Z - within_layer['top'].values[0]
        partial_weight = partial_thickness * within_layer['modulus'].values[0]

        # Append this partial layer to the complete layers
        partial_layer = within_layer.copy()
        partial_layer['thickness'] = partial_thickness
        partial_layer['weight'] = partial_weight
        complete_layers = complete_layers._append(partial_layer, ignore_index=True)

    if complete_layers.empty:
        df['thickness'] = df['bottom'] - df['top']
        df['weight'] = df['thickness'] * df['modulus']
        complete_layers = df.copy()

    # Calculate the weighted average modulus
    weighted_avg_modulus = complete_layers['weight'].sum() / complete_layers['thickness'].sum()

    return weighted_avg_modulus