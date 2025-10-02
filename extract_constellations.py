#!/usr/bin/env python3
"""
Utilitys to extract constellation names and data from JSON files created by  obsSky.py or obsSky_batch.py.

Utilitys include:

- load_json_data(json_file): 
    Load JSON data from file

- extract_constellations(data): 
    Extract unique set of constellations from the loaded JSON data

- extract_sundelta_for_constellation(data, target_constellation): 
    Extract sunDelta values for a specific constellation from the loaded JSON data
- extract_table_for_constellation_sundelta(data, target_constellation, target_sundelta): 
    Extract table data for a specific constellation and sunDelta value

- create_dataframe_for_constellation_sundelta(data, target_constellation, target_sundelta): 
    Create a pandas DataFrame for a specific constellation and sunDelta value

Plotting functions:

- plot_efftot_vs_localtime(data, constellation, sundelta, ax=None, my_color='blue', indexes=[0,1,2]): 
    Plot EffTot_0, EffTot_1, EffTot_2 vs localTime for a given constellation and sunDelta

- plot_efftot_constellations(data, sundelta, efftot_index=0, constellation_pattern="OW_*_6_2", ax=None): 
    Plot one EffTot vs localTime for all constellations matching a pattern

- plot_efftot_multi_index(data, sundelta, constellation_pattern="OW_*_6_2", ax=None, indexes=[2]): 
    Plot multiple EffTot indices vs localTime for all constellations matching a pattern

"""

import json
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def load_json_data(json_file):
    """Load JSON data from file"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"Error: File {json_file} not found!")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON in file {json_file}")
        return None

def extract_constellations(data):
    """Extract unique set of constellations from the loaded JSON data"""
    
    if data is None:
        return set()
    
    # Extract constellations from results section
    constellations = set()
    
    if 'results' in data:
        for result in data['results']:
            if 'CONSTELLATIONS' in result:
                # Add all constellations from this result
                for constellation in result['CONSTELLATIONS']:
                    constellations.add(constellation)
    
    return constellations

def extract_sundelta_for_constellation(data, target_constellation):
    """Extract sunDelta values for a specific constellation from the loaded JSON data"""
    
    if data is None:
        return set()
    
    sundelta_values = set()
    
    if 'results' in data:
        for result in data['results']:
            if 'CONSTELLATIONS' in result and 'sunDelta' in result:
                # Check if this result contains the target constellation
                constellations = result['CONSTELLATIONS']
                if target_constellation in constellations:
                    sundelta_values.add(result['sunDelta'])
    
    return sundelta_values

def extract_table_for_constellation_sundelta(data, target_constellation, target_sundelta):
    """Extract table data for a specific constellation and sunDelta value
    
    Returns a list of dictionaries with columns:
    - localTime = sunAlpha/15 + 12
    - sunElev
    - aircut_0, aircut_1, aircut_2 (3 columns)
    - EffTot_0, EffTot_1, EffTot_2 (3 columns)
    - dens_zenith
    - densobs_zenith
    """
    
    if data is None:
        return []
    
    table_data = []
    
    if 'results' in data:
        for result in data['results']:
            # Check if this result matches our criteria
            if ('CONSTELLATIONS' in result and 
                'sunDelta' in result and 
                'sunAlpha' in result and
                'sunElev' in result):
                
                constellations = result['CONSTELLATIONS']
                sun_delta = result['sunDelta']
                
                if (target_constellation in constellations and 
                    sun_delta == target_sundelta):
                    
                    # Calculate localTime
                    sun_alpha = result['sunAlpha']
                    local_time = sun_alpha / 15.0 + 12.0
                    
                    # Create row dictionary
                    row = {
                        'localTime': local_time,
                        'sunElev': result['sunElev']
                    }
                    
                    # Add aircut columns (if available)
                    if 'aircut' in result:
                        aircut = result['aircut']
                        for i in range(min(3, len(aircut))):
                            row[f'aircut_{i}'] = aircut[i]
                    
                    # Add EffTot columns (if available)  
                    if 'EffTot' in result:
                        efftot = result['EffTot']
                        for i in range(min(3, len(efftot))):
                            row[f'EffTot_{i}'] = efftot[i]
                    
                    # Add dens_zenith and densobs_zenith (if available)
                    if 'dens_zenith' in result:
                        row['dens_zenith'] = result['dens_zenith']
                    
                    if 'densobs_zenith' in result:
                        row['densobs_zenith'] = result['densobs_zenith']
                    
                    table_data.append(row)
    
    # Sort by localTime for better readability
    table_data.sort(key=lambda x: x['localTime'])
    
    return table_data

def create_dataframe_for_constellation_sundelta(data, target_constellation, target_sundelta):
    """Create a pandas DataFrame for a specific constellation and sunDelta value
    
    Returns a DataFrame with columns:
    - localTime = sunAlpha/15 + 12
    - sunElev
    - aircut_0, aircut_1, aircut_2 (3 columns)
    - EffTot_0, EffTot_1, EffTot_2 (3 columns)
    - dens_zenith
    - densobs_zenith
    """
    
    table_data = extract_table_for_constellation_sundelta(data, target_constellation, target_sundelta)
    
    if not table_data:
        return pd.DataFrame()
    
    df = pd.DataFrame(table_data)
    
    # Sort by localTime
    df = df.sort_values('localTime').reset_index(drop=True)
    
    return df

def plot_efftot_vs_localtime(data, constellation, sundelta, ax=None, my_color='blue', indexes=[0,1,2]):
    """Plot EffTot_0, EffTot_1, EffTot_2 vs localTime for a given constellation and sunDelta
    
    Args:
        data: Loaded JSON data
        constellation: Target constellation name
        sundelta: Target sunDelta value
        ax: Matplotlib axis (if None, creates new figure)
    
    Returns:
        matplotlib axis object
    """
    
    # Get the data for this constellation and sunDelta
    df = create_dataframe_for_constellation_sundelta(data, constellation, sundelta)
    
    if df.empty:
        print(f"No data found for constellation {constellation} with sunDelta {sundelta}")
        return None
    
    # Create figure and axis if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # Sort by localTime for proper plotting
    df = df.sort_values('localTime')
    
    # Plot EffTot columns in blue with different alpha values (similar to plot_averageEffects.py)
    alphas = [0.3, 0.5, 0.7]  # Increase alpha for higher indices
    
    for icol in indexes:
        if icol == 99:
            continue
        col = f'EffTot_{icol}'
        if col in df.columns:
            # Create label from aircut value if available
            label = f'{constellation} {df[f"aircut_{icol}"].iloc[0]:.0f}°' if f'aircut_{icol}' in df.columns else f'EffTot_{icol}'
            ax.plot(df['localTime'], df[col], 
                    color=my_color, alpha=alphas[icol], 
                    #marker='o', 
                    linestyle='-', linewidth=1, label=label)
    
    # Also plot zenith densities in red if available
    #if 'dens_zenith' in df.columns:
    #    ax.plot(df['localTime'], df['dens_zenith'], 'r-', alpha=0.7, 
    #           label='Zenith')
    
    
    # Calculate max value for vertical line height
    efftot_cols = [col for col in ['EffTot_0', 'EffTot_1', 'EffTot_2'] if col in df.columns]
    if efftot_cols:
        maxy = np.nanmax(df[efftot_cols].values) * 1.1
    else:
        maxy = 1.0
    
    # Add vertical lines at sun elevation crossings (similar to plot_averageEffects.py)
    if 'sunElev' in df.columns:
        sun_elevations = [0, -6, -12, -18]
        for elev in sun_elevations:
            # Find crossings by looking for sign changes
            crossings = []
            for i in range(len(df) - 1):
                sun_elev_vals = df['sunElev'].iloc[i:i+2].values
                local_time_vals = df['localTime'].iloc[i:i+2].values
                
                if len(sun_elev_vals) == 2 and len(local_time_vals) == 2:
                    if ((sun_elev_vals[0] - elev) * (sun_elev_vals[1] - elev)) <= 0:
                        # Linear interpolation to find exact crossing point
                        x1, y1 = local_time_vals[0], sun_elev_vals[0]
                        x2, y2 = local_time_vals[1], sun_elev_vals[1]
                        if y2 != y1:  # Avoid division by zero
                            x_cross = x1 + (elev - y1) * (x2 - x1) / (y2 - y1)
                            crossings.append(x_cross)
            
            # Plot vertical lines at crossings
            for x_cross in crossings:
                #ax.axvline(x=x_cross, color='blue', alpha=0.3, linestyle='--', linewidth=1)
                ax.fill_betweenx([-maxy**0.05, 1.], 0., x_cross, color='blue', alpha=0.1)

                ax.text(x_cross, 1e-4, f'{elev}°', rotation=90, verticalalignment='bottom', 
                       color='blue', alpha=0.7, fontsize=9)
    
    # Formatting

    ax.set_ylim(-maxy*0.05, maxy)
    ax.set_xlabel('Local Time [hours]')
    ax.set_ylabel('Fraction of frames with trails')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    
    
    # Set reasonable x-axis limits and ticks
    ax.set_xlim(16.95, 24.05)
    hour_ticks = np.arange(17, 25)
    ax.set_xticks(hour_ticks)
    
    # Create second x-axis on top for sun elevation if available
    if 'sunElev' in df.columns and len(df) > 1:
        ax2 = ax.twiny()
        
        # Set the same xlim as the main axis
        ax2.set_xlim(ax.get_xlim())
        
        # Create interpolation for sun elevation
        x_min, x_max = ax.get_xlim()
        hour_ticks = np.arange(int(x_min), int(x_max) + 1)
        
        # Interpolate sun elevation for these hour positions
        if len(df) > 1:
            elev_ticks = np.interp(hour_ticks, df['localTime'].values, df['sunElev'].values)
            
            # Set ticks and labels
            ax2.set_xticks(hour_ticks)
            ax2.set_xticklabels([f'{elev:.0f}°' for elev in elev_ticks])
            ax2.set_xlabel('Sun Elevation')
            
            # Make the top axis labels smaller and less prominent
            ax2.tick_params(axis='x', labelsize=9, colors='gray')
    
    return ax

def plot_efftot_constellations(data, sundelta, efftot_index=0, constellation_pattern="OW_*_6_2", ax=None):
    """Plot one EffTot vs localTime for all constellations matching a pattern
    
    Args:
        data: Loaded JSON data
        sundelta: Target sunDelta value  
        efftot_index: Which EffTot to plot (0, 1, or 2)
        constellation_pattern: Pattern to match constellations (e.g., "OW_*_6_2")
        ax: Matplotlib axis (if None, creates new figure)
    
    Returns:
        matplotlib axis object
    """
    
    # Extract all constellations matching the pattern
    all_constellations = extract_constellations(data)
    
    # Filter constellations matching the pattern
    import re
    # Convert pattern to regex (replace * with \d+ to match numbers)
    regex_pattern = constellation_pattern.replace("*", r"(\d+)")
    matching_constellations = []
    altitude_values = []
    
    for constellation in all_constellations:
        match = re.match(regex_pattern, constellation)
        if match:
            matching_constellations.append(constellation)
            # Extract the number and convert to altitude in km
            altitude_km = int(match.group(1)) * 100
            altitude_values.append(altitude_km)
    
    if not matching_constellations:
        print(f"No constellations found matching pattern: {constellation_pattern}")
        return None
    
    # Sort by altitude for consistent ordering (do this early!)
    sorted_data = sorted(zip(matching_constellations, altitude_values), key=lambda x: x[1])
    matching_constellations = [x[0] for x in sorted_data]
    altitude_values = [x[1] for x in sorted_data]
    
    print(f"Found {len(matching_constellations)} constellations matching '{constellation_pattern}' (sorted by altitude):")
    for const, alt in zip(matching_constellations, altitude_values):
        print(f"  - {const} -> {alt} km")
    
    # Create figure and axis if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create colormap 
    cmap = plt.cm.summer  # 'plasma', 'jet', viridis, etc.
    colors = cmap(np.linspace(0, 1, len(matching_constellations)))
    
    # Create final sorted data with colors
    sorted_data = list(zip(matching_constellations, altitude_values, colors))
    
    maxy = 0  # Track maximum value for plot scaling
    
    for constellation, altitude_km, color in sorted_data:
        # Get the data for this constellation and sunDelta
        df = create_dataframe_for_constellation_sundelta(data, constellation, sundelta)
        
        if df.empty:
            print(f"No data found for constellation {constellation} with sunDelta {sundelta}")
            continue
        
        # Sort by localTime for proper plotting
        df = df.sort_values('localTime')
        
        # Plot the specified EffTot column
        efftot_col = f'EffTot_{efftot_index}'
        if efftot_col in df.columns:
            label = f'{altitude_km} km'
            ax.plot(df['localTime'], df[efftot_col], 
                    color=color, 
                    label=label)
            
            # Update maximum for scaling
            col_max = np.nanmax(df[efftot_col].values)
            if col_max > maxy:
                maxy = col_max
    
    # Set up the plot
    maxy *= 1.1  # Add some margin
    
    # Add vertical lines for sun elevation crossings if we have sun elevation data
    # Use data from the first constellation that has sun elevation data
    for constellation, _, _ in sorted_data:
        df = create_dataframe_for_constellation_sundelta(data, constellation, sundelta)
        if not df.empty and 'sunElev' in df.columns:
            sun_elevations = [0, -6, -12, -18]
            for elev in sun_elevations:
                crossings = []
                for i in range(len(df) - 1):
                    sun_elev_vals = df['sunElev'].iloc[i:i+2].values
                    local_time_vals = df['localTime'].iloc[i:i+2].values
                    
                    if len(sun_elev_vals) == 2 and len(local_time_vals) == 2:
                        if ((sun_elev_vals[0] - elev) * (sun_elev_vals[1] - elev)) <= 0:
                            x1, y1 = local_time_vals[0], sun_elev_vals[0]
                            x2, y2 = local_time_vals[1], sun_elev_vals[1]
                            if y2 != y1:
                                x_cross = x1 + (elev - y1) * (x2 - x1) / (y2 - y1)
                                crossings.append(x_cross)
                
                # Plot vertical lines at crossings
                for x_cross in crossings:
                    ax.fill_betweenx([0, 1], 0., x_cross, color='blue', alpha=0.1)
                    ax.text(x_cross, 2e-5, f'{elev}°', rotation=90,     
                            verticalalignment='bottom', 
                            color='blue', alpha=0.7, fontsize=9)
            break  # Only need to do this once
    
    # Formatting
    ax.set_xlabel('Local Time [hours]')
    ax.set_ylabel('Fraction of frames with trails')
    ax.set_xlim(17.,24.1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', title='Altitude')
    
    # Get aircut value for the title
    aircut_label = ""
    for constellation, _, _ in sorted_data:
        df = create_dataframe_for_constellation_sundelta(data, constellation, sundelta)
        if not df.empty and f'aircut_{efftot_index}' in df.columns:
            aircut_val = df[f'aircut_{efftot_index}'].iloc[0]
            aircut_label = f" (aircut {aircut_val:.0f}°)"
            break
    
    ax.set_title(f'EffTot_{efftot_index}{aircut_label} vs Local Time\nSun Declination: {sundelta}°, Pattern: {constellation_pattern}')
    
    # Set reasonable axis limits
    if maxy > 0:
        ax.set_ylim(0, maxy)
    
    return ax

def plot_efftot_multi_index(data, sundelta, constellation_pattern="OW_*_6_2", ax=None, indexes=[2]):
    """Plot multiple EffTot indices vs localTime for all constellations matching a pattern
    
    Each EffTot index uses a different colormap:
    - EffTot_0: Reds colormap
    - EffTot_1: GnBu colormap  
    - EffTot_2: YlGn colormap
    
    Args:
        data: Loaded JSON data
        sundelta: Target sunDelta value  
        constellation_pattern: Pattern to match constellations (e.g., "OW_*_6_2")
        ax: Matplotlib axis (if None, creates new figure)
    
    Returns:
        matplotlib axis object
    """
    
    # Extract all constellations matching the pattern
    all_constellations = extract_constellations(data)
    
    # Filter constellations matching the pattern
    import re
    # Convert pattern to regex (replace * with \d+ to match numbers)
    regex_pattern = constellation_pattern.replace("*", r"(\d+)")
    matching_constellations = []
    altitude_values = []
    


    for constellation in all_constellations:
        match = re.match(regex_pattern, constellation)
        if match:
            matching_constellations.append(constellation)
            # Extract the number and convert to altitude in km
            altitude_km = int(match.group(1)) * 100
            altitude_values.append(altitude_km)
    
    if not matching_constellations:
        print(f"No constellations found matching pattern: {constellation_pattern}")
        return None
    
    # Sort by altitude for consistent ordering (do this early!)
    sorted_data = sorted(zip(matching_constellations, altitude_values), key=lambda x: x[1])
    matching_constellations = [x[0] for x in sorted_data]
    altitude_values = [x[1] for x in sorted_data]
    
    print(f"Multi-EffTot plot: Found {len(matching_constellations)} constellations matching '{constellation_pattern}' (sorted by altitude):")
    for const, alt in zip(matching_constellations, altitude_values):
        print(f"  - {const} -> {alt} km")
    
    # Create figure and axis if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define colormaps for each EffTot index
    colormaps = {
        0: plt.cm.gnuplot,
        1: plt.cm.gnuplot, 
        2: plt.cm.gnuplot,
        99: plt.cm.copper
    }
    

    linestyles = {
        0: "-",
        1: "-", 
        2: "-",
        99: ":"
    } # - solid, dashed --, dotted :

    linewidths = {
        0: 0.5,
        1: 0.5, 
        2: 1.,  
        99: 1.
    }

    alphas = {
        0: 0.5,
        1: 0.5, 
        2: 0.75,    
        99: 0.75
    }



    maxy = 0  # Track maximum value for plot scaling
    aircut_labels = {}  # Store aircut values for each index
    
    # Iterate through EffTot indices
    for efftot_index in indexes:
        cmap = colormaps[efftot_index]
        colors = cmap(np.linspace(0., 0.7, len(matching_constellations)))  # Start from 0.3 to avoid too light colors
        
        # Get aircut value for this index (from first constellation that has it)
        for constellation, altitude_km in zip(matching_constellations, altitude_values):
            df = create_dataframe_for_constellation_sundelta(data, constellation, sundelta)
            if not df.empty and f'aircut_{efftot_index}' in df.columns:
                aircut_val = df[f'aircut_{efftot_index}'].iloc[0]
                aircut_labels[efftot_index] = aircut_val
                break
        
        # Plot curves for this EffTot index
        for i, (constellation, altitude_km) in enumerate(zip(matching_constellations, altitude_values)):
            # Get the data for this constellation and sunDelta
            df = create_dataframe_for_constellation_sundelta(data, constellation, sundelta)
            
            if df.empty:
                continue
            
            # Sort by localTime for proper plotting
            df = df.sort_values('localTime')
            
            # Plot the specified EffTot column
            efftot_col = f'EffTot_{efftot_index}' if efftot_index < 99 else "dens_zenith"
            if efftot_col in df.columns:
                color = colors[i]
                
                # Create label: include altitude for first curve of each EffTot, EffTot info for first altitude
                if i == 0:  # First altitude - include EffTot info
                    aircut_info = f"{aircut_labels.get(efftot_index, '?')}°" if efftot_index in aircut_labels else ""
                    label = f'{altitude_km} km {aircut_info}'
                    #label = aircut_info
                    if efftot_index == 99:
                        label = f'{altitude_km} km Zenith'
                elif i == len(matching_constellations) - 1: # Last altitude - include altitude only
                    label =  f'{altitude_km} km'
                else:  # Other combinations - no label to avoid cluttering
                    label = None
                
                
                ax.plot(df['localTime'], df[efftot_col], 
                        color=color, 
                        linewidth=linewidths[efftot_index],
                        linestyle=linestyles[efftot_index],
                        label=label, 
                        alpha=alphas[efftot_index] )
                
                # Update maximum for scaling
                col_max = np.nanmax(df[efftot_col].values)
                if col_max > maxy:
                    maxy = col_max
    

    # Set up the plot
    maxy *= 1.1  # Add some margin
    
    # Add vertical lines for sun elevation crossings
    # Use data from the first constellation that has sun elevation data
    for constellation, _ in zip(matching_constellations, altitude_values):
        df = create_dataframe_for_constellation_sundelta(data, constellation, sundelta)
        if not df.empty and 'sunElev' in df.columns:
            sun_elevations = [0, -6, -12, -18]
            for elev in sun_elevations:
                crossings = []
                for i in range(len(df) - 1):
                    sun_elev_vals = df['sunElev'].iloc[i:i+2].values
                    local_time_vals = df['localTime'].iloc[i:i+2].values
                    
                    if len(sun_elev_vals) == 2 and len(local_time_vals) == 2:
                        if ((sun_elev_vals[0] - elev) * (sun_elev_vals[1] - elev)) <= 0:
                            x1, y1 = local_time_vals[0], sun_elev_vals[0]
                            x2, y2 = local_time_vals[1], sun_elev_vals[1]
                            if y2 != y1:
                                x_cross = x1 + (elev - y1) * (x2 - x1) / (y2 - y1)
                                crossings.append(x_cross)
                
                # Plot vertical lines at crossings
                for x_cross in crossings:
                    ax.fill_betweenx([0, 1], 0., x_cross, color='blue', alpha=0.1)
                    ax.text(x_cross, 1e-4, f'{elev}°', rotation=90, verticalalignment='top', 
                           color='blue', alpha=0.7, fontsize=9)
            break  # Only need to do this once
    
    # Formatting
    ax.set_xlabel('Local Time [hours]')
    ax.set_ylabel('Fraction of frames with trails')
    ax.set_xlim(17.,24.1)
    ax.grid(True, alpha=0.3)
    
    ax.legend(loc='upper right', 
              #title='Altitude & elevation cut-off',               
              fontsize=9, title_fontsize=10)
    


    #ax.set_title(f'Multiple EffTot Indices vs Local Time\nSun Declination: {sundelta}°, Pattern: {constellation_pattern}\nReds: EffTot_0, Blue-Green: EffTot_1, Yellow-Green: EffTot_2')

    #inclination
    inclination = "60°" if "_6_" in constellation_pattern else "87.9°"
    antenna = f"{constellation_pattern[-1]}0°"
    ax.text(24., 1e-3, 
            f'Sun Declination: {sundelta}°',
            horizontalalignment='right', 
            verticalalignment='top', 
            color='k', alpha=1.)
    
    # Set reasonable axis limits
    if maxy > 0:
        ax.set_ylim(0, maxy)
    
    return ax


#============================================================================
# Main demonstration function
#============================================================================
def demo():
    json_file = "/home/ohainaut/Dropbox/ohainaut/Documents/ESO/E2E/SatelliteConstellations/SIMULATIONS/Eutelsat/obsSky_batch_results_fullOWeq.json"
    
    print("DEMO")

    # Load JSON data once
    print("Loading JSON data...")
    data = load_json_data(json_file)
    
    if data is None:
        sys.exit(1)
    
    # Extract unique constellations from the loaded data
    constellations = extract_constellations(data)
    
    print(f"Found {len(constellations)} unique constellations in the results section")
    #for constellation in sorted(constellations):
    #     print(f"  - {constellation}")

    if 0:
        # Extract sunDelta values for each constellation
        print("\n" + "="*60)
        print("SUN DELTA VALUES BY CONSTELLATION")
        print("="*60)
        
        # Iterate through each constellation and get its sunDelta values
        for constellation in sorted(constellations):
            sundelta_values = extract_sundelta_for_constellation(data, constellation)
            sundelta_list = sorted(sundelta_values)
            
            print(f"{constellation}: sunDelta values: {sundelta_list}")
        
        # Demonstrate table extraction for one constellation and sunDelta
        print("\n" + "="*60)
        print("TABLE EXTRACTION EXAMPLE")
        print("="*60)
    
    # For one constellation and and 3 sunDelta values:
    if 1:
        #my_constellation = 'OW1gen'
        my_constellation = 'OW_12_9_2'
        sundelta_values = extract_sundelta_for_constellation(data, my_constellation)
        print(f"\nConstellation '{my_constellation}' has {len(sundelta_values)} unique sunDelta values: {sorted(sundelta_values)}"  )

        plt.figure(figsize=(12, 8))

        for i, my_sundelta in enumerate(sorted(sundelta_values, reverse=True)):
            ax = plt.subplot(3, 1, i+1)

            # Set title
            if i == 0:
                ax.set_title(f'Constellation: {my_constellation}')
    
            print(f"\nExtracting table for:")
            print(f"  Constellation: {my_constellation}")
            print(f"  sunDelta: {my_sundelta}")
            
            table_data = extract_table_for_constellation_sundelta(data, my_constellation, my_sundelta)
            
            print(f"\nExtracted {len(table_data)} rows:")
            print("="*80)
            
            # Print header
            if table_data:

                if 0:
                    print(f"{'localTime':>10} {'sunElev':>8} {'aircut_0':>8} {'aircut_1':>8} {'aircut_2':>8} {'EffTot_0':>10} {'EffTot_1':>10} {'EffTot_2':>10} {'dens_zen':>10} {'densobs_zen':>11}")
                    print("-" * 100)
                    
                    # Print first few rows as example
                    for i, row in enumerate(table_data[:5]):  # Show first 5 rows
                        print(f"{row.get('localTime', 'N/A'):>10.1f} "
                            f"{row.get('sunElev', 'N/A'):>8.1f} "
                            f"{row.get('aircut_0', 'N/A'):>8.1f} "
                            f"{row.get('aircut_1', 'N/A'):>8.1f} "
                            f"{row.get('aircut_2', 'N/A'):>8.1f} "
                            f"{row.get('EffTot_0', 'N/A'):>10.6f} "
                            f"{row.get('EffTot_1', 'N/A'):>10.6f} "
                            f"{row.get('EffTot_2', 'N/A'):>10.6f} "
                            f"{row.get('dens_zenith', 'N/A'):>10.6f} "
                            f"{row.get('densobs_zenith', 'N/A'):>11.6f}")
                    
                    if len(table_data) > 5:
                        print(f"... and {len(table_data) - 5} more rows")
                
                
                print(f"Creating plot for {my_constellation} with sunDelta = {my_sundelta}")
                
                # Create the plot
                plot_efftot_vs_localtime(data, my_constellation, my_sundelta, ax)

                ax.set_ylim(5e-5, 5e-2)
                ax.set_yscale('log')
                ax.text(17.,1e-4, f'Sun Declination: {my_sundelta}°')

                
            else:
                print("No data available to create plots.")
        
        # Save the plot
        comparison_filename = f"{my_constellation}_sun.png"
        plt.tight_layout()
        plt.savefig(comparison_filename, dpi=300, bbox_inches='tight')
        print(f"Constellation comparison plot saved as: {comparison_filename}")

        # Optionally show the plot
        plt.show()
                

    if 0:
        my_constellation = sorted(constellations)[1]
        sundelta_values = extract_sundelta_for_constellation(data, my_constellation)
        if sundelta_values:
            my_sundelta = sorted(sundelta_values)[1]
            
            print(f"\nExtracting table for:")
            print(f"  Constellation: {my_constellation}")
            print(f"  sunDelta: {my_sundelta}")
            
            table_data = extract_table_for_constellation_sundelta(data, my_constellation, my_sundelta)
            
            print(f"\nExtracted {len(table_data)} rows:")
            print("="*80)
            if table_data:
                # Demonstrate multi-EffTot plotting
                print("\n" + "="*60)
                print("MULTI-EFFTOT COMPARISON PLOTTING")
                print("="*60)
                
                print(f"Creating multi-EffTot comparison plot for sunDelta = {my_sundelta}")
                print("Pattern: OW_*_6_2 (comparing all EffTot indices with different colormaps)")
                print("  - EffTot_0: Red colormap (solid lines with markers)")
                print("  - EffTot_1: Blue-Green colormap (dashed lines)")
                print("  - EffTot_2: Yellow-Green colormap (dotted lines)")
                
                # Create the multi-EffTot comparison plot
                fig, ax = plt.subplots(figsize=(14, 10))
                plot_efftot_multi_index(data, my_sundelta, 
                                       constellation_pattern="OW_*_6_4", ax=ax)
                
                # Save the plot
                multi_filename = f"OW_multi_EffTot_comparison_sunDelta_{my_sundelta:+.1f}.png"
                plt.tight_layout()
                plt.savefig(multi_filename, dpi=300, bbox_inches='tight')
                print(f"Multi-EffTot comparison plot saved as: {multi_filename}")
                
                # Optionally show the plot
                # plt.show()
                



#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
if __name__ == "__main__":
    #demo()
#if 0:
    json_file = "/home/ohainaut/Dropbox/ohainaut/Documents/ESO/E2E/SatelliteConstellations/SIMULATIONS/Eutelsat/obsSky_batch_results_fullOWeq.json"

    data = load_json_data(json_file)


    icount = 0
    constellation_patterns = []
    for inc in ["6","9"]:
        for mea in ["2","3","4"]:
            constellation_pattern=f"OW_*_{inc}_{mea}"
            sunDeltas = [-23.5, 0.0, 23.5]
            icount += 1
            constellation_patterns.append(constellation_pattern)    
    print(constellation_patterns)


    indexes = [2]
    for constellation_pattern in constellation_patterns:
            print(f"\n[{icount:02d}] Plotting constellation pattern: {constellation_pattern}")
            
            if "*_6_" in constellation_pattern:
                const_label = "OW Eq. Inc=60° "
            else:
                const_label = "OW Eq. Inc=87.9° "
            const_label += f'MEA={constellation_pattern[-1]}0°'
            const_file = f'OW{constellation_pattern[-4:]}'

            plt.figure(figsize=(12, 8))
            for i, sunDelta in enumerate(sunDeltas):
                ax = plt.subplot(3, 1, i+1)
                if i == 0:
                    ax.set_title(const_label, fontsize=14)


                #over plot OW1gen for reference
                extract_table_for_constellation_sundelta(data, 'OW1gen', sunDelta)
                ax = plot_efftot_vs_localtime(data, 'OW1gen', sunDelta, ax=ax, my_color='grey', indexes=indexes)


                plot_efftot_multi_index(data, sunDelta, 
                                        constellation_pattern=constellation_pattern, ax=ax, indexes=indexes)

                if i < 2:
                    ax.set_xlabel("")  # Remove x-label for top plots to avoid clutter

                if i == 0:
                    ax.set_title(const_label, fontsize=14)

                ax.set_ylim(5e-5, 1e-2)
                ax.set_yscale('log')
            plt.tight_layout()

            output_file = f'{const_file}_comparison.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved as {output_file}")         
            #plt.show(block=False)
