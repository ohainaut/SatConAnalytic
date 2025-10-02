#!/usr/bin/env python3
"""
Plot average effects data from averageEffects.dat CSV file

Plots:
- X-axis: local_time converted to hours as h = (local_time+180.)/15.
- Y-axis: 
  - EffTot_* columns in blue
  - TrailTot_* columns in black  
  - dens_zenith and densobs_zenith in red
- Vertical bars where sunElev crosses 0, -6, -12, and -18
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse

def plot_average_effects(ax, csv_file='averageEffects.dat'):
    """Plot the average effects data"""
    
    # Read the CSV file
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: File {csv_file} not found!")
        return
    
    # Convert local_time to hours
    df['hours'] = df['local_time']
    df['hours'][df['local_time'] < 1.] =  df['hours'][df['local_time'] < 1.] + 24.
    
    # Plot EffTot columns in blue
    for icol in [0, 1, 2]:
        col = f'EffTot_{icol}'
        alpha = 0.3 + 0.2 * icol  # Increase alpha for higher indices
        label = df[f"aircut_{icol}"].iloc[0] if f'aircut_{icol}' in df.columns else ''
        label = f'{label:.0f}°'
        if col in df.columns:
            ax.plot(df['hours'], df[col], 'b-', alpha=alpha, label=label)

    maxy = np.nanmax(df[[col for col in ['EffTot_0', 'EffTot_1', 'EffTot_2'] if col in df.columns]].values) * 1.1

    # # Plot TrailTot columns in black
    # for col in ['TrailTot_0', 'TrailTot_1', 'TrailTot_2']:
    #     if col in df.columns:
    #         ax.plot(df['hours'], df[col], 'k-', alpha=0.7, label=col if col == 'TrailTot_0' else '')
    
    # Plot zenith densities in red
    if 'dens_zenith' in df.columns:
        ax.plot(df['hours'], df['dens_zenith'], 'r-', alpha=0.7, label='zenith')
    # if 'densobs_zenith' in df.columns:
    #     ax.plot(df['hours'], df['densobs_zenith'], 'r--', alpha=0.7, label='densobs_zenith')
    
    # Find where sunElev crosses specific values and add vertical lines
    sun_elevations = [0, -6, -12, -18]
    for elev in sun_elevations:
        # Find crossings by looking for sign changes
        crossings = []
        for i in range(len(df) - 1):
            if ((df['sunElev'].iloc[i] - elev) * (df['sunElev'].iloc[i+1] - elev)) <= 0:
                # Linear interpolation to find exact crossing point
                x1, y1 = df['hours'].iloc[i], df['sunElev'].iloc[i]
                x2, y2 = df['hours'].iloc[i+1], df['sunElev'].iloc[i+1]
                if y2 != y1:  # Avoid division by zero
                    x_cross = x1 + (elev - y1) * (x2 - x1) / (y2 - y1)
                    crossings.append(x_cross)
        
        # Plot vertical lines at crossings
        for x_cross in crossings:
            ax.fill_betweenx([0., maxy], 0., x_cross, color='blue', alpha=0.1)
            ax.text(x_cross, maxy, f'{elev}°', rotation=90, verticalalignment='top', color='blue', alpha=0.7)
    
    # Formatting
    ax.set_xlabel('Hours (Local Time)')
    ax.set_ylabel('Fraction of frames with trails')
    ax.text(0.75, 0.9, f'Sun Declination: {df["sunDelta"].iloc[0]:.1f}°', transform=ax.transAxes, ha='right')
    ax.grid(True, alpha=0.3)
    ax.legend( loc='upper right')

    # Set x-axis to show reasonable hour labels
    ax.set_xticks(np.arange(0, 25))
    ax.set_xlim(17.,24.25)
    
    # Create second x-axis on top for sun elevation
    ax2 = ax.twiny()
    
    # Create interpolation function to map hours to sun elevation
    # Sort by hours for proper interpolation
    df_sorted = df.sort_values('hours')
    
    # Remove duplicates and NaN values
    mask = ~df_sorted['hours'].duplicated() & ~df_sorted['sunElev'].isna()
    hours_clean = df_sorted['hours'][mask].values
    sunelev_clean = df_sorted['sunElev'][mask].values
    
    # Set the same xlim as the main axis
    ax2.set_xlim(ax.get_xlim())
    
    # Create tick positions at regular hour intervals within the plot range
    hour_ticks = np.arange(17, 25, 1)  # Every hour from 17 to 24
    
    # Interpolate sun elevation for these hour positions
    elev_ticks = np.interp(hour_ticks, hours_clean, sunelev_clean)
    
    # Set ticks and labels
    ax2.set_xticks(hour_ticks)
    ax2.set_xticklabels([f'{elev:.0f}°' for elev in elev_ticks])
    ax2.set_xlabel('Sun Elevation')
    
    # Make the top axis labels smaller and less prominent
    ax2.tick_params(axis='x', labelsize=9, colors='gray')

    

def plot_one_file():
    parser = argparse.ArgumentParser(description='Plot average effects data')
    parser.add_argument('--file', '-f', default='averageEffects.dat',
                       help='CSV file to plot (default: averageEffects.dat)')
    
    args = parser.parse_args()

    # Create the plot
    plt.figure(figsize=(12, 8))
    ax = plt.subplot(1,1,1)
    plot_average_effects(ax, args.file)
    plt.tight_layout()
    
    # Save the plot
    output_file = args.file.replace('.dat', '_plot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")

    plt.show()


def plot_all_files():
    constellation = 'OW1gen'

    files = ['OW1gen_averageEffects_+23.dat', 'OW1gen_averageEffects_00.dat', 'OW1gen_averageEffects_-23.dat']
    plt.figure(figsize=(12, 8))


    for i, file in enumerate(files):
        ax = plt.subplot(len(files), 1, i+1)
        if i == 0:
            ax.set_title(f'Constellation: {constellation}', fontsize=14)
        plot_average_effects(ax, file)
        
    plt.tight_layout()

    # Save the plot
    output_file = f'{constellation}_averageEffects_plot.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")

    plt.show()

if __name__ == "__main__":
    #plot_one_file()
    plot_all_files()