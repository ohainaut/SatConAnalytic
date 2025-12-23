#!/usr/bin/env python3
"""
Example script showing how to use the refactored obsSky.py main function
"""

# Import the refactored module (assuming it's saved as obsSky_refactored.py)
import obsSky
import conan as caLib
import constellations
import json
import numpy as np
from datetime import datetime

def OW():
    CONSTELLATIONS = constellations.readConstellations()
    results = []

    for c in CONSTELLATIONS.list:  
        print("Constellation:",c)
        for deltaSun in [-23.5, 0.,  23.5]:
            print("  DeltaSun:",deltaSun)
            for alphaSun in range(75,181,3):
                elevSun = alphaSun
                print("    ElevSun:",elevSun)
                args = []
                args.extend(['-T', 'LSST']) 
                args.extend(['-C', c])
                args.extend(['-d', str(deltaSun)])
                args.extend(['-a', str(alphaSun)])
                args.extend(['--noplot'])  # Don't generate plots for batch runs

                results.append(obsSky.main(args))

    # Save results to JSON file
    def numpy_serializer(obj):
        """JSON serializer for numpy objects"""
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    # Prepare metadata
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'total_simulations': len(results),
        'constellations': CONSTELLATIONS.list[:3],
        'parameter_ranges': {
            'deltaSun': [-23.5, 0., 23.5],
            'alphaSun_range': [90, 181, 60],
            'telescope': 'LSST'
        }
    }
    
    # Create output data structure
    output_data = {
        'metadata': metadata,
        'results': results
    }
    
    # Save to JSON file
    output_filename = f"obsSky_batch_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    try:
        with open(output_filename, 'w') as f:
            json.dump(output_data, f, indent=2, default=numpy_serializer)
        print(f"Results successfully saved to {output_filename}")
    except Exception as e:
        print(f"Error saving results to JSON: {e}")
        # Fallback: save a simplified version
        simplified_results = []
        for i, result in enumerate(results):
            if result:
                simplified_results.append({
                    'index': i,
                    'outfileroot': result.get('outfileroot', 'N/A'),
                    'sunAlpha': float(result.get('sunAlpha', 0)),
                    'sunDelta': float(result.get('sunDelta', 0)),
                    'sunElev': float(result.get('sunElev', 0))
                })
        
        simplified_data = {
            'metadata': metadata,
            'simplified_results': simplified_results
        }
        
        fallback_filename = f"obsSky_batch_results_simplified_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(fallback_filename, 'w') as f:
            json.dump(simplified_data, f, indent=2)
        print(f"Simplified results saved to {fallback_filename}")
    
    return results

def run_multiple_simulations():
    """Example of running obsSky with different parameters programmatically"""
    
    # Define parameter sets to run
    parameter_sets = [
        {
            'elevSun': 18,
            'deltaSun': 0,
            'constellations': 'SL',
            'code': 'EFFECT',
            'lat': -24.6
        },
        {
            'elevSun': 12,
            'deltaSun': 0,
            'constellations': 'SL',
            'code': 'EFFECT',
            'lat': -24.6
        },
        {
            'elevSun': 6,
            'deltaSun': 0,
            'constellations': 'SL',
            'code': 'EFFECT',
            'lat': -24.6
        }
    ]
    
    results = []
    
    for i, params in enumerate(parameter_sets):
        print(f"\n=== Running simulation {i+1} ===")
        
        # Convert parameters to command line argument format
        args = []
        for key, value in params.items():
            if key == 'elevSun':
                args.extend(['-e', str(value)])
            elif key == 'deltaSun':
                args.extend(['-d', str(value)])
            elif key == 'constellations':
                args.extend(['-C', str(value)])
            elif key == 'code':
                args.extend(['-T', str(value)])
            elif key == 'lat':
                args.extend(['-l', str(value)])
        
        # Add common arguments
        args.extend(['--noplot'])  # Don't generate plots for batch runs
        
        # Run the simulation
        try:
            result = obsSky.main(args)
            results.append(result)
            print(f"Simulation {i+1} completed successfully")
        except Exception as e:
            print(f"Error in simulation {i+1}: {e}")
            results.append(None)
    
    return results

def run_single_simulation():
    """Example of running a single simulation with specific parameters"""
    
    # Define arguments as you would on command line
    args = [
        '-a', '120',
        '-d', '0',         # Sun declination 0 degrees
        '-C', 'OW1gen',        # Starlink constellation
        '-T', 'LSST'
        #,'--noplot'         # Don't generate plot
    ]
    
    print("Running single simulation...")
    result = obsSky.main(args)
    print("Single simulation completed")
    
    return result

def run_series_simulation():
    """Example of running a single simulation scanning a parameter"""
    
    # Define arguments as you would on command line
    args = [
        '-d', '-23',         # Sun declination 0 degrees
        '-C', 'SLOWGWAK',#'TODAY', #'ALL', #'SL1old',        # Starlink constellation
        '-T', 'FORSimg',
        '-M', 'EFFECT',
        #'--nolabel',
        '--noDots',
        '--noalmuc'
    ]
    
    for a in np.arange(90,271.,.25):
        args.extend(['-a', str(a)])  # Sun elevation from 75 to 180 degrees
        result = obsSky.main(args)
    
    return result

if __name__ == "__main__":

    ###results = OW()

    # Example 1: Run a single simulation
    #single_result = run_single_simulation()
    
    # Example 2: Run multiple simulations
    #multiple_results = run_multiple_simulations()
    
    _ = run_series_simulation()

    print("All simulations completed!")
