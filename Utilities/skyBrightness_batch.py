import numpy as np

import skyBrightness
import SatConAnalytic.conanplot as cpLib


def fillArgs(args, CONST, a, o):
    return [
        '-C', CONST,
        '-d', '0',
        '-u', 'nanoLambert',
        '-a', str(a),
        '-o', o
    ]       


toDo = ['satDens','losses', 'skyDiffuse', 'skyScatter']
CONST = 'SLOWGWAK' # 'SunLight2'

def run_plots():
    for o in toDo:
        args = [  
            '-C', CONST,
            '-T', 'FORSimg',
            '-d', '-23',
            '-e', '24',
            '-u', 'nanoLambertPerM2',
            '-o', o]
        if o == 'skyDiffuse':  args.append(['-m', '0'])
        if o == 'satDens':     args.append(['--dots'])
        print(f"Running simulation with args: {args}")



def run_sequence():
    results_ = []  


    for a in np.arange(100.,181.,3.):

        results_a = {}
        for o in :
            args = [  
                '-C', CONST,
                '-d', '0',
                '-u', 'nanoLambert',
                '-a', str(a),
                '-o', o]
            

            print(f"Running simulation with args: {args}")
            results.append(skyBrightness.main(args))

    # write results to file
    import json
    with open(f"skyBrightness_{CONST}.json", "w") as f:
        json.dump(results, f, indent=4, cls=cpLib.NumpyEncoder) 


if __name__ == "__main__":
    results = run_plots( )

    print("All simulations completed!")
