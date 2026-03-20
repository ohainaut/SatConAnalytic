import json
import numpy as np

import obsSky
import SatConAnalytic.conanplot as cpLib



toDo = ['satDens','losses', 'skyDiffuse', 'skyScattered']



def run_elev(CONST):
    d = -23.
    O = 'trails'
    results = []
    for a in np.arange(70,182,1.):# 3.
        args = [  
            '-C', CONST,
            '-T', 'FORSimg',
            '-d', f'{d:.0f}',
            '-a', f'{a}',
            '-O', O,
            '--noPlot']
    
        results.append(obsSky.main(args))

        if results[-1]['value']['zenith'] == 0. and results[-1]['value']['sun30'] == 0.:
            break

    outfile = f"data_{O}_{CONST}_{d:.0f}.json"
    with open(outfile, "w") as f:
        json.dump(results, f, indent=4, cls=cpLib.NumpyEncoder) 

    print(f"Completed data run in {outfile}")

def run_data(CONST):
    results = []
    d = 0.
    O = 'skyScattered'
    for a in np.arange(90,182,3.):
        args = [  
            '-C', CONST,
            '-T', 'FORSimg',
            '-d', f'{d:.0f}',
            '-a', f'{a}',
            '-u', 'nanoLambert',
            '-O', O,
            '--noPlot']
    
        results.append(obsSky.main(args))

        if results[-1]['value']['zenith'] == 0. and results[-1]['value']['sun30'] == 0.:
            break

    outfile = f"data_{O}_{CONST}_{d:.0f}.json"
    with open(outfile, "w") as f:
        json.dump(results, f, indent=4, cls=cpLib.NumpyEncoder) 

    print(f"Completed data run in {outfile}")


def run_4plots(CONST):
    results = []
    d = -23.
    e = 24.
    for O in ['trails', 'skyDiffuse', 'satDens', 'losses']: #, 'skyScattered',]:
    #for O in ['trails']: #, 'skyScattered',]:
        args = [  
            '-C', CONST,
            '-T', 'FORSimg',
            '-d', f'{d}',
            '-e', f'{e}',
            '-u', 'frac',
            '-O', O]
    
        if CONST == 'SXODC': args += ['--constFile', 'McD.json']
        if O == 'skyDiffuse':  args += ['-m', '0.']
        if O in ['satDens', 'trails']:   
            args += ['--dots']

        print(args)
        results.append(obsSky.main(args))



    print(f"Completed data run for 4 plots with CONST={CONST}")
    
if __name__ == "__main__":


    Cs = [
          'SLOWGWAK', 
          'AST', 
          'AST-243', 'AST-3000',
          'SunLight0', 
          'SunLight2',
          'SXODC'
          ]
    


    for CONST in Cs:
        #run_elev(CONST)
        #run_data(CONST) 
        run_4plots(CONST)
    print("All simulations completed!")
