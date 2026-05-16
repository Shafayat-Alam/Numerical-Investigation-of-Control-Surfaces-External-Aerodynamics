import os
import pandas as pd
import numpy as np

# Enforce active GUI backend loop for pop-up frame renders
import matplotlib
matplotlib.use('QtAgg' if 'QtAgg' in matplotlib.rcsetup.all_backends else 'TkAgg')
import matplotlib.pyplot as plt

def parse_openfoam_file_with_headers(potential_paths):
    """
    Finds the file path, extracts exact column names from the header row,
    strips non-numeric string tokens (like patch names) row-by-row,
    and returns a clean numeric dataframe with matching headers.
    """
    target_path = None
    for path in potential_paths:
        if os.path.exists(path):
            target_path = path
            break
            
    if not target_path:
        return None
        
    raw_header_cols = []
    numeric_data_rows = []
    
    with open(target_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
                
            # Extract column labels from the OpenFOAM header row definition
            if stripped.startswith('#') and 'time' in stripped.lower():
                raw_header_cols = [c for c in stripped.replace('#', '').split() if c]
                continue
            elif stripped.startswith('#'):
                continue
            
            parts = stripped.split()
            if not parts:
                continue
                
            # Isolate and preserve only numeric fields per line
            cleaned_row = []
            for token in parts:
                clean_token = token.replace('(','').replace(')','').replace(',','')
                try:
                    cleaned_row.append(float(clean_token))
                except ValueError:
                    # Skips text components ('Fin', 'GAMG', 'DILUPBiCGStab', etc.)
                    continue
            
            if cleaned_row:
                numeric_data_rows.append(cleaned_row)
                
    if not numeric_data_rows:
        return None
        
    df = pd.DataFrame(numeric_data_rows)
    
    # Strip metadata names out of the header definition list to align with numeric columns
    numeric_headers = []
    text_blacklist = {'patch', 'solver', 'converged'}
    if raw_header_cols:
        numeric_headers = [h for h in raw_header_cols if not any(x in h.lower() for x in text_blacklist)]
        
    # Match names cleanly if vector dimensions align perfectly
    if len(numeric_headers) == df.shape[1]:
        df.columns = numeric_headers
    else:
        df.columns = ['Time'] + [f'Field_Col_{i}' for i in range(1, df.shape[1])]
        
    return df

# Configure file paths and required exact title variables
file_targets = {
    'Residuals': [
        'postProcessing/residuals/0/solverInfo.dat', 
        'solverInfo.dat'
    ],
    'Aerodynamic Force Coefficients': [
        'postProcessing/forceCoeffs/0/coefficient.dat',
        'postProcessing/coefficients/0/coefficient.dat', 
        'coefficient.dat'
    ],
    'y+': [
        'postProcessing/yPlus/0/yPlus.dat', 
        'yPlus.dat'
    ]
}

print("Initializing cleaned OpenFOAM extraction sequence...\n" + "="*70)

for title, paths in file_targets.items():
    df = parse_openfoam_file_with_headers(paths)
    
    if df is not None and 'Time' in df.columns:
        t = df['Time']
        raw_fields = [col for col in df.columns if col != 'Time']
        
        # Isolate exactly what fields to plot based on target configuration requirements
        if title == 'Residuals':
            # Strict filtering: preserve only fields with 'final' in the string index name
            plot_fields = [f for f in raw_fields if 'final' in f.lower()]
        else:
            plot_fields = raw_fields
            
        print(f"Parsed Successfully -> Title: '{title}'")
        print(f"  Steps: {int(t.min())} to {int(t.max())} | Extracted Metrics: {plot_fields}\n")
        
        # Instantiate separate window layouts
        fig, ax = plt.subplots(figsize=(14, 6.5))
        fig.canvas.manager.set_window_title(title)
        
        for field in plot_fields:
            # Set formatting configurations based on file origin types
            lw = 1.6 if any(x in field.lower() for x in ['final', 'max', 'cd', 'cl', 'cs', 'average']) else 1.0
            ls = '-' if not any(x in field.lower() for x in ['min']) else '--'
            
            ax.plot(t, df[field], label=field, alpha=0.85, lw=lw, linestyle=ls)
            
        ax.set_title(title, fontweight='bold', fontsize=14, pad=12)
        ax.set_xlabel('Iteration / Time Step', fontweight='bold', fontsize=11)
        ax.set_ylabel('Magnitude Scale', fontweight='bold', fontsize=11)
        ax.grid(True, which="both", linestyle=':', alpha=0.5)
        ax.set_xlim(t.min(), t.max())
        
        # Force log scaling strictly on the final residual fields
        if title == 'Residuals':
            ax.set_yscale('log')
            ax.set_ylabel('Log Error Magnitude (Final Values Only)', fontweight='bold', fontsize=11)
            
        # Position legends neatly on the side if tracking vectors are dense
        if len(plot_fields) > 8:
            ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize='small', ncol=2)
        else:
            ax.legend(loc='best', fontsize='medium')
            
        plt.tight_layout()
    else:
        print(f"Skipped -> Target path matching '{title}' was empty or missing on disk.\n")

print("="*70 + "\nAll arrays generated. Launching interactive visualization displays...")
plt.show()