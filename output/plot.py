#!/usr/bin/env python3
"""
Plot CMB power spectra from CAMB output files.

Usage: python plot.py [root_name]
Default root_name is 'test'

Data file formats:
- scalCls.dat: l, TT, EE, TE
- vecCls.dat:  l, TT, EE, BB, TE  
- tensCls.dat: l, TT, EE, BB, TE
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_cls_data(filename):
    """Load Cls data from file."""
    try:
        data = np.loadtxt(filename)
        return data
    except FileNotFoundError:
        print(f"Warning: File {filename} not found")
        return None
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def plot_cls(root_name='test', data_dir='data', output_dir='plots'):
    """Plot CMB power spectra."""
    
    # String tension parameter
    Gmu = 2e-7
    normalization = Gmu**2
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # File paths
    scal_file = os.path.join(data_dir, f'{root_name}_scalCls.dat')
    vec_file = os.path.join(data_dir, f'{root_name}_vecCls.dat')
    tens_file = os.path.join(data_dir, f'{root_name}_tensCls.dat')
    
    # Load data
    scal_data = load_cls_data(scal_file)
    vec_data = load_cls_data(vec_file)
    tens_data = load_cls_data(tens_file)
    
    # Set up the plot style
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 12,
        'axes.linewidth': 1.2,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.minor.size': 3,
        'ytick.minor.size': 3,
        'legend.frameon': True,
        'legend.fancybox': True,
        'legend.shadow': True
    })
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'CMB Power Spectra - {root_name}', fontsize=16, fontweight='bold')
    
    # Colors for different spectra types
    colors = {'scalar': 'blue', 'vector': 'red', 'tensor': 'green'}
    
    # Plot TT spectrum
    ax = axes[0, 0]
    if scal_data is not None:
        l_scal = scal_data[:, 0]
        tt_scal = scal_data[:, 1]
        ax.loglog(l_scal, tt_scal * normalization, 
                 color=colors['scalar'], linewidth=2, label='Scalar')
    
    if vec_data is not None:
        l_vec = vec_data[:, 0]
        tt_vec = vec_data[:, 1]
        ax.loglog(l_vec, tt_vec * normalization, 
                 color=colors['vector'], linewidth=2, label='Vector')
    
    if tens_data is not None:
        l_tens = tens_data[:, 0]
        tt_tens = tens_data[:, 1]
        ax.loglog(l_tens, tt_tens * normalization, 
                 color=colors['tensor'], linewidth=2, label='Tensor')
    
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{TT} \cdot G\mu^2/(2\pi)$ [$\mu$K$^2$]')
    ax.set_title('Temperature (TT)', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Plot EE spectrum
    ax = axes[0, 1]
    if scal_data is not None and scal_data.shape[1] > 2:
        ee_scal = scal_data[:, 2]
        ax.loglog(l_scal, ee_scal * normalization, 
                 color=colors['scalar'], linewidth=2, label='Scalar')
    
    if vec_data is not None and vec_data.shape[1] > 2:
        ee_vec = vec_data[:, 2]
        ax.loglog(l_vec, ee_vec * normalization, 
                 color=colors['vector'], linewidth=2, label='Vector')
    
    if tens_data is not None and tens_data.shape[1] > 2:
        ee_tens = tens_data[:, 2]
        ax.loglog(l_tens, ee_tens * normalization, 
                 color=colors['tensor'], linewidth=2, label='Tensor')
    
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{EE} \cdot G\mu^2/(2\pi)$ [$\mu$K$^2$]')
    ax.set_title('E-mode Polarization (EE)', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Plot BB spectrum
    ax = axes[1, 0]
    if vec_data is not None and vec_data.shape[1] > 3:
        bb_vec = vec_data[:, 3]
        # Only plot positive values for log scale
        mask_vec = bb_vec > 0
        if np.any(mask_vec):
            ax.loglog(l_vec[mask_vec], bb_vec[mask_vec] * normalization, 
                     color=colors['vector'], linewidth=2, label='Vector')
    
    if tens_data is not None and tens_data.shape[1] > 3:
        bb_tens = tens_data[:, 3]
        mask_tens = bb_tens > 0
        if np.any(mask_tens):
            ax.loglog(l_tens[mask_tens], bb_tens[mask_tens] * normalization, 
                     color=colors['tensor'], linewidth=2, label='Tensor')
    
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{BB} \cdot G\mu^2/(2\pi)$ [$\mu$K$^2$]')
    ax.set_title('B-mode Polarization (BB)', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Plot TE spectrum
    ax = axes[1, 1]
    if scal_data is not None and scal_data.shape[1] > 3:
        te_scal = scal_data[:, 3]
        # Handle both positive and negative values
        ax.semilogx(l_scal, te_scal * normalization, 
                   color=colors['scalar'], linewidth=2, label='Scalar')
    
    if vec_data is not None and vec_data.shape[1] > 4:
        te_vec = vec_data[:, 4]
        ax.semilogx(l_vec, te_vec * normalization, 
                   color=colors['vector'], linewidth=2, label='Vector')
    
    if tens_data is not None and tens_data.shape[1] > 4:
        te_tens = tens_data[:, 4]
        ax.semilogx(l_tens, te_tens * normalization, 
                   color=colors['tensor'], linewidth=2, label='Tensor')
    
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{TE} \cdot G\mu^2/(2\pi)$ [$\mu$K$^2$]')
    ax.set_title('Temperature-E Cross-correlation (TE)', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Save the plot
    output_file = os.path.join(output_dir, f'{root_name}_cls_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_file}")
    
    # Also save as PDF
    output_file_pdf = os.path.join(output_dir, f'{root_name}_cls_comparison.pdf')
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Plot saved as: {output_file_pdf}")
    
    plt.show()
    
    # Create individual plots for each spectrum type
    create_individual_plots(scal_data, vec_data, tens_data, root_name, output_dir, colors, normalization)

def create_individual_plots(scal_data, vec_data, tens_data, root_name, output_dir, colors, normalization):
    """Create individual plots for TT, EE, BB, and TE spectra."""
    
    spectra = ['TT', 'EE', 'BB', 'TE']
    col_indices = {'scalar': [1, 2, None, 3], 'vector': [1, 2, 3, 4], 'tensor': [1, 2, 3, 4]}
    
    for i, spectrum in enumerate(spectra):
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot each component if available
        for data_type, data, color in [('scalar', scal_data, colors['scalar']), 
                                       ('vector', vec_data, colors['vector']), 
                                       ('tensor', tens_data, colors['tensor'])]:
            if data is not None and col_indices[data_type][i] is not None:
                col_idx = col_indices[data_type][i]
                if data.shape[1] > col_idx:
                    l = data[:, 0]
                    cl = data[:, col_idx]
                    
                    if spectrum == 'BB' and data_type == 'scalar':
                        continue  # Scalar doesn't contribute to BB
                    
                    if spectrum in ['TT', 'EE']:
                        # Log-log plot for TT and EE
                        ax.loglog(l, cl * normalization, 
                                 color=color, linewidth=2, label=data_type.capitalize())
                    elif spectrum == 'BB':
                        # Log-log for BB, but handle positive values only
                        mask = cl > 0
                        if np.any(mask):
                            ax.loglog(l[mask], cl[mask] * normalization, 
                                     color=color, linewidth=2, label=data_type.capitalize())
                    else:  # TE
                        # Semi-log for TE (can be negative)
                        ax.semilogx(l, cl * normalization, 
                                   color=color, linewidth=2, label=data_type.capitalize())
        
        ax.set_xlabel(r'$\ell$', fontsize=14)
        ylabel = f'$\\ell(\\ell+1)C_\\ell^{{{spectrum}}} \\cdot G\\mu^2/(2\\pi)$ [$\\mu$K$^2$]'
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(f'{spectrum} Power Spectrum - {root_name}', fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=12)
        
        if spectrum == 'TE':
            ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        
        # Save individual plot
        output_file = os.path.join(output_dir, f'{root_name}_{spectrum}_spectrum.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Individual plot saved as: {output_file}")

def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) > 1:
        root_name = sys.argv[1]
    else:
        root_name = 'test'
    
    print(f"Plotting CMB power spectra for root name: {root_name}")
    plot_cls(root_name)

if __name__ == '__main__':
    main() 
