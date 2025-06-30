#!/usr/bin/env python3
"""
Plot eigenvectors from UETC eigenvalue decomposition.

This script reads eigenvector data files and creates plots showing the 
first 10 eigenvectors as functions of k*tau values.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys

def load_ktau_data(filename):
    """Load k*tau values from file."""
    try:
        ktau = np.loadtxt(filename)
        return ktau
    except Exception as e:
        print(f"Error loading k*tau data from {filename}: {e}")
        return None

def load_eigenvector_data(filename):
    """Load eigenvector matrix data from file."""
    try:
        # Read all values from the file, handling Fortran formatting
        with open(filename, 'r') as f:
            content = f.read()
        
        # Split by whitespace and convert to float
        values = []
        for line in content.split('\n'):
            if line.strip():
                # Split each line by whitespace and add to values
                line_values = line.split()
                for val in line_values:
                    if val.strip():
                        try:
                            values.append(float(val))
                        except ValueError:
                            continue
        
        values = np.array(values)
        
        # Determine matrix size - should be n_evec x n_ktau
        # For 256 k*tau points, we expect n_evec rows of 256 values each
        n_ktau = 256  # Based on the ktau file
        n_total = len(values)
        n_evec = n_total // n_ktau
        
        if n_evec * n_ktau != n_total:
            print(f"Warning: {filename} has {n_total} values, expected multiple of {n_ktau}")
            # Truncate to fit
            n_evec = n_total // n_ktau
            values = values[:n_evec * n_ktau]
        
        # Reshape into matrix (each row is an eigenvector)
        evec_matrix = values.reshape((n_evec, n_ktau))
        
        return evec_matrix
        
    except Exception as e:
        print(f"Error loading eigenvector data from {filename}: {e}")
        return None

def plot_eigenvectors(ktau, evec_matrix, n_evec=10, output_dir='plots', root_name='test', evec_type='ss00'):
    """Plot the first n_evec eigenvectors on the same plot."""
    
    # Limit to available eigenvectors
    n_available = min(n_evec, evec_matrix.shape[0])
    print(f"Plotting first {n_available} eigenvectors out of {evec_matrix.shape[0]} available")
    
    # Create single plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # For many eigenvectors, use a colormap and reduce alpha
    if n_available <= 20:
        # Use distinct colors for small number of eigenvectors
        colors = plt.cm.tab20(np.linspace(0, 1, n_available))
        alpha = 0.8
        linewidth = 2
        show_legend = True
    else:
        # Use continuous colormap for many eigenvectors
        colors = plt.cm.viridis(np.linspace(0, 1, n_available))
        alpha = 0.6
        linewidth = 1
        show_legend = False
    
    # Plot each eigenvector
    for i in range(n_available):
        evec = evec_matrix[i, :]
        
        # Plot eigenvector with different color
        if show_legend and i < 10:  # Only show legend for first 10 if showing legend
            label = f'Eigenvector {i+1}'
        else:
            label = None
            
        ax.semilogx(ktau, evec, color=colors[i], linewidth=linewidth, 
                   label=label, alpha=alpha)
    
    ax.set_xlabel('k τ')
    ax.set_ylabel('Eigenvector amplitude')
    ax.set_title(f'First {n_available} Eigenvectors - {evec_type.upper()}')
    ax.grid(True, alpha=0.3)
    
    if show_legend:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        # Add colorbar for many eigenvectors
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                                  norm=plt.Normalize(vmin=1, vmax=n_available))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Eigenvector number')
    
    # Add zero line for reference
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{root_name}_eigenvectors_{evec_type}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved {output_file}")
    
    # Also save as PDF
    output_file_pdf = os.path.join(output_dir, f'{root_name}_eigenvectors_{evec_type}.pdf')
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved {output_file_pdf}")
    
    plt.show()

def plot_eigenvalue_comparison(ktau, evec_matrices, evec_types, n_evec=5, output_dir='plots', root_name='test'):
    """Create comparison plot of eigenvectors from different UETC types."""
    
    fig, axes = plt.subplots(1, n_evec, figsize=(4*n_evec, 6))
    if n_evec == 1:
        axes = [axes]
    
    fig.suptitle(f'Eigenvector Comparison - First {n_evec} Modes', fontsize=16)
    
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    
    for i in range(n_evec):
        ax = axes[i]
        
        for j, (evec_type, evec_matrix) in enumerate(zip(evec_types, evec_matrices)):
            if evec_matrix is not None and i < evec_matrix.shape[0]:
                evec = evec_matrix[i, :]
                color = colors[j % len(colors)]
                ax.semilogx(ktau, evec, color=color, linewidth=2, 
                           label=evec_type.upper(), alpha=0.8)
        
        ax.set_xlabel('k τ')
        ax.set_ylabel(f'Eigenvector {i+1}')
        ax.set_title(f'Mode {i+1}')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        
        if i == 0:  # Only show legend on first subplot
            ax.legend()
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{root_name}_eigenvectors_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved {output_file}")
    
    # Also save as PDF
    output_file_pdf = os.path.join(output_dir, f'{root_name}_eigenvectors_comparison.pdf')
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved {output_file_pdf}")
    
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot eigenvectors from UETC eigenvalue decomposition')
    parser.add_argument('--root-name', default='test_uetc_1', 
                       help='Root name for data files (default: test_uetc_1)')
    parser.add_argument('--data-dir', default='data', 
                       help='Directory containing data files (default: data)')
    parser.add_argument('--output-dir', default='plots', 
                       help='Directory for output plots (default: plots)')
    parser.add_argument('--n-evec', type=int, default=100,
                       help='Number of eigenvectors to plot (default: 100)')
    parser.add_argument('--evec-type', default='ss00',
                       help='Eigenvector type to plot (default: ss00)')
    parser.add_argument('--comparison', action='store_true',
                       help='Create comparison plot of different eigenvector types')
    
    args = parser.parse_args()
    
    # Load k*tau data
    ktau_file = os.path.join(args.data_dir, 'test_uetc_ktau.dat')
    print(f"Loading k*tau data from {ktau_file}")
    ktau = load_ktau_data(ktau_file)
    if ktau is None:
        print("Failed to load k*tau data")
        return 1
    
    print(f"Loaded {len(ktau)} k*tau values from {ktau_file}")
    
    if args.comparison:
        # Load multiple eigenvector types for comparison
        evec_types = ['ss00', 'ss', 'vv', 'tt']  # Add more as needed
        evec_matrices = []
        available_types = []
        
        for evec_type in evec_types:
            evec_file = os.path.join(args.data_dir, f'{args.root_name}_{evec_type}_evec.dat')
            if os.path.exists(evec_file):
                print(f"Loading eigenvector data from {evec_file}")
                evec_matrix = load_eigenvector_data(evec_file)
                if evec_matrix is not None:
                    evec_matrices.append(evec_matrix)
                    available_types.append(evec_type)
                    print(f"Loaded {evec_matrix.shape[0]} eigenvectors of length {evec_matrix.shape[1]} for {evec_type}")
                else:
                    print(f"Failed to load eigenvector data from {evec_file}")
            else:
                print(f"Eigenvector file not found: {evec_file}")
        
        if evec_matrices:
            plot_eigenvalue_comparison(ktau, evec_matrices, available_types, 
                                     min(args.n_evec, 5), args.output_dir, 'test_uetc_1')
        else:
            print("No eigenvector data loaded for comparison")
            return 1
    
    else:
        # Load single eigenvector type
        evec_file = os.path.join(args.data_dir, f'{args.root_name}_{args.evec_type}_evec.dat')
        if not os.path.exists(evec_file):
            print(f"Eigenvector file not found: {evec_file}")
            return 1
        
        print(f"Loading eigenvector data from {evec_file}")
        evec_matrix = load_eigenvector_data(evec_file)
        if evec_matrix is None:
            print(f"Failed to load eigenvector data from {evec_file}")
            return 1
        
        print(f"Loaded {evec_matrix.shape[0]} eigenvectors of length {evec_matrix.shape[1]}")
        
        # Plot eigenvectors
        plot_eigenvectors(ktau, evec_matrix, args.n_evec, args.output_dir, 
                         'test_uetc_1', args.evec_type)
    
    print(f"\nAll eigenvector plots saved to {args.output_dir}/")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
