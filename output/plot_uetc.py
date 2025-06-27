#!/usr/bin/env python3
"""
Plot UETC (Unequal Time Correlator) data from CAMB output files.

This script reads UETC data files and creates plots showing the correlators
as functions of k*tau values.
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

def load_uetc_matrix(filename):
    """Load UETC matrix data from file."""
    try:
        # Read all values from the file, ignoring line structure
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
        
        # Determine matrix size based on number of values
        n_total = len(values)
        n_size = int(np.sqrt(n_total))
        
        if n_size * n_size != n_total:
            print(f"Warning: {filename} has {n_total} values, not a perfect square")
            # Try to find the closest square
            n_size = int(np.sqrt(n_total))
            values = values[:n_size*n_size]
        
        # Reshape into matrix (assuming Fortran column-major order)
        matrix = values.reshape((n_size, n_size), order='F')
        
        return matrix
        
    except Exception as e:
        print(f"Error loading UETC data from {filename}: {e}")
        return None

def plot_uetc_1d(ktau, uetc_data, title, ax, normalization=1.0):
    """Plot 1D UETC data."""
    uetc_normalized = uetc_data * normalization
    ax.loglog(ktau, np.abs(uetc_normalized), 'b-', linewidth=1.5, alpha=0.8)
    ax.set_xlabel('k τ')
    ax.set_ylabel(f'|UETC| · Gμ²')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(ktau[0], ktau[-1])

def plot_uetc_2d(ktau, uetc_matrix, title, ax, normalization=1.0):
    """Plot 2D UETC matrix data as heatmap."""
    uetc_normalized = uetc_matrix * normalization
    
    if uetc_matrix.ndim == 2:
        # Create k*tau grid for 2D plot
        ktau_grid_x, ktau_grid_y = np.meshgrid(ktau, ktau)
        
        # Plot as contour/heatmap
        im = ax.contourf(ktau_grid_x, ktau_grid_y, uetc_normalized, levels=50, cmap='viridis')
        ax.set_xlabel('k τ₁')
        ax.set_ylabel('k τ₂')
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.colorbar(im, ax=ax, label=f'UETC · Gμ²')
    else:
        # Fall back to 1D plot
        plot_uetc_1d(ktau, uetc_matrix, title, ax, normalization)
    
    ax.set_title(title)

def plot_uetc_diagonal(ktau, uetc_matrix, title, ax, normalization=1.0):
    """Plot diagonal elements of UETC matrix."""
    if uetc_matrix.ndim == 2 and uetc_matrix.shape[0] == uetc_matrix.shape[1]:
        diagonal = np.diag(uetc_matrix) * normalization
        ax.loglog(ktau, np.abs(diagonal), 'r-', linewidth=2, label='Diagonal')
        ax.set_xlabel('k τ')
        ax.set_ylabel(f'|UETC diagonal| · Gμ²')
        ax.set_title(f'{title} (Diagonal)')
        ax.grid(True, alpha=0.3)
        ax.legend()
    else:
        plot_uetc_1d(ktau, uetc_matrix, title, ax, normalization)

def plot_uetc_comparison(ktau, uetc_data, output_dir='plots', root_name='test'):
    """Create comparison plot of all UETC types."""
    
    # Apply normalization factor Gμ² = (2×10⁻⁷)²
    #normalization = (2e-7)**2
    normalization = 1.0
    print(f"Using normalization factor Gμ² = {normalization:.2e}")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('UETC Comparison', fontsize=16)
    
    # Define plot parameters for each UETC type
    uetc_types = ['ss00', 'ss', 'tt']
    colors = ['blue', 'red', 'green']
    labels = ['SS00', 'SS', 'TT']
    
    # Plot diagonal elements (i=j) for each matrix
    for i, (uetc_type, color, label) in enumerate(zip(uetc_types, colors, labels)):
        if uetc_type in uetc_data and uetc_data[uetc_type] is not None:
            matrix = uetc_data[uetc_type] * normalization
            diagonal = np.diag(matrix)
            
            # Plot in first subplot
            axes[0, 0].loglog(ktau, np.abs(diagonal), color=color, label=f'{label} diagonal', linewidth=2)
    
    axes[0, 0].set_xlabel('k τ')
    axes[0, 0].set_ylabel('|UETC(k τ, k τ)| · Gμ²')
    axes[0, 0].set_title('Diagonal Elements')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()
    
    # Plot off-diagonal slice (fixed i, varying j)
    fixed_i = len(ktau) // 4  # Use 1/4 point along k*tau
    for i, (uetc_type, color, label) in enumerate(zip(uetc_types, colors, labels)):
        if uetc_type in uetc_data and uetc_data[uetc_type] is not None:
            matrix = uetc_data[uetc_type] * normalization
            slice_data = matrix[fixed_i, :]
            
            axes[0, 1].semilogx(ktau, slice_data, color=color, label=f'{label} slice', linewidth=2)
    
    axes[0, 1].set_xlabel('k τ')
    axes[0, 1].set_ylabel(f'UETC(k τ_{{{fixed_i}}}, k τ) · Gμ²')
    axes[0, 1].set_title(f'Cross-correlation (fixed k τ_{{{fixed_i}}})')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    
    # Plot 2D heatmap for SS00 if available
    if 'ss00' in uetc_data and uetc_data['ss00'] is not None:
        matrix = uetc_data['ss00'] * normalization
        n_points = len(ktau)
        im = axes[1, 0].imshow(matrix, aspect='auto', origin='lower', cmap='RdBu_r')
        axes[1, 0].set_xlabel('k τ')
        axes[1, 0].set_ylabel('k τ')
        axes[1, 0].set_title('SS00 UETC Matrix')
        
        # Set custom tick locations and labels
        tick_indices = [0, n_points//4, n_points//2, 3*n_points//4, n_points-1]
        tick_labels = [f'{ktau[i]:.2f}' for i in tick_indices]
        axes[1, 0].set_xticks(tick_indices)
        axes[1, 0].set_xticklabels(tick_labels)
        axes[1, 0].set_yticks(tick_indices)
        axes[1, 0].set_yticklabels(tick_labels)
        
        plt.colorbar(im, ax=axes[1, 0], label='UETC · Gμ²')
    
    # Plot correlation function vs separation
    if 'ss00' in uetc_data and uetc_data['ss00'] is not None:
        matrix = uetc_data['ss00'] * normalization
        n_points = len(ktau)
        separations = []
        correlations = []
        
        # Calculate correlation as function of separation |i-j|
        for sep in range(0, n_points//4):  # Only go to 1/4 of the range
            corr_values = []
            for i in range(n_points - sep):
                if i + sep < n_points:
                    corr_values.append(matrix[i, i + sep])
            if corr_values:
                separations.append(sep)
                correlations.append(np.mean(corr_values))
        
        if separations:
            axes[1, 1].semilogy(separations, np.abs(correlations), 'b-', linewidth=2, label='SS00')
            axes[1, 1].set_xlabel('Index separation |i - j|')
            axes[1, 1].set_ylabel('|⟨UETC⟩| · Gμ²')
            axes[1, 1].set_title('Correlation vs Separation')
            axes[1, 1].grid(True, alpha=0.3)
            axes[1, 1].legend()
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{root_name}_uetc_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved {output_file}")
    
    # Also save as PDF
    output_file_pdf = os.path.join(output_dir, f'{root_name}_uetc_comparison.pdf')
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved {output_file_pdf}")
    
    plt.show()

def plot_individual_uetc(ktau, matrix, uetc_type, output_dir='plots', root_name='test'):
    """Create individual plot for a specific UETC type."""
    
    # Apply normalization factor Gμ² = (2×10⁻⁷)²
    #normalization = (2e-7)**2
    normalization = 1.0
    matrix_norm = matrix * normalization
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'UETC {uetc_type.upper()} Analysis', fontsize=16)
    
    # 1. Diagonal elements
    diagonal = np.diag(matrix_norm)
    axes[0, 0].loglog(ktau, np.abs(diagonal), 'b-', linewidth=2)
    axes[0, 0].set_xlabel('k τ')
    axes[0, 0].set_ylabel('|UETC(k τ, k τ)| · Gμ²')
    axes[0, 0].set_title('Diagonal Elements')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Several cross-sections
    n_points = len(ktau)
    indices = [n_points//8, n_points//4, n_points//2, 3*n_points//4]
    colors = ['red', 'green', 'blue', 'orange']
    
    for idx, color in zip(indices, colors):
        if idx < n_points:
            slice_data = matrix_norm[idx, :]
            axes[0, 1].semilogx(ktau, slice_data, color=color, 
                               label=f'k τ = {ktau[idx]:.2f}', linewidth=2)
    
    axes[0, 1].set_xlabel('k τ')
    axes[0, 1].set_ylabel('UETC · Gμ²')
    axes[0, 1].set_title('Cross-correlations')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    
    # 3. 2D heatmap
    n_points = len(ktau)
    im = axes[1, 0].imshow(matrix_norm, aspect='auto', origin='lower', cmap='RdBu_r')
    axes[1, 0].set_xlabel('k τ')
    axes[1, 0].set_ylabel('k τ')
    axes[1, 0].set_title('Full UETC Matrix')
    
    # Set custom tick locations and labels
    tick_indices = [0, n_points//4, n_points//2, 3*n_points//4, n_points-1]
    tick_labels = [f'{ktau[i]:.2f}' for i in tick_indices]
    axes[1, 0].set_xticks(tick_indices)
    axes[1, 0].set_xticklabels(tick_labels)
    axes[1, 0].set_yticks(tick_indices)
    axes[1, 0].set_yticklabels(tick_labels)
    
    plt.colorbar(im, ax=axes[1, 0], label='UETC · Gμ²')
    
    # 4. Correlation function vs separation
    separations = []
    correlations = []
    
    for sep in range(0, n_points//4):
        corr_values = []
        for i in range(n_points - sep):
            if i + sep < n_points:
                corr_values.append(matrix_norm[i, i + sep])
        if corr_values:
            separations.append(sep)
            correlations.append(np.mean(corr_values))
    
    if separations:
        axes[1, 1].semilogy(separations, np.abs(correlations), 'b-', linewidth=2)
        axes[1, 1].set_xlabel('Index separation |i - j|')
        axes[1, 1].set_ylabel('|⟨UETC⟩| · Gμ²')
        axes[1, 1].set_title('Correlation vs Separation')
        axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{root_name}_uetc_{uetc_type}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved {output_file}")
    
    # Also save as PDF
    output_file_pdf = os.path.join(output_dir, f'{root_name}_uetc_{uetc_type}.pdf')
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved {output_file_pdf}")
    
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot UETC data from CAMB output')
    parser.add_argument('--root-name', default='test', 
                       help='Root name for data files (default: test)')
    parser.add_argument('--data-dir', default='data', 
                       help='Directory containing data files (default: data)')
    parser.add_argument('--output-dir', default='plots', 
                       help='Directory for output plots (default: plots)')
    parser.add_argument('--individual', action='store_true',
                       help='Create individual plots for each UETC type')
    
    args = parser.parse_args()
    
    # File paths
    ktau_file = os.path.join(args.data_dir, f'{args.root_name}_uetc_ktau.dat')
    
    # Load k*tau data
    print(f"Loading k*tau data from {ktau_file}")
    ktau = load_ktau_data(ktau_file)
    if ktau is None:
        print("Failed to load k*tau data")
        return 1
    
    print(f"Loaded {len(ktau)} k*tau values from {ktau_file}")
    
    # Load UETC data
    uetc_types = ['ss00', 'ss', 'sscross', 'vv', 'tt']
    uetc_data = {}
    
    for uetc_type in uetc_types:
        uetc_file = os.path.join(args.data_dir, f'{args.root_name}_uetc_{uetc_type}.dat')
        if os.path.exists(uetc_file):
            print(f"Loading UETC data from {uetc_file}")
            matrix = load_uetc_matrix(uetc_file)
            if matrix is not None:
                uetc_data[uetc_type] = matrix
                print(f"Loaded {matrix.shape[0]}×{matrix.shape[1]} UETC matrix for {uetc_type}")
            else:
                print(f"Failed to load UETC data from {uetc_file}")
        else:
            print(f"UETC file not found: {uetc_file}")
    
    if not uetc_data:
        print("No UETC data loaded successfully")
        return 1
    
    # Create comparison plot
    plot_uetc_comparison(ktau, uetc_data, args.output_dir, args.root_name)
    
    # Create individual plots if requested
    if args.individual:
        for uetc_type, matrix in uetc_data.items():
            plot_individual_uetc(ktau, matrix, uetc_type, args.output_dir, args.root_name)
    
    print(f"\nAll UETC plots saved to {args.output_dir}/")
    
    return 0

if __name__ == '__main__':
    sys.exit(main()) 
