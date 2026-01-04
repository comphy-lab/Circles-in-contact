#!/usr/bin/env python3
"""
Generate initial conditions for drop/bubble coalescence simulations.

Creates smooth interface geometries for unequal-sized drops/bubbles
with a fillet transition at the contact point.
"""

import argparse
import os
import subprocess as sp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

plt.rc('text', usetex=True)


def get_circles(delta, Rr):
    """
    Generate interface coordinates for two coalescing drops/bubbles.

    Parameters
    ----------
    delta : float
        Neck radius / minimum length scale. Separation = 2*delta.
    Rr : float
        Radius ratio (large drop / small drop).

    Returns
    -------
    Interface : pd.DataFrame
        DataFrame with 'x' and 'y' columns for interface coordinates.
    X1, Y1 : arrays
        Circle 1 (main drop) coordinates.
    Xf, Yf : arrays
        Fillet circle coordinates.
    X2, Y2 : arrays
        Circle 2 (second drop) coordinates.
    """
    X1c = -(1 + delta)

    # phic is the angle where y coordinate is 2*delta
    phic1 = np.arcsin(2 * delta)
    phi1 = np.linspace(np.pi, phic1, int(1e3))
    X1 = X1c + np.cos(phi1)
    Y1 = np.sin(phi1)

    # Fillet circle
    Yfc = (1 + delta) * np.tan(phic1)
    Rf = (1 + delta) / np.cos(phic1) - 1

    phifStart = np.pi / 2 - phic1
    theta = np.arcsin(Yfc / (Rf + Rr))
    phifEnd = np.pi / 2 - theta

    phif = np.linspace(phifStart, -phifEnd, int(1e3))
    Xf = -Rf * np.sin(phif)
    Yf = Yfc - Rf * np.cos(phif)

    phic2 = np.pi - theta
    X2c = (Rf + Rr) * np.cos(theta)

    phi2 = np.linspace(phic2, 0, int(1e3))
    X2 = X2c + Rr * np.cos(phi2)
    Y2 = Rr * np.sin(phi2)

    # Combine all segments
    X = np.concatenate([X1, Xf, X2])
    Y = np.concatenate([Y1, Yf, Y2])

    Interface = pd.DataFrame({'x': X, 'y': Y})

    return Interface, X1, Y1, Xf, Yf, X2, Y2


def get_facets(L0):
    """
    Run Basilisk to generate facets from the interface data.

    Parameters
    ----------
    L0 : int
        Domain size for Basilisk simulation.

    Returns
    -------
    segs : list
        List of line segments for plotting.
    """
    exe = ["./InitialCondition", str(L0)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if len(temp2) > 1e2:
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
            else:
                if not skip:
                    temp4 = temp2[n1 + 1].split(" ")
                    r1, z1 = float(temp3[1]), float(temp3[0])
                    r2, z2 = float(temp4[1]), float(temp4[0])
                    segs.append(((z1, r1), (z2, r2)))
                    skip = True
    return segs


def plot_interfaces(Interface, X1, Y1, Xf, Yf, X2, Y2, delta, Rr, image_name, show=False):
    """Plot interface without Basilisk comparison."""
    plt.close()
    fig, axs = plt.subplots(2, 2, figsize=(24, 12))
    ax1, ax2, ax3, ax4 = axs.flatten()

    # Window definitions
    NWidth = 10
    xW1, xW2, yW1, yW2 = -NWidth*delta, NWidth*delta, 0, NWidth*delta
    NWidth = 5
    xW3, xW4, yW3, yW4 = -NWidth*delta, NWidth*delta, 0, NWidth*delta
    NWidth = 1
    xW5, xW6, yW5, yW6 = -NWidth*delta, NWidth*delta, delta, NWidth*delta+delta

    # Full view
    ax1.plot(X1, Y1, '-', lw=2, color='#fdae61')
    ax1.plot(Xf, Yf, '-', lw=2, color='#d7191c')
    ax1.plot(X2, Y2, '-', lw=2, color='#2c7bb6')
    ax1.plot([xW1, xW2, xW2, xW1, xW1], [yW1, yW1, yW2, yW2, yW1], 'k-', lw=2)
    ax1.axis('square')
    ax1.set_xlim(-(delta+2)-0.1, delta+2*Rr+0.1)
    ax1.set_ylim(0, 1.1*Rr)
    ax1.set_xlabel(r'$\mathcal{Z}$', fontsize=24)
    ax1.set_ylabel(r'$\mathcal{R}$', fontsize=24)
    ax1.tick_params(axis='both', which='major', labelsize=24)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)

    # Zoom 1
    ax2.plot([xW1, xW2, xW2, xW1, xW1], [yW1, yW1, yW2, yW2, yW1], 'k-', lw=4)
    ax2.plot(X1, Y1, '-', lw=4, color='#fdae61')
    ax2.plot(Xf, Yf, '-', lw=4, color='#d7191c')
    ax2.plot(X2, Y2, '-', lw=4, color='#2c7bb6')
    ax2.plot([xW3, xW4, xW4, xW3, xW3], [yW3, yW3, yW4, yW4, yW3], '-', lw=2, color='gray')
    ax2.axis('square')
    ax2.set_xlim(xW1, xW2)
    ax2.set_ylim(yW1, yW2)
    ax2.axis('off')

    # Zoom 2
    ax3.plot([xW3, xW4, xW4, xW3, xW3], [yW3, yW3, yW4, yW4, yW3], '-', lw=4, color='gray')
    ax3.plot(X1, Y1, '-', lw=4, color='#fdae61')
    ax3.plot(Xf, Yf, '-', lw=4, color='#d7191c')
    ax3.plot(X2, Y2, '-', lw=4, color='#2c7bb6')
    ax3.plot([0, 0], [0, yW4], '-', lw=2, color='gray')
    ax3.plot([xW5, xW6, xW6, xW5, xW5], [yW5, yW5, yW6, yW6, yW5], '-', lw=2, color='#abdda4')
    ax3.axis('square')
    ax3.set_xlim(xW3, xW4)
    ax3.set_ylim(yW3, yW4)
    ax3.axis('off')

    # Zoom 3
    ax4.plot([xW5, xW6, xW6, xW5, xW5], [yW5, yW5, yW6, yW6, yW5], '-', lw=4, color='#abdda4')
    ax4.plot(X1, Y1, '-', lw=4, color='#fdae61')
    ax4.plot(Xf, Yf, '-', lw=4, color='#d7191c')
    ax4.plot(X2, Y2, '-', lw=4, color='#2c7bb6')
    ax4.axis('square')
    ax4.set_xlim(xW5, xW6)
    ax4.set_ylim(yW5, yW6)
    ax4.axis('off')

    plt.savefig(image_name, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()


def plot_interfaces_basilisk(Interface, facets, delta, Rr, image_name, show=False):
    """Plot interface with Basilisk comparison overlay."""
    plt.close()
    fig, axs = plt.subplots(2, 2, figsize=(24, 12))
    ax1, ax2, ax3, ax4 = axs.flatten()

    # Window definitions
    NWidth = 10
    xW1, xW2, yW1, yW2 = -NWidth*delta, NWidth*delta, 0, NWidth*delta
    NWidth = 5
    xW3, xW4, yW3, yW4 = -NWidth*delta, NWidth*delta, 0, NWidth*delta
    NWidth = 1
    xW5, xW6, yW5, yW6 = -NWidth*delta, NWidth*delta, delta, NWidth*delta+delta

    for ax in [ax1, ax2, ax3, ax4]:
        ax.plot(Interface['x'], Interface['y'], '-', lw=4, color='#fdae61')
        line_segments = LineCollection(facets, linewidths=4, colors='#1a9641', linestyle='solid')
        ax.add_collection(line_segments)

    # Full view
    ax1.plot([xW1, xW2, xW2, xW1, xW1], [yW1, yW1, yW2, yW2, yW1], 'k-', lw=2)
    ax1.axis('square')
    ax1.set_xlim(-(delta+2)-0.1, delta+2*Rr+0.1)
    ax1.set_ylim(0, 1.1*Rr)
    ax1.set_xlabel(r'$\mathcal{Z}$', fontsize=24)
    ax1.set_ylabel(r'$\mathcal{R}$', fontsize=24)
    ax1.tick_params(axis='both', which='major', labelsize=24)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)

    # Zoom levels
    ax2.plot([xW1, xW2, xW2, xW1, xW1], [yW1, yW1, yW2, yW2, yW1], 'k-', lw=4)
    ax2.plot([xW3, xW4, xW4, xW3, xW3], [yW3, yW3, yW4, yW4, yW3], '-', lw=2, color='gray')
    ax2.axis('square')
    ax2.set_xlim(xW1, xW2)
    ax2.set_ylim(yW1, yW2)
    ax2.axis('off')

    ax3.plot([xW3, xW4, xW4, xW3, xW3], [yW3, yW3, yW4, yW4, yW3], '-', lw=4, color='gray')
    ax3.plot([0, 0], [0, yW4], '-', lw=2, color='gray')
    ax3.plot([xW5, xW6, xW6, xW5, xW5], [yW5, yW5, yW6, yW6, yW5], '-', lw=2, color='#abdda4')
    ax3.axis('square')
    ax3.set_xlim(xW3, xW4)
    ax3.set_ylim(yW3, yW4)
    ax3.axis('off')

    ax4.plot([xW5, xW6, xW6, xW5, xW5], [yW5, yW5, yW6, yW6, yW5], '-', lw=4, color='#abdda4')
    ax4.axis('square')
    ax4.set_xlim(xW5, xW6)
    ax4.set_ylim(yW5, yW6)
    ax4.axis('off')

    plt.savefig(image_name, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()


def generate_initial_conditions(delta, Rr_list, data_folder, image_folder=None,
                                 run_basilisk=False, L0=8, show_plots=False):
    """
    Generate initial conditions for multiple radius ratios.

    Parameters
    ----------
    delta : float
        Neck radius / minimum length scale.
    Rr_list : list
        List of radius ratios to generate.
    data_folder : str
        Output folder for data files.
    image_folder : str, optional
        Output folder for images. If None, no images are saved.
    run_basilisk : bool
        Whether to run Basilisk for verification.
    L0 : int
        Domain size for Basilisk (if run_basilisk=True).
    show_plots : bool
        Whether to display plots interactively.
    """
    # Create output folders
    os.makedirs(data_folder, exist_ok=True)
    if image_folder:
        os.makedirs(image_folder, exist_ok=True)

    for Rr in Rr_list:
        print(f"Generating Rr = {Rr}, delta = {delta}")

        Interface, X1, Y1, Xf, Yf, X2, Y2 = get_circles(delta, Rr)

        # Save data file
        filename = f"InitialConditionRr-{Rr:.2f}.dat"
        filepath = os.path.join(data_folder, filename)
        Interface.to_csv(filepath, index=False, header=False, sep=' ')
        print(f"  Saved: {filepath}")

        # Also save to f_Testing.dat for Basilisk
        Interface.to_csv('f_Testing.dat', index=False, header=False, sep=' ')

        if image_folder:
            # Plot without Basilisk
            img_name = os.path.join(image_folder, f"Interface_Rr{Rr}.pdf")
            plot_interfaces(Interface, X1, Y1, Xf, Yf, X2, Y2, delta, Rr, img_name, show_plots)

            if run_basilisk:
                # Run Basilisk and plot comparison
                facets = get_facets(L0)
                img_name = os.path.join(image_folder, f"InterfaceBasilisk_Rr{Rr}.pdf")
                plot_interfaces_basilisk(Interface, facets, delta, Rr, img_name, show_plots)


def parse_float_list(s):
    """Parse comma-separated float list."""
    return [float(x.strip()) for x in s.split(',')]


def main():
    parser = argparse.ArgumentParser(
        description='Generate initial conditions for drop/bubble coalescence.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --delta 0.01 --Rr 1,2,4,8
  %(prog)s --delta 0.01 --Rr 1.5,2.0,2.5 --data-folder Data_New --images
  %(prog)s --delta 0.01 --Rr 2 --basilisk --L0 8
        """
    )

    parser.add_argument('--delta', type=float, required=True,
                        help='Neck radius / minimum length scale')
    parser.add_argument('--Rr', type=str, required=True,
                        help='Comma-separated list of radius ratios (e.g., "1,2,4,8")')
    parser.add_argument('--data-folder', type=str, default=None,
                        help='Output folder for data files (default: DataFiles_delta{delta})')
    parser.add_argument('--image-folder', type=str, default=None,
                        help='Output folder for images (default: no images)')
    parser.add_argument('--images', action='store_true',
                        help='Generate images in ImageFiles_delta{delta} folder')
    parser.add_argument('--basilisk', action='store_true',
                        help='Run Basilisk for verification')
    parser.add_argument('--L0', type=int, default=8,
                        help='Domain size for Basilisk (default: 8)')
    parser.add_argument('--show', action='store_true',
                        help='Display plots interactively')

    args = parser.parse_args()

    Rr_list = parse_float_list(args.Rr)

    data_folder = args.data_folder
    if not data_folder:
        data_folder = f"DataFiles_delta{args.delta}"

    image_folder = args.image_folder
    if args.images and not image_folder:
        image_folder = f"ImageFiles_delta{args.delta}"

    generate_initial_conditions(
        delta=args.delta,
        Rr_list=Rr_list,
        data_folder=data_folder,
        image_folder=image_folder,
        run_basilisk=args.basilisk,
        L0=args.L0,
        show_plots=args.show
    )


if __name__ == '__main__':
    main()
