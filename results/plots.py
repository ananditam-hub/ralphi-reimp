import matplotlib.pyplot as plt
import numpy as np

def plot_results(results1, results2, results3, labels, title="Results Plot",
                 xlabel="X-axis", ylabel="Y-axis", ylim=None, save_name="results_plot.png"):
    """
    Plots comparison results of the RALPHI algorithm against other state-of-the-art methods.

    Parameters:
        results1 (list): Result values for dataset 1.
        results2 (list): Result values for dataset 2.
        results3 (list): Result values for dataset 3.
        labels (list): Labels for the x-axis corresponding to each method/tool.
        title (str): Title of the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        ylim (tuple): Y-axis limits (ymin, ymax) to ensure consistent scaling.
        save_name (str): Filename for saving the plot.
    """
    plt.figure(figsize=(10, 6))

    plt.plot(labels, results1, marker='o', linestyle='-', color='g', label='chr14')
    plt.plot(labels, results2, marker='o', linestyle='-', color='r', label='chr22')
    plt.plot(labels, results3, marker='o', linestyle='-', color='b', label='chr4' if len(results3) == 4 else 'chr1')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45)
    if ylim:
        plt.ylim(ylim)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(save_name, dpi=300)
    plt.show()

def main():
    # Define fixed y-axis limits (based on min/max across all datasets)
    y_min = 0
    y_max = 3  # Can be adjusted depending on expected value range

    # 5x coverage results
    chr14 = [0.15, 2.25, 0.35, 0.20]
    chr22 = [0.28, 2.71, 0.44, 0.28]
    chr4  = [0.13, 1.83, 0.28, 0.14]
    labels_5x = ['ralphi', 'WhatsHap', 'LongPhase', 'HapCUT2']

    plot_results(
        results1=chr14,
        results2=chr22,
        results3=chr4,
        labels=labels_5x,
        title="Switch Rate Comparison for 5x Coverage (NA12878: Illumina Platinum Genomes)",
        xlabel="Tools",
        ylabel="Switch Rate",
        ylim=(y_min, y_max),
        save_name="switch_rate_5x.png"
    )

    # 30x coverage results
    chr14_30 = [0.16, 0.97, 0.16]
    chr22_30 = [0.14, 0.71, 0.13]
    chr1_30  = [0.22, 0.89, 0.18]
    labels_30x = ['ralphi', 'WhatsHap', 'HapCUT2']

    plot_results(
        results1=chr14_30,
        results2=chr22_30,
        results3=chr1_30,
        labels=labels_30x,
        title="Switch Rate Comparison for 30x Coverage (NA12878: Illumina Platinum Genomes)",
        xlabel="Tools",
        ylabel="Switch Rate",
        ylim=(y_min, y_max),
        save_name="switch_rate_30x.png"
    )

if __name__ == "__main__":
    main()
