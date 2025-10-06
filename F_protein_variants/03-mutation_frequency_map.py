#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches # For drawing site bars
import matplotlib.lines as mlines # For antibody sites
from matplotlib.colors import BoundaryNorm, ListedColormap # For custom color mapping
import sys
from pathlib import Path
import re # For sorting
import numpy as np # For finding unique tick positions

import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'  # Ensures text remains as text

# Set font to Arial
matplotlib.rcParams['font.family'] = 'Arial'

sns.set_context("notebook", font_scale=1.1) # Adjusted base font scale slightly
sns.set_style("ticks")

# --- Configuration ---
RESULTS_DIR = Path("results")
INPUT_FREQUENCY_TABLE = RESULTS_DIR / "rsv_F_mutation_frequency_by_country_mafft_extracted_final.tsv"
OUTPUT_HEATMAP_FILE_A = RESULTS_DIR / "rsv_A_mutation_heatmap_UAE_vs_Public_gt1pct_or_gt2strains_gt30samples_v4.svg"
OUTPUT_HEATMAP_FILE_B = RESULTS_DIR / "rsv_B_mutation_heatmap_UAE_vs_Public_gt1pct_or_gt2strains_gt30samples_v4.svg"
# Per-type UAE frequency thresholds
UAE_FREQUENCY_THRESHOLD_A = 0.01  # >1% for RSV-A
UAE_FREQUENCY_THRESHOLD_B = 0.01  # >1% for RSV-B
MIN_UAE_STRAIN_COUNT = 2          # Keep at least 3 UAE strains to avoid singletons
MIN_SAMPLE_THRESHOLD = 30         # Minimum samples for public data countries to be included
ROW_MIN_FREQUENCY_FILTER_THRESHOLD = 0.95 # Filter out mutations if all countries have freq > this
F_PROTEIN_LENGTH = 574 # Max length of F protein for the site axis
PUBLIC_DATA_YEAR_RANGE = "2021-2024" # For plot title

# --- Helper function for sorting mutations ---
def get_pos_from_mutation(mut_str):
    """Extracts position number from mutation string for sorting."""
    match = re.search(r'\d+', str(mut_str))
    return int(match.group()) if match else 0

# --- Function to generate and save heatmap ---
def generate_heatmap(data_df, virus_type, output_filename):
    """Pivots data, generates and saves a heatmap with site annotations."""
    print(f"\n--- Generating Heatmap for RSV-{virus_type} ---")
    print(f"Data passed to generate_heatmap for RSV-{virus_type} has {len(data_df)} rows.", file=sys.stderr)
    if not data_df.empty:
        print(f"Countries in data_df for RSV-{virus_type} before pivot: {data_df['Country'].unique()}", file=sys.stderr)
    else:
        print(f"No data available for RSV-{virus_type} after filtering. Skipping heatmap.", file=sys.stderr)
        return

    # Get total samples per country for x-tick labels before pivoting
    country_totals = data_df.groupby('Country')['Total_Samples'].first().to_dict()

    # Pivot the data for the heatmap
    print(f"Pivoting data for RSV-{virus_type} heatmap...")
    try:
        heatmap_data = data_df.pivot_table(
            index='Mutation',
            columns='Country',
            values='Frequency',
            fill_value=0
        )
    except Exception as e:
        print(f"ERROR: Could not pivot data for RSV-{virus_type}: {e}", file=sys.stderr)
        return

    if heatmap_data.empty or heatmap_data.isnull().all().all():
        print(f"Pivoted data for RSV-{virus_type} is empty or all NaN. Skipping heatmap.", file=sys.stderr)
        return

    # Sort the heatmap rows (mutations) by position
    heatmap_data['Pos'] = heatmap_data.index.map(get_pos_from_mutation)
    heatmap_data = heatmap_data.sort_values(by='Pos').drop(columns=['Pos'])

    # Ensure UAE is the first column if present
    if 'UAE' in heatmap_data.columns:
        uae_col = heatmap_data.pop('UAE')
        heatmap_data.insert(0, 'UAE', uae_col)

    # Filter out mutations where all country frequencies are very high (e.g., > 0.95)
    if not heatmap_data.empty:
        rows_to_drop = heatmap_data[heatmap_data.gt(ROW_MIN_FREQUENCY_FILTER_THRESHOLD).all(axis=1)].index
        if not rows_to_drop.empty:
            print(f"Filtering out {len(rows_to_drop)} mutations where all country frequencies are > {ROW_MIN_FREQUENCY_FILTER_THRESHOLD}", file=sys.stderr)
            heatmap_data = heatmap_data.drop(index=rows_to_drop)

    if heatmap_data.empty:
        print(f"Heatmap data became empty after filtering high-frequency rows for RSV-{virus_type}. Skipping heatmap.", file=sys.stderr)
        return

    print(f"Heatmap data prepared for RSV-{virus_type} with {heatmap_data.shape[0]} mutations and {heatmap_data.shape[1]} countries.", file=sys.stderr)

    print(f"Generating heatmap for RSV-{virus_type}...")
    fig_height = max(8, heatmap_data.shape[0] * 0.45) + 2
    fig_width = max(5, heatmap_data.shape[1] * 0.65 + 3.5)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    boundaries = [0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.80, 0.85, 0.90, 0.95, 1.01]
    num_colors = len(boundaries) - 1
    blues_cmap_base = matplotlib.colormaps.get_cmap('Blues')

    color_list = []
    discrete_colors_end_proportion = 0.7
    num_discrete_bins_low_freq = 7

    for i in range(num_discrete_bins_low_freq):
        color_list.append(blues_cmap_base( (i / (num_discrete_bins_low_freq)) * discrete_colors_end_proportion ) )

    num_gradient_bins_high_freq = 4
    for i in range(num_gradient_bins_high_freq):
        color_list.append(blues_cmap_base( discrete_colors_end_proportion + ( (i+1) / num_gradient_bins_high_freq ) * (1-discrete_colors_end_proportion) ) )

    if len(color_list) != num_colors:
        print("Warning: Color list length mismatch. Using default linspace sampling.", file=sys.stderr)
        color_list = [blues_cmap_base(x) for x in np.linspace(0.1, 1, num_colors)]

    custom_cmap = ListedColormap(color_list)
    norm = BoundaryNorm(boundaries, custom_cmap.N)

    cbar_left = 0.88
    cbar_bottom = 0.20
    cbar_width = 0.015
    cbar_height = 0.7
    cbar_ax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])

    sns.heatmap(
        heatmap_data,
        annot=False,
        fmt=".2f",
        cmap=custom_cmap,
        norm=norm,
        linewidths=.5,
        linecolor='lightgray',
        cbar_ax=cbar_ax,
        cbar_kws={'label': 'Mutation Frequency', 'ticks': boundaries[:-1], 'spacing': 'proportional', 'format': '%.2f'},
        ax=ax
    )
    cbar_ax.set_yticklabels([f'{b*100:.0f}%' for b in boundaries[:-1]])
    cbar_ax.yaxis.label.set_size(14)
    cbar_ax.tick_params(labelsize=12)

    # Updated title to reflect new UAE filtering criteria (per-type thresholds)
    uae_filter_desc = f"UAE: >{UAE_FREQUENCY_THRESHOLD_A*100:.0f}% (A) / >{UAE_FREQUENCY_THRESHOLD_B*100:.0f}% (B) and ≥{MIN_UAE_STRAIN_COUNT} strains"
    public_filter_desc = f"Public N > {MIN_SAMPLE_THRESHOLD}"
    row_filter_desc = f"Not all Freq > {ROW_MIN_FREQUENCY_FILTER_THRESHOLD*100:.0f}%"
    full_filter_description = f"({uae_filter_desc}; {public_filter_desc}; {row_filter_desc})"

    ax.set_title(f'RSV-{virus_type} F Mutation Frequency: UAE vs Public Data ({PUBLIC_DATA_YEAR_RANGE})\n{full_filter_description}', fontsize=16, loc='center', pad=20)
    ax.set_xlabel('Country/Region', fontsize=16)
    ax.set_ylabel('Mutation', fontsize=16)

    current_xticklabels_texts = [label.get_text() for label in ax.get_xticklabels()]
    new_xticklabels = []
    for country_name in current_xticklabels_texts:
        total_n = country_totals.get(country_name, 'N/A')
        new_xticklabels.append(f"{country_name} (N={total_n})")

    ax.set_xticklabels(new_xticklabels, rotation=90, ha='right', rotation_mode='anchor')
    ax.tick_params(axis='x', labelsize=24)
    ax.tick_params(axis='y', labelsize=20, rotation=0)

    # --- Add Site Annotation Panel ---
    site_panel_left = 0.76
    site_panel_width = 0.08
    site_panel_bottom = 0.20
    site_panel_height = 0.75

    ax_sites = fig.add_axes([site_panel_left, site_panel_bottom, site_panel_width, site_panel_height])
    ax_sites.set_ylim(F_PROTEIN_LENGTH + 1, 0)
    ax_sites.set_xlim(0, 1)
    ax_sites.set_xticks([])
    ax_sites.set_yticks([])
    ax_sites.axis('off')

    antigenic_sites_data = {
        "Site Ø": {"ranges": [(62, 96), (195, 227)], "color": "#FFB6C1"},
        "Site I": {"ranges": [(27, 45), (312, 318), (378, 389)], "color": "#ADD8E6"},
        "Site II": {"ranges": [(254, 277)], "color": "#90EE90"},
        "Site III": {"ranges": [(46, 54), (301, 311), (345, 352), (367, 378)], "color": "#FFD700"},
        "Site IV": {"ranges": [(422, 471)], "color": "#DA70D6"},
        "Site V": {"ranges": [(55, 61), (146, 194), (287, 300)], "color": "#FAA460"},
        "P27 (FP)": {"ranges": [(110, 136)], "color": "#778899"},
    }
    antibody_sites_data = {
        "Nirsevimab": {"ranges": [(62, 69), (196, 212)], "color": "darkviolet", "linestyle":'-', "linewidth": 2.5, "label": "Nirsevimab", "label_fontsize": 12, "label_x_offset": 0.65},
        "Palivizumab": {"ranges": [(262, 275)], "color": "darkcyan", "linestyle":'--', "linewidth": 2.5, "label": "Palivizumab", "label_fontsize": 12, "label_x_offset": 0.65},
    }

    site_bar_x_start = 0.05
    site_bar_width = 0.3
    aa_label_x = site_bar_x_start + site_bar_width + 0.35

    all_site_boundaries = {1, F_PROTEIN_LENGTH}
    legend_handles_antigenic = []

    for site_name, site_info in antigenic_sites_data.items():
        color = site_info["color"]
        legend_handles_antigenic.append(patches.Patch(facecolor=color, edgecolor='black', label=site_name, linewidth=0.5))

        for start_aa, end_aa in site_info["ranges"]:
            all_site_boundaries.add(start_aa)
            all_site_boundaries.add(end_aa)
            rect = patches.Rectangle((site_bar_x_start, start_aa - 0.5), site_bar_width, (end_aa - start_aa) + 1,
                                     linewidth=0.5, edgecolor='black', facecolor=color, alpha=0.7,
                                     clip_on=False)
            ax_sites.add_patch(rect)

    antibody_line_x = site_bar_x_start + site_bar_width + 0.05
    for ab_name, ab_info in antibody_sites_data.items():
        label_y_positions = []
        for r_start, r_end in ab_info["ranges"]:
            all_site_boundaries.add(r_start)
            all_site_boundaries.add(r_end)
            line = mlines.Line2D([antibody_line_x, antibody_line_x + 0.1], [r_start, r_start],
                                 color=ab_info["color"], linestyle=ab_info["linestyle"],
                                 linewidth=ab_info["linewidth"], clip_on=False)
            ax_sites.add_line(line)
            line = mlines.Line2D([antibody_line_x, antibody_line_x + 0.1], [r_end, r_end],
                                 color=ab_info["color"], linestyle=ab_info["linestyle"],
                                 linewidth=ab_info["linewidth"], clip_on=False)
            ax_sites.add_line(line)
            line = mlines.Line2D([antibody_line_x + 0.05, antibody_line_x + 0.05], [r_start, r_end],
                                 color=ab_info["color"], linestyle=ab_info["linestyle"],
                                 linewidth=ab_info["linewidth"], clip_on=False)
            ax_sites.add_line(line)
            label_y_positions.append((r_start + r_end) / 2)

        avg_y_pos = sum(label_y_positions) / len(label_y_positions) if label_y_positions else 0
        ax_sites.text(ab_info["label_x_offset"], avg_y_pos, ab_info["label"],
                      horizontalalignment='left', verticalalignment='center',
                      fontsize=ab_info["label_fontsize"], color=ab_info["color"], fontweight='bold', clip_on=False)

    sorted_boundaries = sorted(list(all_site_boundaries))
    tick_positions = [1]
    for i in range(1, len(sorted_boundaries)):
        if sorted_boundaries[i] - tick_positions[-1] >= 20:
            tick_positions.append(sorted_boundaries[i])
    if F_PROTEIN_LENGTH not in tick_positions:
         tick_positions.append(F_PROTEIN_LENGTH)
    tick_positions = sorted(list(set(tick_positions)))

    for aa_pos in tick_positions:
        ax_sites.text(aa_label_x, aa_pos, str(aa_pos),
                      horizontalalignment='left', verticalalignment='center',
                      fontsize=10, clip_on=False)

    ax_sites.set_title("F Protein Sites", fontsize=14, loc='center', pad=10)

    fig.legend(handles=legend_handles_antigenic,
               loc='lower center',
               bbox_to_anchor=(1.1, 0.02), # Adjusted bbox_to_anchor; may need fine-tuning
               ncol=len(legend_handles_antigenic) // 2 or 1,
               fontsize=12, title="Antigenic Sites", title_fontsize=14)

    plt.subplots_adjust(left=0.1, right=0.72, top=0.92, bottom=0.22, wspace=0.5)

    try:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Heatmap for RSV-{virus_type} saved successfully to: {output_filename}", file=sys.stderr)
    except Exception as e:
        print(f"ERROR: Could not save heatmap file {output_filename}: {e}", file=sys.stderr)
    plt.close(fig)


# --- Main Script ---
print(f"Loading frequency data from: {INPUT_FREQUENCY_TABLE}")
try:
    freq_df = pd.read_csv(INPUT_FREQUENCY_TABLE, sep='\t')
    print(f"Successfully loaded {len(freq_df)} frequency records.", file=sys.stderr)
except FileNotFoundError:
    print(f"ERROR: Input frequency file not found: {INPUT_FREQUENCY_TABLE}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"ERROR: Could not read input file {INPUT_FREQUENCY_TABLE}: {e}", file=sys.stderr)
    sys.exit(1)

required_cols = ['Mutation', 'Country', 'Frequency', 'Type', 'Total_Samples']
if not all(col in freq_df.columns for col in required_cols):
    print(f"ERROR: Input file is missing one or more required columns: {required_cols}", file=sys.stderr)
    print(f"Available columns: {freq_df.columns.tolist()}", file=sys.stderr)
    sys.exit(1)

# Convert relevant columns to numeric, handling potential errors
for col in ['Frequency', 'Total_Samples']:
    freq_df[col] = pd.to_numeric(freq_df[col], errors='coerce')

# Drop rows where Frequency or Total_Samples could not be converted (became NaN)
# or where Total_Samples is zero, as it would make strain count calculation problematic.
initial_rows = len(freq_df)
freq_df.dropna(subset=['Frequency', 'Total_Samples'], inplace=True)
freq_df = freq_df[freq_df['Total_Samples'] > 0]
if len(freq_df) < initial_rows:
    print(f"Warning: Dropped {initial_rows - len(freq_df)} rows due to missing/invalid Frequency or Total_Samples, or zero Total_Samples.", file=sys.stderr)


print(f"Identifying mutations with frequency > {UAE_FREQUENCY_THRESHOLD_A*100:.0f}% (A) / {UAE_FREQUENCY_THRESHOLD_B*100:.0f}% (B) AND found in at least {MIN_UAE_STRAIN_COUNT} strains in UAE...")
uae_df = freq_df[freq_df['Country'] == 'UAE'].copy()

if uae_df.empty:
    print("No data found for UAE. Cannot identify mutations for filtering.", file=sys.stderr)
    uae_high_freq_mutations = [] # Ensure it's an empty list, not undefined
else:
    # Calculate the number of UAE strains with each mutation
    uae_df['Num_UAE_Strains_With_Mutation'] = uae_df['Frequency'] * uae_df['Total_Samples']

    # Per-type frequency threshold
    freq_thresh_map = {'A': UAE_FREQUENCY_THRESHOLD_A, 'B': UAE_FREQUENCY_THRESHOLD_B}
    uae_df['FreqThresh'] = uae_df['Type'].map(freq_thresh_map).fillna(1.0)

    # Condition 1: Frequency > per-type threshold
    cond1_freq = uae_df['Frequency'] > uae_df['FreqThresh']
    # Condition 2: At least MIN_UAE_STRAIN_COUNT UAE strains have the mutation
    cond2_count = uae_df['Num_UAE_Strains_With_Mutation'] >= MIN_UAE_STRAIN_COUNT

    # Combine conditions (AND)
    uae_selected_mutations_df = uae_df[cond1_freq & cond2_count]
    uae_high_freq_mutations = uae_selected_mutations_df['Mutation'].unique()

if len(uae_high_freq_mutations) == 0:
    print(f"No mutations found in UAE with frequency > {UAE_FREQUENCY_THRESHOLD*100:.0f}% OR in at least {MIN_UAE_STRAIN_COUNT} strains. Cannot generate heatmaps.", file=sys.stderr)
    sys.exit(0) # Exit gracefully if no mutations are selected

print(f"Found {len(uae_high_freq_mutations)} mutations meeting the criteria in UAE (combined A and B).", file=sys.stderr)
# print(f"Selected UAE mutations: {list(uae_high_freq_mutations)}", file=sys.stderr) # Optional: for debugging

print("Filtering data for selected mutations across all countries...")
filtered_df_mutations = freq_df[freq_df['Mutation'].isin(uae_high_freq_mutations)].copy()

if filtered_df_mutations.empty:
     print(f"Filtered dataframe is empty after selecting mutations. Check mutation list and input file.", file=sys.stderr)
     sys.exit(1)
print(f"Data filtered for {len(uae_high_freq_mutations)} mutations. Rows: {len(filtered_df_mutations)}", file=sys.stderr)

print(f"Identifying Country/Type combinations with > {MIN_SAMPLE_THRESHOLD} total samples...")
# Use max() for Total_Samples as it should be consistent for a given Country/Type
country_type_totals = freq_df.groupby(['Country', 'Type'])['Total_Samples'].max()

countries_to_keep = country_type_totals[country_type_totals > MIN_SAMPLE_THRESHOLD].index.tolist()
countries_types_set = set(countries_to_keep)
# Ensure UAE is always considered for plotting if it has data for the selected mutations,
# regardless of its Total_Samples meeting MIN_SAMPLE_THRESHOLD (which is for public data)
countries_types_set.add(('UAE', 'A'))
countries_types_set.add(('UAE', 'B'))
print(f"Found {len(countries_types_set)} Country/Type combinations meeting sample threshold (or is UAE).", file=sys.stderr)

print("Filtering data based on country sample size (or if country is UAE)...")
final_filtered_df = filtered_df_mutations[
    filtered_df_mutations.apply(lambda row: (row['Country'], row['Type']) in countries_types_set, axis=1)
].copy()

if final_filtered_df.empty:
     print(f"Filtered dataframe is empty after filtering by country sample size. Check thresholds and input file.", file=sys.stderr)
     sys.exit(1)
print(f"Data filtered down to {len(final_filtered_df)} rows after country sample size filter.", file=sys.stderr)

df_a = final_filtered_df[final_filtered_df['Type'] == 'A'].copy()
df_b = final_filtered_df[final_filtered_df['Type'] == 'B'].copy()

generate_heatmap(df_a, 'A', OUTPUT_HEATMAP_FILE_A)
generate_heatmap(df_b, 'B', OUTPUT_HEATMAP_FILE_B)

print("\nScript finished.")
