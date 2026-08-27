#!/bin/bash
# stages/06_visualize.sh — Stage 6: Generate interactive HTML visualizations
#
# Produces: Plotly (sunburst, treemap, pie), Bokeh (rank-abundance curve,
# phylum bar chart, top-taxa heatmap), and Krona (interactive taxonomy wheel).
#
# Usage: bash stages/06_visualize.sh [options]
#   --sample-name NAME  Sample label (must match column name in feature-table.tsv)
#   --output-dir DIR    Directory containing exported-data/ and for HTML output
#   --config FILE       Load configuration from file
#   --help              Show this message
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"

show_help "$0" "$@"

[ -f "$SCRIPT_DIR/../pipeline.conf" ] && load_config "$SCRIPT_DIR/../pipeline.conf"
parse_args "$@"

# Verify inputs
for f in "${OUTPUT_DIR}/exported-data/feature-table.tsv" "${OUTPUT_DIR}/exported-data/taxonomy.tsv"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: $f not found. Run stage 05_export.sh first." >&2
        exit 1
    fi
done

# Check and install dependencies before generating script
install_dependencies python
install_python_packages pandas plotly bokeh

echo "STEP 6: Creating visualizations..."

# --------------------------------------------------------------------------
# Krona setup — install KronaTools if not already available
# --------------------------------------------------------------------------
KRONA_AVAILABLE=true
if ! command -v ktImportText &>/dev/null; then
    echo "  KronaTools not found. Attempting to install via conda..."
    if command -v conda &>/dev/null; then
        if conda install -y -c bioconda krona 2>/dev/null; then
            # KronaTools needs its taxonomy database initialized
            ktUpdateTaxonomy.sh 2>/dev/null || true
        else
            echo "  Warning: could not install KronaTools via conda. Skipping Krona chart." >&2
            KRONA_AVAILABLE=false
        fi
    else
        echo "  Warning: conda not available. Skipping Krona chart." >&2
        echo "  Install manually: conda install -c bioconda krona" >&2
        KRONA_AVAILABLE=false
    fi
fi

PLOT_SCRIPT="${OUTPUT_DIR}/create_plots.py"

# Always regenerate the script to stay in sync with pipeline updates.
# Heredoc uses 'EOF' (single-quoted) so Python variables are NOT shell-expanded.
echo "  Writing Python plot script: $PLOT_SCRIPT"
cat << 'EOF' > "$PLOT_SCRIPT"
import argparse
import os
import sys
import pandas as pd
import plotly.express as px
from bokeh.io import save, output_file
from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, HoverTool, NumeralTickFormatter,
                          BasicTicker, ColorBar, LinearColorMapper)
from bokeh.palettes import Spectral11, Turbo256, Category20
from bokeh.transform import transform
import math

parser = argparse.ArgumentParser(description='Generate eDNA taxonomy visualizations')
parser.add_argument('--sample-name', required=True,
                    help='Sample label as it appears as a column in feature-table.tsv')
parser.add_argument('--data-dir', default='exported-data',
                    help='Path to exported-data directory (default: exported-data)')
parser.add_argument('--output-dir', default='.',
                    help='Directory to write HTML output files (default: current directory)')
args = parser.parse_args()

SAMPLE_COL = args.sample_name
DATA_DIR = args.data_dir
OUT_DIR = args.output_dir

# ==========================================================================
# Load and Merge Data
# ==========================================================================
print("Loading exported data...")
abundances = pd.read_csv(f'{DATA_DIR}/feature-table.tsv', sep='\t', skiprows=1)

if SAMPLE_COL not in abundances.columns:
    available = [c for c in abundances.columns if c != '#OTU ID']
    print(f"Error: sample column '{SAMPLE_COL}' not found in feature-table.tsv", file=sys.stderr)
    print(f"Available sample columns: {available}", file=sys.stderr)
    print(f"Hint: use --sample-name with one of the above values.", file=sys.stderr)
    sys.exit(1)

abundances.rename(columns={'#OTU ID': 'ASV_ID', SAMPLE_COL: 'Abundance'}, inplace=True)
taxonomy = pd.read_csv(f'{DATA_DIR}/taxonomy.tsv', sep='\t')
taxonomy.rename(columns={'Feature ID': 'ASV_ID', 'Taxon': 'Taxonomy'}, inplace=True)
df = pd.merge(abundances, taxonomy, on='ASV_ID')

# ==========================================================================
# Clean and Prepare Taxonomic Data
# ==========================================================================
print("Preparing taxonomic data for plotting...")
tax_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
tax_split = df['Taxonomy'].str.split(';', expand=True)
for i, level in enumerate(tax_levels):
    if i < tax_split.shape[1]:
        df[level] = tax_split[i].str.replace(r'^[dpcofgs]__', '', regex=True).str.strip()
    else:
        df[level] = 'Unassigned'
df.replace('', pd.NA, inplace=True)
for level in tax_levels:
    df[level] = df[level].fillna('Unassigned')

df = df[df['Abundance'] > 0]

if df.empty:
    print("Warning: no ASVs with non-zero abundance found. Check your sample name.", file=sys.stderr)
    sys.exit(1)

total_reads = df['Abundance'].sum()
df['RelativeAbundance'] = df['Abundance'] / total_reads

# ==========================================================================
# Plotly — Sunburst, Treemap, Pie (existing)
# ==========================================================================
print("Generating Plotly sunburst chart...")
sunburst_fig = px.sunburst(
    df, path=['Domain', 'Phylum', 'Class'], values='Abundance',
    title='Hierarchical Taxonomic Composition', color='Phylum')
sunburst_fig.write_html(f'{OUT_DIR}/taxonomic_sunburst.html')

print("Generating Plotly treemap...")
treemap_fig = px.treemap(
    df, path=['Domain', 'Phylum', 'Class'], values='Abundance',
    title='Hierarchical Taxonomic Composition (Treemap View)', color='Phylum')
treemap_fig.write_html(f'{OUT_DIR}/taxonomic_treemap.html')

print("Generating Plotly pie chart...")
phylum_df = df[df['Phylum'] != 'Unassigned']
pie_fig = px.pie(
    phylum_df.groupby('Phylum')['Abundance'].sum().reset_index(),
    names='Phylum', values='Abundance',
    title='Taxonomic Composition at Phylum Level')
pie_fig.write_html(f'{OUT_DIR}/taxonomic_pie_chart.html')

# ==========================================================================
# Bokeh — Rank-Abundance Curve
# ==========================================================================
print("Generating Bokeh rank-abundance curve...")

ranked = df[['ASV_ID', 'Abundance', 'Phylum', 'Class', 'Genus']].copy()
ranked = ranked.sort_values('Abundance', ascending=False).reset_index(drop=True)
ranked['Rank'] = ranked.index + 1
ranked['LogAbundance'] = ranked['Abundance'].apply(lambda x: math.log10(x) if x > 0 else 0)
ranked['Label'] = ranked.apply(
    lambda r: f"{r['Genus']}" if r['Genus'] != 'Unassigned'
    else f"{r['Class']}" if r['Class'] != 'Unassigned'
    else r['Phylum'], axis=1)

source_rank = ColumnDataSource(ranked)

p_rank = figure(
    title="Rank-Abundance Curve (ASVs)",
    x_axis_label="ASV Rank",
    y_axis_label="Read Count (log₁₀ scale)",
    width=900, height=500,
    tools="pan,wheel_zoom,box_zoom,reset,save")
p_rank.scatter('Rank', 'LogAbundance', source=source_rank, size=5,
               color='#2b83ba', alpha=0.7)
p_rank.line('Rank', 'LogAbundance', source=source_rank,
            color='#2b83ba', alpha=0.4)
p_rank.add_tools(HoverTool(tooltips=[
    ("Rank", "@Rank"),
    ("ASV", "@ASV_ID"),
    ("Reads", "@Abundance{0,0}"),
    ("Taxon", "@Label"),
    ("Phylum", "@Phylum"),
]))

output_file(f'{OUT_DIR}/rank_abundance_curve.html',
            title="Rank-Abundance Curve")
save(p_rank)

# ==========================================================================
# Bokeh — Phylum Composition Stacked Bar
# ==========================================================================
print("Generating Bokeh phylum composition bar chart...")

phylum_counts = df.groupby('Phylum')['Abundance'].sum().sort_values(ascending=False)

# Keep top 15 phyla, group the rest as "Other"
TOP_N_PHYLA = 15
top_phyla = phylum_counts.head(TOP_N_PHYLA).index.tolist()
df['PhylumGrouped'] = df['Phylum'].where(df['Phylum'].isin(top_phyla), 'Other')

phylum_summary = df.groupby('PhylumGrouped')['Abundance'].sum().sort_values(ascending=False)
phylum_rel = (phylum_summary / phylum_summary.sum() * 100).reset_index()
phylum_rel.columns = ['Phylum', 'Percentage']

n_phyla = len(phylum_rel)
if n_phyla <= 20:
    palette = Category20[max(3, min(n_phyla, 20))][:n_phyla]
else:
    step = max(1, len(Turbo256) // n_phyla)
    palette = [Turbo256[i * step % len(Turbo256)] for i in range(n_phyla)]

phylum_rel['Color'] = palette
source_bar = ColumnDataSource(phylum_rel)

p_bar = figure(
    x_range=phylum_rel['Phylum'].tolist(),
    title=f"Phylum Composition — {SAMPLE_COL}",
    x_axis_label="Phylum",
    y_axis_label="Relative Abundance (%)",
    width=1000, height=500,
    tools="pan,wheel_zoom,box_zoom,reset,save")
p_bar.vbar(x='Phylum', top='Percentage', width=0.8, source=source_bar,
           color='Color', alpha=0.85)
p_bar.xaxis.major_label_orientation = 0.8
p_bar.add_tools(HoverTool(tooltips=[
    ("Phylum", "@Phylum"),
    ("Abundance (%)", "@Percentage{0.2f}"),
]))

output_file(f'{OUT_DIR}/phylum_composition_bar.html',
            title="Phylum Composition Bar Chart")
save(p_bar)

# ==========================================================================
# Bokeh — Top Taxa Heatmap (Phylum × Class)
# ==========================================================================
print("Generating Bokeh taxonomy heatmap...")

cross = df.groupby(['Phylum', 'Class'])['Abundance'].sum().reset_index()
cross = cross[cross['Abundance'] > 0]

# Top 12 phyla and top 20 classes by total abundance
top_phyla_hm = cross.groupby('Phylum')['Abundance'].sum().nlargest(12).index.tolist()
top_classes_hm = cross.groupby('Class')['Abundance'].sum().nlargest(20).index.tolist()
cross = cross[cross['Phylum'].isin(top_phyla_hm) & cross['Class'].isin(top_classes_hm)]

if not cross.empty:
    cross['LogAbundance'] = cross['Abundance'].apply(lambda x: math.log10(x + 1))

    source_hm = ColumnDataSource(cross)

    phyla_list = cross.groupby('Phylum')['Abundance'].sum().sort_values(ascending=False).index.tolist()
    class_list = cross.groupby('Class')['Abundance'].sum().sort_values(ascending=False).index.tolist()

    mapper = LinearColorMapper(
        palette=list(reversed(Turbo256)),
        low=cross['LogAbundance'].min(),
        high=cross['LogAbundance'].max())

    p_hm = figure(
        title="Taxonomy Heatmap — Phylum × Class (log₁₀ abundance)",
        x_range=phyla_list,
        y_range=list(reversed(class_list)),
        width=900,
        height=max(400, len(class_list) * 28 + 100),
        tools="hover,save",
        tooltips=[
            ("Phylum", "@Phylum"),
            ("Class", "@Class"),
            ("Reads", "@Abundance{0,0}"),
        ])
    p_hm.rect(x='Phylum', y='Class', width=1, height=1, source=source_hm,
              fill_color=transform('LogAbundance', mapper), line_color=None)
    p_hm.xaxis.major_label_orientation = 0.8

    color_bar = ColorBar(color_mapper=mapper, location=(0, 0),
                         ticker=BasicTicker())
    p_hm.add_layout(color_bar, 'right')

    output_file(f'{OUT_DIR}/taxonomy_heatmap.html',
                title="Taxonomy Heatmap")
    save(p_hm)
else:
    print("  Skipping heatmap — not enough cross-classified data.")

# ==========================================================================
# Krona — Generate tab-separated input file for ktImportText
# ==========================================================================
print("Generating Krona input file...")

krona_rows = []
for _, row in df.iterrows():
    levels = [row[lvl] for lvl in tax_levels if row[lvl] != 'Unassigned']
    krona_rows.append([str(int(row['Abundance']))] + levels)

krona_path = f'{OUT_DIR}/krona_input.txt'
with open(krona_path, 'w') as fh:
    for row in krona_rows:
        fh.write('\t'.join(row) + '\n')

print(f"  Wrote {krona_path} ({len(krona_rows)} ASVs)")

# ==========================================================================
# Summary
# ==========================================================================
total_asvs = len(df)
total_phyla = df['Phylum'].nunique()
dominant_phylum = df.groupby('Phylum')['Abundance'].sum().idxmax()
dominant_pct = df.groupby('Phylum')['Abundance'].sum().max() / total_reads * 100

print(f"\n--- Sample Summary ---")
print(f"  Total ASVs:        {total_asvs}")
print(f"  Total reads:       {int(total_reads):,}")
print(f"  Phyla detected:    {total_phyla}")
print(f"  Dominant phylum:   {dominant_phylum} ({dominant_pct:.1f}%)")
print(f"\nPlotly outputs:")
print(f"  taxonomic_sunburst.html")
print(f"  taxonomic_treemap.html")
print(f"  taxonomic_pie_chart.html")
print(f"Bokeh outputs:")
print(f"  rank_abundance_curve.html")
print(f"  phylum_composition_bar.html")
print(f"  taxonomy_heatmap.html")
print(f"Krona input:")
print(f"  krona_input.txt")
EOF

python "$PLOT_SCRIPT" \
  --sample-name "$SAMPLE_NAME" \
  --data-dir "${OUTPUT_DIR}/exported-data" \
  --output-dir "$OUTPUT_DIR"

# --------------------------------------------------------------------------
# Krona HTML generation (runs outside Python, requires ktImportText)
# --------------------------------------------------------------------------
if [[ "$KRONA_AVAILABLE" == true ]] && command -v ktImportText &>/dev/null; then
    echo "  Generating Krona interactive chart..."
    ktImportText \
        "${OUTPUT_DIR}/krona_input.txt" \
        -o "${OUTPUT_DIR}/krona_taxonomy.html" \
        -n "$SAMPLE_NAME"
    echo "  Created: krona_taxonomy.html"
else
    echo "  Skipped Krona HTML (ktImportText not available)."
    echo "  To generate manually: ktImportText ${OUTPUT_DIR}/krona_input.txt -o krona_taxonomy.html"
fi

echo "========================================="
echo "          PIPELINE COMPLETE!             "
echo "========================================="
echo "Check $OUTPUT_DIR for your final outputs."
