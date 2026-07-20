import os
import sys
import glob
import re
import pandas as pd
from collections import defaultdict

def find_all_overlapping_svs_with_precise_breakdowns(input_dir="cell_beds", output_report="overlapping_sv_precise_report.tsv"):
    """
    Reads all per-cell BED files, groups overlapping/nested events spatially,
    and reports cluster footprints ONLY if they contain at least 2 independently
    occurring or unique structural variant events (ignoring clonal duplicates).
    """
    if not os.path.exists(input_dir):
        print(f"Error: Directory '{input_dir}' does not exist.")
        return

    # Find all cell BED files in the folder
    bed_files = glob.glob(os.path.join(input_dir, "cell_*.bed"))
    if not bed_files:
        print(f"No cell BED files found in '{input_dir}'.")
        return

    print(f"Found {len(bed_files)} cell profile files. Parsing data from '{input_dir}'...")

    # Step 1: Collect unique events per cell first to avoid intra-cell file line duplication
    all_records = []

    for file_path in bed_files:
        match = re.search(r"cell_(\d+)\.bed", os.path.basename(file_path))
        cell_id = int(match.group(1)) if match else os.path.basename(file_path)

        try:
            df = pd.read_csv(
                file_path,
                sep=r'\s+',
                header=None,
                names=['chrom', 'start', 'end', 'event_type'],
                usecols=[0, 1, 2, 3]
            )
            # Drop identical duplicates within the exact same file if any exist
            df = df.drop_duplicates()

            for _, row in df.iterrows():
                all_records.append({
                    'chrom': str(row['chrom']),
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'event_type': str(row['event_type']),
                    'cell_id': cell_id
                })
        except pd.errors.EmptyDataError:
            continue  # Skip empty files safely

    if not all_records:
        print("No structural variant data found inside the files.")
        return

    master_df = pd.DataFrame(all_records)

    # Step 2: Sort records strictly by genomic coordinates
    master_df = master_df.sort_values(by=['chrom', 'start', 'end']).reset_index(drop=True)

    print("Clustering overlapping genomic windows & stripping clonal duplicates...")

    # Step 3: Single-pass interval merging protocol
    clustered_events = []

    def create_new_cluster(row):
        initial_event = {
            'start': row['start'],
            'end': row['end'],
            'event_type': row['event_type'],
            'cell_id': row['cell_id']
        }
        return {
            'chrom': row['chrom'],
            'start': row['start'],
            'end': row['end'],
            'all_cells': {row['cell_id']},
            'raw_events': [initial_event]
        }

    # Initialize the first cluster
    current_cluster = create_new_cluster(master_df.iloc[0])

    for i in range(1, len(master_df)):
        row = master_df.iloc[i]

        # Overlap criteria: same chromosome and overlapping boundaries
        if row['chrom'] == current_cluster['chrom'] and row['start'] <= current_cluster['end']:
            current_cluster['end'] = max(current_cluster['end'], row['end'])
            current_cluster['all_cells'].add(row['cell_id'])
            current_cluster['raw_events'].append({
                'start': row['start'],
                'end': row['end'],
                'event_type': row['event_type'],
                'cell_id': row['cell_id']
            })
        else:
            # Evaluate the completed cluster block before saving
            if evaluate_cluster_validity(current_cluster):
                clustered_events.append(current_cluster)

            current_cluster = create_new_cluster(row)

    # Catch final remaining active cluster block
    if evaluate_cluster_validity(current_cluster):
        clustered_events.append(current_cluster)

    if not clustered_events:
        print("\n--- Analysis Complete ---")
        print("No true multi-event overlapping hotspots found after filtering out clonal lineage duplicates.")
        return

    # Step 4: Map filtered cluster metrics to clean structured output columns
    report_data = []
    for cluster in clustered_events:
        chrom = cluster['chrom']
        breakdown_parts = []

        # Group cell IDs sharing the exact same coordinate mutation to make output compact
        # Maps: (start, end, event_type) -> list of cell_ids
        distinct_mutation_tracks = defaultdict(list)
        for ev in cluster['raw_events']:
            distinct_mutation_tracks[(ev['start'], ev['end'], ev['event_type'])].append(ev['cell_id'])

        # Sort variants by start coordinate for sequential printing
        sorted_mutations = sorted(distinct_mutation_tracks.items(), key=lambda x: (x[0][0], x[0][1]))

        for (start, end, e_type), cells in sorted_mutations:
            cells_str = ",".join(map(str, sorted(cells)))
            # Format: chr:start:end:type(cells_sharing_this_exact_event)
            breakdown_parts.append(f"{chrom}:{start}:{end}:{e_type}(cells:{cells_str})")

        event_breakdown_str = " | ".join(breakdown_parts)

        report_data.append({
            'chrom': chrom,
            'cluster_start': cluster['start'],
            'cluster_end': cluster['end'],
            'total_unique_cells': len(cluster['all_cells']),
            'distinct_events_count': len(distinct_mutation_tracks),
            'event_type_breakdown': event_breakdown_str,
            'all_cell_ids': ",".join(map(str, sorted(cluster['all_cells'])))
        })

    report_df = pd.DataFrame(report_data)

    # Sort by the number of distinct overlapping genomic mutations, then by coordinate
    report_df = report_df.sort_values(by=['distinct_events_count', 'chrom', 'cluster_start'], ascending=[False, True, True])

    # Export report to disk
    report_df.to_csv(output_report, sep='\t', index=False)

    print(f"\n--- Analysis Complete ---")
    print(f"Identified {len(report_df)} true multi-variant hotspots.")
    print(f"Output summary exported to: '{output_report}'")
    print("\nTop 10 absolute genomic hotspot clusters with precise breakdowns:")

    pd.set_option('display.max_colwidth', None)
    print(report_df.drop(columns=['all_cell_ids']).head(10))


def evaluate_cluster_validity(cluster):
    """
    Returns True if the cluster contains at least 2 distinct genomic mutations.
    Collapses identical shared inherited events across cells down into 1 occurrence.
    """
    # Create a unique set of coordinates and event types inside this cluster window
    unique_mutations = set()
    for ev in cluster['raw_events']:
        unique_mutations.add((ev['start'], ev['end'], ev['event_type']))

    # If the unique mutations set size >= 2, we have a true overlapping variant hotspot
    return len(unique_mutations) >= 2


if __name__ == "__main__":
    # Default values fallback
    target_input_dir = "20_sim_cell_svs"
    target_output_report = "overlapping_sv_precise_report.tsv"

    # If parameters are supplied as command line arguments (e.g., python script.py my_input_folder output_name.tsv)
    if len(sys.argv) > 1:
        target_input_dir = sys.argv[1]
    if len(sys.argv) > 2:
        target_output_report = sys.argv[2]

    find_all_overlapping_svs_with_precise_breakdowns(
        input_dir=target_input_dir, 
        output_report=target_output_report
    )