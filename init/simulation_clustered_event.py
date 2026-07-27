import os
import sys
import glob
import re
import random
import shutil
import pandas as pd
from collections import defaultdict, Counter

def prepare_noisy_directory(src_dir, dest_dir="tsim_cplxev", noise_factor=0.0):
    """
    Copies files from `src_dir` to `dest_dir`. Replaces a given fraction
    defibned by (noise_factor) of files with contents from other randomly selected files.
    """
    if not os.path.exists(src_dir):
        print(f"Error: Source directory '{src_dir}' does not exist.")
        sys.exit(1)

    # Ensure clean destination directory
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    os.makedirs(dest_dir)

    all_files = sorted(glob.glob(os.path.join(src_dir, "*")))
    if not all_files:
        print(f"Error: No files found in source directory '{src_dir}'.")
        sys.exit(1)

    total_files = len(all_files)
    num_noisy_files = int(total_files * noise_factor)

    # Select random files to be substituted
    noisy_indices = set(random.sample(range(total_files), num_noisy_files))

    print(f"Copying {total_files} files from '{src_dir}' to '{dest_dir}'...")
    if num_noisy_files > 0:
        print(f"Applying noise factor {noise_factor}: {num_noisy_files} file(s) will be substituted.")

    for idx, src_file_path in enumerate(all_files):
        filename = os.path.basename(src_file_path)
        dest_file_path = os.path.join(dest_dir, filename)

        if idx in noisy_indices and total_files > 1:
            # Pick a donor file distinct from the current file
            donor_candidates = [f for f in all_files if f != src_file_path]
            donor_file_path = random.choice(donor_candidates)
            shutil.copy2(donor_file_path, dest_file_path)
        else:
            shutil.copy2(src_file_path, dest_file_path)

    return dest_dir


def find_all_overlapping_svs_with_precise_breakdowns(
    input_dir="tsim_cplxev",
    output_report="overlapping_sv_precise_report.tsv",
    min_cells_with_multiple_events=1
):
    """
    Reads all per-cell BED files, groups overlapping/nested events spatially,
    and reports cluster footprints only if at least `min_cells_with_multiple_events`
    contain multiple events within the cluster.
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
        print("No structural variant data found inside the files. Exciting!")
        return

    master_df = pd.DataFrame(all_records)

    # Step 2: Sort records strictly by genomic coordinates
    master_df = master_df.sort_values(by=['chrom', 'start', 'end']).reset_index(drop=True)

    print(f"Clustering overlapping genomic windows (requiring >={min_cells_with_multiple_events} cell(s) with multiple events)...")

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
            if evaluate_cluster_validity(current_cluster, min_cells_with_multiple_events):
                clustered_events.append(current_cluster)

            current_cluster = create_new_cluster(row)

    # Catch final remaining active cluster block
    if evaluate_cluster_validity(current_cluster, min_cells_with_multiple_events):
        clustered_events.append(current_cluster)

    if not clustered_events:
        print("\n--- Analysis Complete ---")
        print("No hotspot clusters met the criteria for multi-event cell density. Repeat the simulation!")
        return

    # Step 4: Map filtered cluster metrics to clean structured output columns
    report_data = []
    for cluster in clustered_events:
        chrom = cluster['chrom']
        breakdown_parts = []

        distinct_mutation_tracks = defaultdict(list)
        for ev in cluster['raw_events']:
            distinct_mutation_tracks[(ev['start'], ev['end'], ev['event_type'])].append(ev['cell_id'])

        sorted_mutations = sorted(distinct_mutation_tracks.items(), key=lambda x: (x[0][0], x[0][1]))

        for (start, end, e_type), cells in sorted_mutations:
            cells_str = ",".join(map(str, sorted(cells)))
            breakdown_parts.append(f"{chrom}:{start}:{end}:{e_type}(cells:{cells_str})")

        event_breakdown_str = " | ".join(breakdown_parts)

        # Count cells having multiple events in this specific cluster
        cell_event_counts = Counter(ev['cell_id'] for ev in cluster['raw_events'])
        cells_with_multi = sum(1 for count in cell_event_counts.values() if count >= 2)

        report_data.append({
            'chrom': chrom,
            'cluster_start': cluster['start'],
            'cluster_end': cluster['end'],
            'total_unique_cells': len(cluster['all_cells']),
            'cells_with_multi_events': cells_with_multi,
            'distinct_events_count': len(distinct_mutation_tracks),
            'event_type_breakdown': event_breakdown_str,
            'all_cell_ids': ",".join(map(str, sorted(cluster['all_cells'])))
        })

    report_df = pd.DataFrame(report_data)

    report_df = report_df.sort_values(
        by=['cells_with_multi_events', 'distinct_events_count', 'chrom', 'cluster_start'],
        ascending=[False, False, True, True]
    )

    report_df.to_csv(output_report, sep='\t', index=False)

    print(f"\n--- Analysis Complete ---")
    print(f"Identified {len(report_df)} multi-variant hotspots matching criteria.")
    print(f"Output summary exported to: '{output_report}'")
    print("\nTop 10 hotspot clusters:")

    pd.set_option('display.max_colwidth', None)
    print(report_df.drop(columns=['all_cell_ids']).head(10))


def evaluate_cluster_validity(cluster, min_cells_threshold):
    """
    Returns True if at least `min_cells_threshold` cells contain
    2 or more raw events within this cluster window.
    """
    cell_counts = Counter(ev['cell_id'] for ev in cluster['raw_events'])
    cells_with_multiple_events = sum(1 for count in cell_counts.values() if count >= 2)

    return cells_with_multiple_events >= min_cells_threshold


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])

    # Validate positional arguments (requires 4 arguments)
    if len(sys.argv) < 5:
        print(f"Usage: python {script_name} [input_folder] [output_report.txt] [min_cells] [noise_factor]")
        print("Example: python simulation_clustered_event.py 10_sim_run output_report.tsv 1 0.2")
        sys.exit(1)

    source_input_dir = sys.argv[1]
    target_output_report = sys.argv[2]
    min_cells = int(sys.argv[3])

    try:
        noise_factor = float(sys.argv[4])
        if not (0.0 <= noise_factor <= 1.0):
            raise ValueError
    except ValueError:
        print("Error: noise_factor must be a fraction between 0.0 and 1.0. Keep it 0 or a small value!")
        sys.exit(1)

    # 1. Copy files and apply noise substitutions to te destination folder
    working_dir = prepare_noisy_directory(
        src_dir=source_input_dir,
        dest_dir="tsim_cplxev",
        noise_factor=noise_factor
    )

    # 2. Run analysis on the newly populated 'tsim_cplxev' directory
    find_all_overlapping_svs_with_precise_breakdowns(
        input_dir=working_dir,
        output_report=target_output_report,
        min_cells_with_multiple_events=min_cells
    )
