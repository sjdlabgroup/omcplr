import sys
import pandas as pd
import re
from collections import defaultdict
import string

def parse_and_label_hotspots(
    input_report="overlapping_sv_precise_report.tsv", 
    output_summary="simulation_complex_sv_wrapper_report.tsv"
):
    """
    Reads the overlapping SV precise report, maps unique breakpoints per cluster line
    to lowercase alphabetical pairs (e.g., a:b, c:d), and reports:
      1. Precise chromosomal coordinates for each breakpoint pair.
      2. Individual breakpoint frequencies across all cells.
      3. Multi-breakpoint combination frequencies per cell.
    """
    try:
        df = pd.read_csv(input_report, sep='\t')
    except Exception as e:
        print(f"Error reading file {input_report}: {e}")
        return

    def get_alphabet_labels():
        letters = string.ascii_lowercase
        for l in letters:
            yield l
        for l in letters:
            yield l * 2

    processed_rows = []

    for index, row in df.iterrows():
        chrom = str(row['chrom'])
        cluster_start = row['cluster_start']
        cluster_end = row['cluster_end']
        breakdown_str = row['event_type_breakdown']

        if pd.isna(breakdown_str):
            continue

        events = [e.strip() for e in breakdown_str.split('|')]

        # Step 1: Assign alphabetic pair labels to unique coordinate ranges
        breakpoint_to_label = {}
        label_to_coord = {}
        label_generator = get_alphabet_labels()

        # Step 2: Track list of labels per cell
        cell_to_labels = defaultdict(list)

        for event in events:
            match = re.match(r"([^:]+):(\d+):(\d+):([^(]+)\(cells:([^)]+)\)", event)
            if match:
                ev_chrom, start, end, ev_type, cells_str = match.groups()
                coord_key = f"{start}:{end}"

                if coord_key not in breakpoint_to_label:
                    lbl1 = next(label_generator)
                    lbl2 = next(label_generator)
                    assigned_label = f"{lbl1}:{lbl2}"
                    
                    breakpoint_to_label[coord_key] = assigned_label
                    # Store exact chromosome breakpoint mapping: e.g., a:b -> chr1:1234:chr1:2345
                    label_to_coord[assigned_label] = f"{ev_chrom}:{start}:{ev_chrom}:{end}"

                assigned_label = breakpoint_to_label[coord_key]

                cell_ids = [c.strip() for c in cells_str.split(',')]
                for cell in cell_ids:
                    cell_to_labels[cell].append(assigned_label)

        # Step 3: Format mapping summary string (a:b=chr1:1234:chr1:2345 | c:d=chr1:1256:chr1:2351)
        breakpoint_coords_str = " | ".join([
            f"{label}={coord}" for label, coord in sorted(label_to_coord.items(), key=lambda x: x[0])
        ])

        breakpoint_mapping_str = " | ".join([f"{k} -> {v}" for k, v in breakpoint_to_label.items()])
        all_labels_used = ", ".join(breakpoint_to_label.values())

        # Step 4: Compute single counts vs cell-level combination counts
        single_counts = defaultdict(int)
        combination_counts = defaultdict(int)

        for cell, labels in cell_to_labels.items():
            unique_cell_labels = set(labels)
            for lbl in unique_cell_labels:
                single_counts[lbl] += 1
            
            sorted_combo = ", ".join(sorted(unique_cell_labels))
            combination_counts[sorted_combo] += 1

        # Format outputs (omitting zero counts automatically)
        single_counts_str = "; ".join([
            f"{k}={v}" for k, v in sorted(single_counts.items(), key=lambda x: x[0]) if v > 0
        ])

        combination_counts_str = "; ".join([
            f"{k}={v}" for k, v in sorted(combination_counts.items(), key=lambda x: x[0]) if v > 0
        ])

        processed_rows.append({
            'chrom': chrom,
            'cluster_start': cluster_start,
            'cluster_end': cluster_end,
            'breakpoint_coordinates': breakpoint_coords_str,
            'total_unique_breakpoints': len(breakpoint_to_label),
            'breakpoint_labels_set': all_labels_used,
            'breakpoint_mapping_details': breakpoint_mapping_str,
            'single_breakpoint_counts': single_counts_str,
            'combination_breakpoint_counts': combination_counts_str
        })

    output_df = pd.DataFrame(processed_rows)
    output_df.to_csv(output_summary, sep='\t', index=False)

    print(f"\n--- Analysis Complete ---")
    print(f"Successfully processed {len(output_df)} hotspots.")
    print(f"Summary report saved to: '{output_summary}'\n")

    pd.set_option('display.max_colwidth', None)
    print(output_df[['chrom', 'cluster_start', 'cluster_end', 'breakpoint_coordinates', 'single_breakpoint_counts', 'combination_breakpoint_counts']].head(5).to_string(index=False))


if __name__ == "__main__":
    target_input = "overlapping_sv_precise_report.tsv"
    target_output = "simulation_complex_sv_wrapper_report.tsv"

    if len(sys.argv) > 1:
        target_input = sys.argv[1]
    if len(sys.argv) > 2:
        target_output = sys.argv[2]

    parse_and_label_hotspots(
        input_report=target_input, 
        output_summary=target_output
    )