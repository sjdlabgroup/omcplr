import pandas as pd
import re
from collections import defaultdict
import string

def parse_and_label_hotspots(input_report="simulation_complex_events_report.tsv", output_summary="simulation_complex_sv_wrapper_report.tsv"):
    """
    Reads the overlapping SV precise report, maps unique breakpoints per cluster line
    to lowercase alphabetical pairs (e.g., a:b, c:d), lists all available breakpoints
    per line, and counts the combination profile of these labeled events across individual cells.
    """
    # Load the TSV file
    try:
        df = pd.read_csv(input_report, sep='\t')
    except Exception as e:
        print(f"Error reading file {input_report}: {e}")
        return

    # To generate labels: a, b, c, ..., z, aa, bb, cc... if we run out of letters
    def get_alphabet_labels():
        letters = string.ascii_lowercase
        for l in letters:
            yield l
        # Fallback if a cluster has more than 26 distinct breakpoints
        for l in letters:
            yield l * 2

    processed_rows = []

    # Loop through each genomic cluster row
    for index, row in df.iterrows():
        chrom = row['chrom']
        cluster_start = row['cluster_start']
        cluster_end = row['cluster_end']
        breakdown_str = row['event_type_breakdown']

        if pd.isna(breakdown_str):
            continue

        # Split the individual event parts by the pipe separator
        events = [e.strip() for e in breakdown_str.split('|')]

        # Step 1: Map each unique breakpoint (start:end) in this line to a character pair label
        breakpoint_to_label = {}
        label_generator = get_alphabet_labels()

        # Step 2: Track which cells have which event labels
        cell_to_labels = defaultdict(list)

        for event in events:
            # RegEx to capture: chrom:start:end:type(cells:1,2,3)
            match = re.match(r"([^:]+):(\d+):(\d+):([^(]+)\(cells:([^)]+)\)", event)
            if match:
                ev_chrom, start, end, ev_type, cells_str = match.groups()
                coord_key = f"{start}:{end}"

                # Assign labels dynamically (e.g., 'a:b' for the first unique coordinate pair)
                if coord_key not in breakpoint_to_label:
                    lbl1 = next(label_generator)
                    lbl2 = next(label_generator)
                    breakpoint_to_label[coord_key] = f"{lbl1}:{lbl2}"

                assigned_label = breakpoint_to_label[coord_key]

                # Parse out the cell IDs sharing this mutation
                cell_ids = [c.strip() for c in cells_str.split(',')]
                for cell in cell_ids:
                    cell_to_labels[cell].append(assigned_label)

        # Step 3: Compile the total sets of breakpoints mapped to their characters
        # e.g., "10012518:10019440 -> a:b, 10014124:10014534 -> c:d"
        breakpoint_mapping_str = " | ".join([f"{k} -> {v}" for k, v in breakpoint_to_label.items()])
        all_labels_used = ", ".join(breakpoint_to_label.values())

        # Step 4: Count the cell frequency for each unique combination of breakpoint labels
        combination_counts = defaultdict(int)
        for cell, labels in cell_to_labels.items():
            # Sort labels so 'a:b, c:d' is counted the same as 'c:d, a:b'
            sorted_combination = ",".join(sorted(labels))
            combination_counts[sorted_combination] += 1

        # Format combinations cleanly as: "a:b=35; a:b,c:d=1"
        combination_profile_str = "; ".join([f"{k}={v}" for k, v in sorted(combination_counts.items())])

        processed_rows.append({
            'chrom': chrom,
            'cluster_start': cluster_start,
            'cluster_end': cluster_end,
            'total_unique_breakpoints': len(breakpoint_to_label),
            'breakpoint_labels_set': all_labels_used,
            'breakpoint_mapping_details': breakpoint_mapping_str,
            'cell_combinations_frequency': combination_profile_str
        })

    # Convert results into a clean dataframe and export
    output_df = pd.DataFrame(processed_rows)
    output_df.to_csv(output_summary, sep='\t', index=False)

    print(f"\n--- Analysis Complete ---")
    print(f"Successfully processed {len(output_df)} hotspots.")
    print(f"Combinations summary report saved to: '{output_summary}'\n")

    # Preview top rows for immediate inspection
    pd.set_option('display.max_colwidth', None)
    print(output_df[['chrom', 'cluster_start', 'cluster_end', 'breakpoint_labels_set', 'cell_combinations_frequency']].head(5).to_string(index=False))

if __name__ == "__main__":
    parse_and_label_hotspots()
