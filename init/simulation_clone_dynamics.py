import os
import random
import sys
import pandas as pd
import numpy as np
from collections import Counter

class Cell:
    def __init__(self, cell_id, parent_id=None, genome_config=None, rates=None):
        self.cell_id = cell_id
        self.parent_id = parent_id
        self.rates = rates.copy() if rates else {}

        # Open DNA breaks active during the current time-step
        if genome_config:
            self.chromosomes = {chrom: [] for chrom in genome_config.keys()}
        else:
            self.chromosomes = {}

        # Inherited or newly acquired structural variants: list of tuples (chrom, start, end, sv_type)
        self.sv_records = []
        self.num_drivers = 0

    def sample_dna_breaks(self, genome_config):
        """Samples DNA breaks across chromosomes based on a single DNA break rate."""
        break_rate = self.rates.get('dna_break_rate', 0.01)

        for chrom, size in genome_config.items():
            num_breaks = np.random.poisson(break_rate * (size / 1e6)) # rate per Mb
            for _ in range(num_breaks):
                coord = random.randint(0, size)
                self.chromosomes[chrom].append(coord)
                offset = int(np.clip(np.random.poisson(5000), 100, 1000000))
                coord = min(coord + offset, size)
                self.chromosomes[chrom].append(coord)

    def process_sv_events(self, driver_prob, genome_config=None):
        """Processes background SVs and localized single step clustered events (e.g. chromothripsis)."""
        one_step_clustered_event_prob = self.rates.get('one_step_clustered_event_probability', 0.0)
        one_step_clustered_event_complexity = self.rates.get('one_step_clustered_event_complexity', 0.0)

        # ----------------------------------------------------
        # 1. Single Step Clustered Event (e.g. Chromothripsis) Generation
        # ----------------------------------------------------
        if random.random() < one_step_clustered_event_prob and genome_config:
            # Randomly select a target chromosome for localized shatter/rearrangement
            target_chrom = random.choice(list(genome_config.keys()))
            chrom_length = genome_config[target_chrom]

            # Number of localized rearrangements drawn from Poisson distribution
            num_rearrangements = np.random.poisson(one_step_clustered_event_complexity)

            if num_rearrangements > 0:
                # Define a localized region (e.g., spanning 1-10% of the chromosome)
                region_size = int(chrom_length * random.uniform(0.01, 0.10))
                region_start = random.randint(0, max(0, chrom_length - region_size))
                region_end = region_start + region_size

                # Generate localized breakpoints within the shattered region
                breakpoints = sorted([
                    random.randint(region_start, region_end)
                    for _ in range(num_rearrangements * 2)
                ])

                # Pair breakpoints to form localized complex rearrangements
                for i in range(0, len(breakpoints) - 1, 2):
                    bp1 = breakpoints[i]
                    bp2 = breakpoints[i + 1]
                    if bp2 > bp1:
                        sv_type = 'one_step_clustered_sv_events'
                        self.sv_records.append((target_chrom, bp1, bp2, sv_type))

                        if random.random() < driver_prob:
                            self.num_drivers += 1
                            self.update_phenotype()

        # ----------------------------------------------------
        # 2. Standard Background SV Events
        # ----------------------------------------------------
        for chrom, breaks in self.chromosomes.items():
            if len(breaks) >= 2:
                breaks.sort()

                while len(breaks) >= 2:
                    bp1 = breaks.pop(0)
                    bp2 = breaks.pop(0)

                    sv_type = 'sv'

                    # Record the genomic alteration coordinates, after rejecting excessively large events
                    if bp2 - bp1 < 1000000:
                        self.sv_records.append((chrom, bp1, bp2, sv_type))

                        # Determine if it's a driver or passenger event
                        is_driver = random.random() < driver_prob
                        if is_driver:
                            self.num_drivers += 1
                            self.update_phenotype()

        # Clear temporary daily open breaks
        for chrom in self.chromosomes:
            self.chromosomes[chrom] = []

    def update_phenotype(self):
        """Adjusts birth or death rate by up to 20% per driver event."""
        modifier = 1.0 + random.uniform(0.0, 0.20)
        self.rates['birth_rate'] = max(0.0, self.rates['birth_rate'] * modifier)

    def clone(self, next_id):
        """Creates a daughter cell inheriting the structural variants and adjusted rates."""
        daughter = Cell(cell_id=next_id, parent_id=self.cell_id, rates=self.rates)
        daughter.num_drivers = self.num_drivers
        # Pass down clonal genomic alterations explicitly
        daughter.sv_records = list(self.sv_records)
        daughter.chromosomes = {chrom: [] for chrom in self.chromosomes.keys()}
        return daughter


def load_inputs(bed_path, rates_path):
    # Robust reading using whitespace regex to catch spaces or tabs
    bed_df = pd.read_csv(bed_path, sep=r'\s+', header=None, names=['chrom', 'start', 'end'], usecols=[0, 1, 2], comment='#')
    bed_df = bed_df.dropna(subset=['chrom', 'start', 'end'])

    if bed_df.empty:
        print(f"Error: No valid genomic tracks found in {bed_path}. Check delimiters.")
        sys.exit(1)

    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)
    genome_config = {row['chrom']: int(row['end'] - row['start']) for _, row in bed_df.iterrows()}

    rates_df = pd.read_csv(rates_path, sep=r'\s+')
    rates_df = rates_df.dropna(subset=['parameter', 'value'])

    rates = {}
    for _, row in rates_df.iterrows():
        param = row['parameter']
        val = row['value']
        try:
            rates[param] = float(val)
        except ValueError:
            rates[param] = str(val)

    # Extract engine execution variables dynamically from parameters.txt with fallback defaults
    engine_config = {
        'output_dir': str(rates.get('output_dir', 'sim_cell_svs')),
        'max_time_steps': int(float(rates.get('max_time_steps', 5))),
        'max_cells': int(float(rates.get('max_cells', 10000))),
    }

    # Robust boolean parsing for report_all_timepoints parameter
    raw_report = rates.get('report_all_timepoints', 'False')
    if isinstance(raw_report, str):
        engine_config['report_all_timepoints'] = raw_report.strip().lower() in ['true', '1', 'yes']
    else:
        engine_config['report_all_timepoints'] = bool(raw_report)

    return genome_config, rates, engine_config


def export_alterations_to_bed(population, output_bed_path):
    """Aggregates matching SVs across the living population to output a custom bulk BED file."""
    all_svs = []
    for cell in population:
        all_svs.extend(cell.sv_records)

    sv_counts = Counter(all_svs)

    bed_rows = []
    for (chrom, start, end, sv_type), count in sv_counts.items():
        bed_rows.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'event_type': sv_type,
            'cell_count': count
        })

    if bed_rows:
        output_df = pd.DataFrame(bed_rows)
        output_df = output_df.sort_values(by=['chrom', 'start', 'end'])
        output_df.to_csv(output_bed_path, sep='\t', index=False, header=False)
        print(f"Successfully exported bulk summary to: {output_bed_path}")
    else:
        print("No structural variants occurred in the surviving population to export.")


def export_per_cell_beds(population, output_dir, t=None):
    """Creates an output folder and exports individual BED files for each living cell."""
    os.makedirs(output_dir, exist_ok=True)

    cells_written = 0
    for cell in population:
        if not cell.sv_records:
            continue

        cell_rows = []
        for (chrom, start, end, sv_type) in cell.sv_records:
            cell_rows.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'event_type': sv_type
            })

        cell_df = pd.DataFrame(cell_rows)
        cell_df = cell_df.sort_values(by=['chrom', 'start', 'end'])

        cell_file_path = os.path.join(output_dir, f"cell_{cell.cell_id}.bed")
        cell_df.to_csv(cell_file_path, sep='\t', index=False, header=False)
        cells_written += 1

    if t is not None:
        print(f"Time {t:02d}: Exported per-cell events to the folder: {output_dir}.")
    else:
        print(f"Final Step: Exported per-cell events to the folder: {output_dir}.")


def run_simulation(bed_path, rates_path, output_bed_path="simulated_alterations.bed"):
    genome_config, rates, engine_config = load_inputs(bed_path, rates_path)

    output_dir = engine_config['output_dir']
    max_time_steps = engine_config['max_time_steps']
    max_cells = engine_config['max_cells']
    report_all_timepoints = engine_config['report_all_timepoints']

    cell_counter = 1
    root_cell = Cell(cell_id=cell_counter, genome_config=genome_config, rates=rates)
    population = [root_cell]

    driver_prob = rates.get('driver_probability', 0.01)

    print(f"Starting simulation with 1 cell across {len(genome_config)} chromosomes...")
    print(f"Config details -> Max Steps: {max_time_steps} | Max Cells: {max_cells} | Report All: {report_all_timepoints}")

    for t in range(1, max_time_steps + 1):
        if not population or len(population) > max_cells:
            break

        next_generation = []

        for cell in population:
            # 1. Death cycle
            if random.random() < cell.rates.get('death_rate', 0.1):
                continue

            next_generation.append(cell)

            # 2. Birth cycle
            if random.random() < cell.rates.get('birth_rate', 0.2):
                cell_counter += 1
                daughter = cell.clone(cell_counter)

                # 3. Apply structural variants upon creation
                daughter.sample_dna_breaks(genome_config)
                daughter.process_sv_events(driver_prob, genome_config=genome_config)

                next_generation.append(daughter)

        population = next_generation

        if report_all_timepoints:
            timepoint_dir = f"{output_dir}_{t}"
            export_per_cell_beds(population, timepoint_dir, t)
        else:
            print(f"Time {t:02d}: Population size = {len(population)}")

    export_alterations_to_bed(population, output_bed_path)

    if not report_all_timepoints:
        export_per_cell_beds(population, output_dir)

    return population


if __name__ == "__main__":
    run_simulation(
        bed_path="chromosomes.bed",
        rates_path="parameters.txt",
        output_bed_path="simulated_alterations.bed"
    )
