import os
from tqdm import tqdm
from main.ztfdatafiltering import *

# Base paths for your directories and files
base_path = "data/"  # Common base directory

path_raw = os.path.join(base_path, "ztffp_raw_data/")  # Where the raw data is located
path_processed = os.path.join(base_path, "ztffp_processed_data/")  # Processed files location
path_final = os.path.join(base_path, "ztffp_final_data/")  # Final processed and filtered DataFrames
path_plots = os.path.join(base_path, "plots/")  # Directory for saving plots
coordinate_path = os.path.join("/Users/nilsongoncalves/Mestrado/env/ZTF-data-filtering/src/coordinates.txt")  # File with RA/Dec coordinates

results_table_path = os.path.join(os.path.dirname(__file__), "variability_results.csv")
results = []

def main():
    # Step 1: Read RA and Dec from the coordinates file
    ra_list, dec_list = read_coordinates_from_file(coordinate_path)
    
    # Step 2: Process each source
    for ra, dec in tqdm(zip(ra_list, dec_list), total=len(ra_list), desc="Processing"):
        source_name = create_source_name(ra, dec)

        # Step 3: Create the ZTF table (uncomment this line if needed)
        create_ztf_table(path=path_raw, output_directory=path_processed, ra=ra, dec=dec, source_name=source_name, ext_cor=True)

        # Step 4: Plot light curves for each source and save to `path_plots`
        #plot_lc(output_directory=path_processed, path_final=path_final, path_plots=path_plots, source_name=source_name)

        # Step 5: Calculate F_var and store variability results in lists
        #result = plot_hist(output_directory=path_processed, source_name=source_name, ra=ra, dec=dec, results_table_path=results_table_path)
        #results.append(result)

output_directory = os.path.dirname(__file__)

# Generate histograms and statistics
plot_fvar_statistics(results_table_path, output_directory)


if __name__ == "__main__":
    main()
