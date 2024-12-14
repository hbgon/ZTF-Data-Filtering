import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import FK5, SkyCoord
from astroquery.ipac.irsa.irsa_dust import IrsaDust

# Convert magnitude to flux
def mag_to_flux(ab_mag, wl):
    fnu = (10.0 ** (-0.4 * (48.6 + ab_mag)))
    return (2.99792458e+18 * fnu) / (wl ** 2.0)

# Convert magnitude error to flux error
def dmag_to_df(dmag, flux):
    return flux * 0.4 * np.log(10) * dmag

# Get E(B-V) Galactic extinction
def get_ebv(ra, dec):
    coo = SkyCoord(ra=float(ra), dec=float(dec), unit="degree", frame=FK5)
    table = IrsaDust.get_query_table(coo, section="ebv")
    return table["ext SandF mean"][0]

def create_sdss_name(ra_dec_str):
    """
    Converts a source name (RA and Dec) into an SDSS-formatted name.

    Parameters:
    - ra_dec_str: Source name string in the format 'RA_Dec' (e.g., '150.0_2.2').

    Returns:
    - name_sdss: Formatted SDSS name as a string.
    """
    try:
        # Split coordinates
        ra_str, dec_str = ra_dec_str.split('_')
        
        # Convert coordinates from string to float
        ra, dec = float(ra_str), float(dec_str)

        # Format the name: SDSS JHHMMSS.SS+DDMMSS.S
        ra_hms = coords.Angle(ra, unit=u.deg).to_string(unit=u.hour, sep='', precision=3, pad=True)
        dec_dms = coords.Angle(dec, unit=u.deg).to_string(unit=u.deg, sep='', precision=2, alwayssign=True, pad=True)
        name_sdss = f"SDSS J{ra_hms}{dec_dms}"

        return name_sdss
    except ValueError:
        print(f"Error: Invalid input format. Expected 'RA_Dec', got '{ra_dec_str}'.")
        return None

def read_coordinates_from_file(path):
    """
    Reads RA and Dec coordinates from a file with space-separated values.

    Parameters:
    - path: Path to the file containing coordinates.

    Returns:
    - ra_list: List of right ascension (RA) values in degrees.
    - dec_list: List of declination (Dec) values in degrees.
    """
    ra_list, dec_list = [], []
    with open(path, 'r') as file:
        for line in file:
            try:
                ra, dec = map(float, line.strip().split())
                ra_list.append(ra)
                dec_list.append(dec)
            except ValueError:
                print(f"Skipping invalid line: {line.strip()}")
    return ra_list, dec_list

def create_source_name(ra, dec):
    """
    Creates a source name from RA and Dec.

    Parameters:
    - ra: Right ascension (float).
    - dec: Declination (float).
    Returns:
    - source_name: Formatted source name as a string.
    """
    ra_str = f"{ra}"
    dec_str = f"{dec}"
    return f"{ra_str}_{dec_str}"

def create_ztf_table(path, output_directory, source_name,ra,dec, ext_cor=True):
    if os.path.exists(os.path.join(path, source_name + '.txt')):
        #print('Doing forced photometry on ZTF data...')
        
        # Correct for galactic extinction or not?
        if not ext_cor:
            ext_corr = np.array([3.74922728, 2.64227246, 1.92969698]) * 0.0
        if ext_cor:
            ext_corr = np.array([3.74922728, 2.64227246, 1.92969698]) * get_ebv(ra,dec)
        
        
        # Loading the giant table. It has to be saved with the name ztf_"source_name".txt
        index, field, ccdid, qid, filter, pid, infobitssci, sciinpseeing, scibckgnd, scisigpix, zpmaginpsci, zpmaginpsciunc, zpmaginpscirms, clrcoeff, clrcoeffunc, ncalmatches, exptime, adpctdif1, adpctdif2, diffmaglim, zpdiff, programid, jd, rfid, forcediffimflux, forcediffimfluxunc, forcediffimsnr, forcediffimchisq, forcediffimfluxap, forcediffimfluxuncap, forcediffimsnrap, aperturecorr, dnearestrefsrc, nearestrefmag, nearestrefmagunc, nearestrefchi, nearestrefsharp, refjdstart, refjdend, procstatus \
            = np.genfromtxt(os.path.join(path, source_name + '.txt'), delimiter=" ", comments='#', dtype=str, skip_header=57,
                            unpack=True)

        # why/how I'm doing this is described in the tutorial
        nearestrefmag[nearestrefmag == 'null'] = 0
        nearestrefmagunc[nearestrefmagunc == 'null'] = 0
        forcediffimflux[forcediffimflux == 'null'] = 0
        forcediffimfluxunc[forcediffimfluxunc == 'null'] = 0

        nearestrefmag = np.array(nearestrefmag, dtype=float)
        nearestrefmagunc = np.array(nearestrefmagunc, dtype=float)
        zpdiff = np.array(zpdiff, dtype=float)
        forcediffimfluxunc = np.array(forcediffimfluxunc, dtype=float)
        forcediffimflux = np.array(forcediffimflux, dtype=float)
        mjd = np.array(jd, dtype=float) - 2400000.5 # changing from Julian Date to Modified Julian Date

        nearestrefflux = 10 ** (0.4 * (zpdiff - nearestrefmag))
        nearestrefflux[nearestrefmag == 0] = 0
        nearestreffluxunc = nearestrefmagunc * nearestrefflux / 1.0857 
        nearestreffluxunc[nearestrefmag == 0] = 0

        Flux_tot = forcediffimflux + nearestrefflux  # Adding the refernce flux with the diferene photoemtry flux
        Fluxunc_tot = np.sqrt(abs(forcediffimfluxunc ** 2 - nearestreffluxunc ** 2)) # Uncertanty on the flux
        SNR_tot = Flux_tot / Fluxunc_tot    
    
        mag_obs = np.zeros(np.shape(SNR_tot)) 
        mag_err_obs = np.zeros(np.shape(SNR_tot))


        # Calculating photometry in magnitudes of the AB system for both g and r filters
        for i in range(len(SNR_tot)):
            if SNR_tot[i] > 5:
                if filter[i] == 'ZTF_g':
                    mag_obs[i] = zpdiff[i] - 2.5 * np.log10(Flux_tot[i]) - ext_corr[0]
                    mag_err_obs[i] = 1.0857 / SNR_tot[i]
                if filter[i] == 'ZTF_r':
                    mag_obs[i] = zpdiff[i] - 2.5 * np.log10(Flux_tot[i]) - ext_corr[1]
                    mag_err_obs[i] = 1.0857 / SNR_tot[i]
                if filter[i] == 'ZTF_i':
                    mag_obs[i] = zpdiff[i] - 2.5 * np.log10(Flux_tot[i]) - ext_corr[2]
                    mag_err_obs[i] = 1.0857 / SNR_tot[i]
            else:
                if filter[i] == 'ZTF_g':
                    mag_obs[i] = zpdiff[i] - 2.5 * np.log10(3 * Fluxunc_tot[i]) - ext_corr[0]
                    mag_err_obs[i] = np.nan
                if filter[i] == 'ZTF_r':
                    mag_obs[i] = zpdiff[i] - 2.5 * np.log10(3 * Fluxunc_tot[i]) - ext_corr[1]
                    mag_err_obs[i] = np.nan
                if filter[i] == 'ZTF_i':
                    mag_obs[i] = zpdiff[i] - 2.5 * np.log10(3 * Fluxunc_tot[i]) - ext_corr[2]
                    mag_err_obs[i] = np.nan
        
        out_directory = os.path.join(output_directory, f'{ra}_{dec}')
        os.makedirs(out_directory, exist_ok=True)  # Create the folder if it does not exist
        
        # Creating a table with photometry in the g-band
        g_cent_wl =  4722.74
        is_g = filter == 'ZTF_g'
        ztf_g = open(os.path.join(out_directory, f'{source_name}_ztf_g.txt'), 'w')
        #ztf_g.write('# if ab_mag_err = nan, the measurement is an upper limit\n')
        ztf_g.write('mjd' + '\t' + 'ab_mag' + '\t' + 'ab_mag_err' + '\t' + 'flux_dens' + '\t' + 'flux_dens_err' + '\n')
        for yy in range(len(mjd[is_g])):
            ztf_g.write('{:.2f}'.format(mjd[is_g][yy]) + '\t' + '{:.2f}'.format(mag_obs[is_g][yy]) + '\t' + '{:.2f}'.format(
                mag_err_obs[is_g][yy]) + '\t' + '{:.2e}'.format(mag_to_flux(mag_obs[is_g][yy], g_cent_wl)) + '\t' + '{:.2e}'.format(dmag_to_df(mag_err_obs[is_g][yy], mag_to_flux(mag_obs[is_g][yy], g_cent_wl))) + '\n')
        ztf_g.close()

        # Creating a table with photometry in the r-band
        r_cent_wl =  6339.61
        is_r = filter == 'ZTF_r'
        ztf_r = open(os.path.join(out_directory, f'{source_name}_ztf_r.txt'), 'w')
        #ztf_r.write('# if ab_mag_err = nan, the measurement is an upper limit\n')
        ztf_r.write('mjd' + '\t' + 'ab_mag' + '\t' + 'ab_mag_err' + '\t' + 'flux_dens' + '\t' + 'flux_dens_err' + '\n')
        for yy in range(len(mjd[is_r])):
            ztf_r.write('{:.2f}'.format(mjd[is_r][yy]) + '\t' + '{:.2f}'.format(mag_obs[is_r][yy]) + '\t' + '{:.2f}'.format(
                mag_err_obs[is_r][yy]) + '\t' + '{:.2e}'.format(
                mag_to_flux(mag_obs[is_r][yy], r_cent_wl)) + '\t' + '{:.2e}'.format(
                dmag_to_df(mag_err_obs[is_r][yy], mag_to_flux(mag_obs[is_r][yy], r_cent_wl))) + '\n')
        ztf_r.close()

        # Creating a table with photometry in the i-band
        i_cent_wl =  7829.03
        is_i = filter == 'ZTF_i'
        ztf_i = open(os.path.join(out_directory, f'{source_name}_ztf_i.txt'), 'w')
        #ztf_i.write('# if ab_mag_err = nan, the measurement is an upper limit\n')
        ztf_i.write('mjd' + '\t' + 'ab_mag' + '\t' + 'ab_mag_err' + '\t' + 'flux_dens' + '\t' + 'flux_dens_err' + '\n')
        for yy in range(len(mjd[is_i])):
            ztf_i.write('{:.2f}'.format(mjd[is_i][yy]) + '\t' + '{:.2f}'.format(mag_obs[is_i][yy]) + '\t' + '{:.2f}'.format(
                mag_err_obs[is_i][yy]) + '\t' + '{:.2e}'.format(
                mag_to_flux(mag_obs[is_i][yy], i_cent_wl)) + '\t' + '{:.2e}'.format(
                dmag_to_df(mag_err_obs[is_i][yy], mag_to_flux(mag_obs[is_i][yy], i_cent_wl))) + '\n')
        ztf_i.close()

def plot_lc(output_directory, source_name, path_final, path_plots, offset=1.15, figsize=(10, 13), min_observation_points=100):
    """
    Plots light curves for different bands applying filters and saves filtered DataFrames and plots.
    """
    def calculate_f_var(df):
        """Calculate F_var and its error."""
        mean_flux = df['flux'].mean()
        S_squared = np.var(df['flux'], ddof=1)
        mean_fluxe_squared = np.mean(df['fluxe'] ** 2)
        sigma_rms_squared = S_squared - mean_fluxe_squared
        f_var = np.sqrt(sigma_rms_squared) / mean_flux if sigma_rms_squared > 0 else 0
        N = len(df)
        error_f_var = np.sqrt((1 / (2 * N)) * (mean_fluxe_squared / (mean_flux ** 2 * f_var ** 2)) + (mean_fluxe_squared / mean_flux ** 2))
        return f_var, error_f_var

    def filter_rmsd(df):
        """Apply rolling mean and RMSD filtering."""
        df['flux_rolling_mean'] = df['flux'].rolling(window=30, center=True).mean()
        df['flux_rolling_rmsd'] = np.sqrt((df['flux'] - df['flux_rolling_mean']) ** 2).rolling(window=30, center=True).mean()
        return df[np.abs(df['flux'] - df['flux_rolling_mean']) <= 4 * df['flux_rolling_rmsd']]

    # Initialize
    band_list = ['ztf_g', 'ztf_r', 'ztf_i']
    color_list = ['r', 'g', 'darkorange']
    mjd_combined = pd.Series(dtype=float)

    # Create output directories if not existing
    os.makedirs(path_plots, exist_ok=True)
    source_directory = os.path.join(path_final, source_name)
    os.makedirs(source_directory, exist_ok=True)

    plt.close('all')
    fig, axes = plt.subplots(len(band_list), 1, figsize=figsize, sharex=True)

    for i, band in enumerate(band_list):
        file_path = os.path.join(output_directory, source_name, f'{source_name}_{band}.txt')
        if not os.path.exists(file_path):
            print(f"File for {band} not found, skipping...")
            continue

        try:
            # Load data
            mjd, abmag, abmage, flux, fluxe = np.loadtxt(file_path, skiprows=1, unpack=True)
            if np.isnan(flux).all() or np.isnan(fluxe).all():
                print(f"All flux values NaN for {band}, skipping...")
                continue

            # Create and filter DataFrame
            df_band = pd.DataFrame({'mjd': mjd, 'flux': flux, 'fluxe': fluxe}).dropna()
            df_band['flux'] = df_band.groupby('mjd')['flux'].transform('mean')
            df_band['fluxe'] = df_band.groupby('mjd')['fluxe'].transform('mean')
            df_band = df_band.drop_duplicates('mjd')
            df_band = filter_rmsd(df_band)

            if len(df_band) < min_observation_points:
                print(f"Insufficient observations for {band}, skipping...")
                continue

            # Calculate F_var
            f_var, error_f_var = calculate_f_var(df_band)
            if f_var <= 0:
                print(f"F_var <= 0 for {band}, skipping...")
                continue

            # Save filtered DataFrame
            filtered_file_path = os.path.join(source_directory, f'{source_name}_{band}_filtered.txt')
            df_band.to_csv(filtered_file_path, index=False)

            # Plot the light curve
            ax = axes[i]
            ax.errorbar(df_band['mjd'], df_band['flux'], yerr=df_band['fluxe'], color=color_list[i], fmt='o',
                        ecolor='darkgrey', markeredgecolor='darkgrey', markersize=7, label=band)
            ax.set_ylabel('Flux', fontsize=15)
            ax.legend(title=f"F_var: {f_var:.4f}\nError F_var: {error_f_var:.4f}\nn: {len(df_band)}")
            ax.grid(True)

            # Update x-axis limits for combined plot
            mjd_combined = pd.concat([mjd_combined, df_band['mjd']])
        except Exception as e:
            print(f"Error processing {band}: {e}")
            continue

    plt.xlabel('MJD', fontsize=15)
    plt.xlim(mjd_combined.min() - 100, mjd_combined.max() + 100)
    plt.tight_layout()
    plt.savefig(os.path.join(path_plots, f'{source_name}.png'))
    plt.close()

def plot_hist(output_directory, source_name, ra, dec, results_table_path, min_observation_points=100):
    """
    Computes F_var for each band, saves the results in a table, and plots histograms for all bands.

    Parameters:
    - output_directory: Directory containing the ZTF band files.
    - source_name: Name of the source.
    - ra: Right Ascension of the source.
    - dec: Declination of the source.
    - results_table_path: Path to save the consolidated results table.
    - min_observation_points: Minimum number of observations to calculate F_var.

    Returns:
    - result_row: A list containing [ra, dec, f_var_g, f_var_r, f_var_i, error_g, error_r, error_i].
    """
    def calculate_f_var(df):
        """Calculate F_var and its error."""
        mean_flux = df['flux'].mean()
        S_squared = np.var(df['flux'], ddof=1)
        mean_fluxe_squared = np.mean(df['fluxe'] ** 2)
        sigma_rms_squared = S_squared - mean_fluxe_squared
        f_var = np.sqrt(sigma_rms_squared) / mean_flux if sigma_rms_squared > 0 else 0
        N = len(df)
        error_f_var = np.sqrt((1 / (2 * N)) * (mean_fluxe_squared / (mean_flux ** 2 * f_var ** 2)) + (mean_fluxe_squared / mean_flux ** 2))
        return f_var, error_f_var

    def filter_rmsd(df):
        """Apply rolling mean and RMSD filtering."""
        df['flux_rolling_mean'] = df['flux'].rolling(window=30, center=True).mean()
        df['flux_rolling_rmsd'] = np.sqrt((df['flux'] - df['flux_rolling_mean']) ** 2).rolling(window=30, center=True).mean()
        return df[np.abs(df['flux'] - df['flux_rolling_mean']) <= 4 * df['flux_rolling_rmsd']]

    # Initialize results
    f_var_g, f_var_r, f_var_i = None, None, None
    error_g, error_r, error_i = None, None, None
    f_var_hist = []  # Collect F_var values for histogram

    # Process each band
    band_list = ['ztf_g', 'ztf_r', 'ztf_i']
    for band in band_list:
        file_path = os.path.join(output_directory, source_name, f'{source_name}_{band}.txt')
        if not os.path.exists(file_path):
            print(f"File for {band} not found, skipping...")
            continue

        try:
            # Load data
            mjd, abmag, abmage, flux, fluxe = np.loadtxt(file_path, skiprows=1, unpack=True)
            if np.isnan(flux).all() or np.isnan(fluxe).all():
                print(f"All flux values NaN for {band}, skipping...")
                continue

            # Create and filter DataFrame
            df_band = pd.DataFrame({'mjd': mjd, 'flux': flux, 'fluxe': fluxe}).dropna()
            df_band['flux'] = df_band.groupby('mjd')['flux'].transform('mean')
            df_band['fluxe'] = df_band.groupby('mjd')['fluxe'].transform('mean')
            df_band = df_band.drop_duplicates('mjd')
            df_band = filter_rmsd(df_band)

            if len(df_band) < min_observation_points:
                print(f"Insufficient observations for {band}, skipping...")
                continue

            # Calculate F_var
            f_var, error_f_var = calculate_f_var(df_band)
            if f_var <= 0:
                print(f"F_var <= 0 for {band}, skipping...")
                continue

            # Append to histogram list
            f_var_hist.append(f_var)

            # Assign results to the correct band
            if band == 'ztf_g':
                f_var_g, error_g = f_var, error_f_var
            elif band == 'ztf_r':
                f_var_r, error_r = f_var, error_f_var
            elif band == 'ztf_i':
                f_var_i, error_i = f_var, error_f_var

        except Exception as e:
            print(f"Error processing {band}: {e}")
            continue

    # Save results to a consolidated table
    result_row = [ra, dec, f_var_g, f_var_r, f_var_i, error_g, error_r, error_i]
    save_results_to_table(result_row, results_table_path)

    return result_row

def save_results_to_table(result_row, results_table_path):
    """Save results to a consolidated table."""
    # Create DataFrame if the file doesn't exist
    if not os.path.exists(results_table_path):
        df = pd.DataFrame(columns=['RA', 'DEC', 'F_var_g', 'F_var_r', 'F_var_i', 'Error_g', 'Error_r', 'Error_i'])
    else:
        df = pd.read_csv(results_table_path)

    # Append the new row and save the table
    new_row = pd.DataFrame([result_row], columns=df.columns)
    df = pd.concat([df, new_row], ignore_index=True)
    df.to_csv(results_table_path, index=False)

def plot_fvar_statistics(results_table_path, output_directory):
    """
    Reads the variability results CSV, plots histograms for F_var for g, r, and i bands,
    and displays basic statistics.

    Parameters:
    - results_table_path: Path to the variability results CSV.
    - output_directory: Directory to save the output histogram.

    Outputs:
    - Saves the histogram as 'fvar_statistics.png' in the output directory.
    """
    # Read the results table
    if not os.path.exists(results_table_path):
        print(f"Error: Results table not found at {results_table_path}")
        return

    df = pd.read_csv(results_table_path)

    # Extract F_var columns
    f_var_g = df['F_var_g'].dropna()
    f_var_r = df['F_var_r'].dropna()
    f_var_i = df['F_var_i'].dropna()

    # Compute statistics for each band
    stats = {}
    for band, f_var in zip(['g', 'r', 'i'], [f_var_g, f_var_r, f_var_i]):
        stats[band] = {
            'mean': np.mean(f_var),
            'median': np.median(f_var),
            'std': np.std(f_var),
            'count': len(f_var),
        }

    # Plot histograms
    plt.figure(figsize=(10, 6))
    plt.hist(f_var_g, bins=30, alpha=0.7, color='green', label=f'g-band (n={stats["g"]["count"]})')
    plt.hist(f_var_r, bins=30, alpha=0.7, color='red', label=f'r-band (n={stats["r"]["count"]})')
    plt.hist(f_var_i, bins=30, alpha=0.7, color='orange', label=f'i-band (n={stats["i"]["count"]})')
    plt.xlabel(r'$F_{var}$', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title(r'Histogram of $F_{var}$ for ZTF Bands', fontsize=16)
    plt.legend()
    plt.grid(alpha=0.6)
    plt.tight_layout()


    # Save the histogram
    output_path = os.path.join(output_directory, 'fvar_statistics.png')
    plt.savefig(output_path)

    # Print statistics
    print("\nStatistics for F_var:")
    for band in ['g', 'r', 'i']:
        print(f"{band}-band:")
        print(f"  Mean: {stats[band]['mean']:.4f}")
        print(f"  Median: {stats[band]['median']:.4f}")
        print(f"  Std Dev: {stats[band]['std']:.4f}")
        print(f"  Count: {stats[band]['count']}")
        print("-" * 30)
