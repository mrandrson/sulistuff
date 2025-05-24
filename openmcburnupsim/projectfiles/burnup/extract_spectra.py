import openmc
import openmc.deplete
import pandas as pd
import numpy as np
import glob

def extract_data(statepoint_file, depletion_file, output_csv):
    h5_files = sorted(glob.glob("openmc_simulation_n*.h5"))  
    if not h5_files:
        raise FileNotFoundError("No OpenMC statepoint files found.")

    latest_statepoint = h5_files[-1]  
    print(f"Using statepoint data from: {latest_statepoint}")

    # Load neutron and gamma flux from the statepoint file
    sp = openmc.StatePoint(latest_statepoint)

    neutron_tally = sp.get_tally(name='Neutron flux')
    gamma_tally = sp.get_tally(name='Gamma spectrum')

    neutron_df = neutron_tally.get_pandas_dataframe()
    gamma_df = gamma_tally.get_pandas_dataframe()

    neutron_flux = neutron_df.groupby("energy low [eV]")["mean"].sum().values
    gamma_flux = gamma_df.groupby("energy low [eV]")["mean"].sum().values

    # Load burnup data from depletion_results.h5
    try:
        results = openmc.deplete.Results(depletion_file)
        burnup_time, u235_atoms = results.get_atoms("1", "U235")
        burnup = burnup_time[-1]  # Final burnup value
        print(f"Final Burnup Level: {burnup:.2f} MWd/kgU")
    except Exception as e:
        print(f"Warning: Could not load depletion data ({e})")
        burnup = np.nan  # Assign NaN if burnup data is missing

    df = pd.DataFrame([np.concatenate((neutron_flux, gamma_flux, [burnup]))])
    df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    extract_data("statepoint.100.h5", "depletion_results.h5", "spectral_data.csv")

