import itertools
import subprocess

enrichments = [3.0, 5.0, 10.0]
powers = [5e6, 10e6, 20e6]
burn_times = [[1, 10, 30, 50], [5, 20, 40, 60]]

for enrich, power, times in itertools.product(enrichments, powers, burn_times):
    output_dir = f"burnup_results_{enrich}_{power}"
    subprocess.run(["python", "burnup_simulation.py", str(enrich), str(power), str(times), output_dir])
    subprocess.run(["python", "extract_spectra.py", f"{output_dir}/depletion_results.h5", "full_spectral_dataset.csv"])

