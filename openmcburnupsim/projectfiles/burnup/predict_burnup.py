import joblib
import pandas as pd

model = joblib.load("burnup_predictor.pkl")

new_spectrum = pd.read_csv("new_pebble_spectrum.csv").values
predicted_burnup = model.predict(new_spectrum)

print(f"Predicted Burnup: {predicted_burnup[0]:.2f} MWd/kgU")

