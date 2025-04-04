#!/bin/bash

# Create necessary directories
mkdir -p plots
mkdir -p sample_outputs

# Install dependencies
pip install -r requirements.txt

# Run the simulation to generate sample plots
python src/epigenetic_oscillator.py

echo "Setup complete! You can now explore the simulation in the notebook: notebooks/epigenetic_oscillator.ipynb" 