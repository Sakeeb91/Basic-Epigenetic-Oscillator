#!/usr/bin/env python3
"""
Script to open the epigenetic oscillator notebook in Google Colab.

This script can be executed directly or imported and its functions used.
"""

import os
import webbrowser
from urllib.parse import quote

def open_in_colab(notebook_path, github_username="USERNAME", github_repo="Basic-Epigenetic-Oscillator", branch="main"):
    """
    Open a Jupyter notebook in Google Colab.
    
    Args:
        notebook_path (str): Path to the notebook relative to the repository root.
        github_username (str): GitHub username that hosts the repository.
        github_repo (str): Name of the GitHub repository.
        branch (str): Repository branch (default is "main").
        
    Returns:
        str: URL to the notebook in Google Colab.
    """
    # Construct the GitHub URL to the notebook
    github_path = f"{github_username}/{github_repo}/{branch}/{notebook_path}"
    
    # Construct the Colab URL
    colab_url = f"https://colab.research.google.com/github/{github_path}"
    
    # Open the URL in a browser
    webbrowser.open(colab_url)
    
    return colab_url

if __name__ == "__main__":
    # Path to the notebook relative to the repository root
    notebook_path = "notebooks/epigenetic_oscillator.ipynb"
    
    # Open the notebook in Colab
    url = open_in_colab(notebook_path)
    print(f"Opening notebook in Colab: {url}") 