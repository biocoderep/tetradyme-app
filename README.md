# Tetradyme - ARMS-PCR Designer

Tetradyme is a web-based tool for designing inner primers for Tetra-Primer ARMS-PCR based on SNP information and target melting temperatures.

## Architecture

This project has been migrated from Streamlit to a **FastAPI backend** with a **Vanilla HTML/JS/CSS frontend**. 
The core logic has been isolated into `core.py`, making the application highly modular, extensible, and performant.

### Core Logic Enhancements
- **Expanded Primer Search Space**: The algorithm now evaluates primer lengths from 15 to 35 bp (up from the restrictive 18-25 bp in the legacy app). This drastically improves the design success rate.
- **Graceful Tm Tolerance Fallback**: Instead of outright failing if no primer combination strictly meets the Tm tolerance, the algorithm now returns the *best available* primer pair if it falls within a slight relaxed tolerance window (+2.0°C). This ensures users still get viable candidate primers for difficult SNPs instead of a complete design failure.

## Prerequisites

- Python 3.9+ 
- Or Docker (optional)

## Setup & Run (Local Development)

1. **Install dependencies**
   Open your terminal (e.g. PowerShell or Command Prompt), navigate to the project directory, and run:
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the server**
   Use Python to run the Uvicorn server:
   ```bash
   python -m uvicorn main:app --port 8000
   ```
   *(Note: If you have uvicorn in your system PATH, you can also just run `uvicorn main:app --port 8000`)*

3. **Open the app**
   Once the terminal says "Application startup complete", open your web browser and navigate to:
   [http://127.0.0.1:8000](http://127.0.0.1:8000)

## Run via Docker

1. **Build the image**
   ```bash
   docker build -t tetradyme-app .
   ```

2. **Run the container**
   ```bash
   docker run -p 8000:8000 tetradyme-app
   ```
   The application will be accessible at [http://localhost:8000](http://localhost:8000).

## Usage

1. Upload your FASTA file.
2. Upload your SNP Data CSV file (must include columns: `ID, pos, ref, alt, outer_tm`).
3. Set your Reaction Conditions (Na+, Mg2+, dNTPs).
4. Specify Tm tolerance and optional Outer Primers.
5. Click **Design Primers**.
6. View the summary, detailed results table, gel electrophoresis visualization, and download the results CSV.
