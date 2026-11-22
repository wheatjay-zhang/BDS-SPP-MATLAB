# BDS Single Point Positioning (SPP) with MATLAB

A MATLAB implementation of BeiDou Navigation Satellite System (BDS) Single Point Positioning (SPP) using broadcast ephemeris and dual-frequency pseudorange observations.

## 1. Introduction

This project implements a BeiDou SPP algorithm using MATLAB. It processes RINEX observation and navigation files to calculate satellite positions and estimate receiver coordinates through weighted least squares adjustment.

The program includes several important error correction models and positioning procedures.

## 2. Main Features

The main features and corrections implemented in this project include:

- **Dual-frequency ionosphere-free combination:** Uses BDS B1I/B3I pseudorange observations to eliminate the first-order ionospheric delay.
- **Satellite clock correction:** Applies satellite clock bias and relativistic correction.
- **Signal transmission time correction:** Calculates satellite transmission time using pseudorange and satellite clock information.
- **Earth rotation correction:** Corrects satellite coordinates for Earth rotation during signal propagation using an iterative procedure.
- **Satellite orbit calculation:** Computes satellite coordinates from broadcast ephemeris, supporting GEO and MEO satellites.
- **Tropospheric delay correction:** Applies the Saastamoinen model with a simple elevation-dependent mapping function.
- **Weighted least squares adjustment:** Estimates receiver coordinates and receiver clock offset using elevation-dependent weights.
- **Outlier rejection:** Iteratively removes satellites with large post-fit residuals.

## 3. Requirements

- MATLAB with support for RINEX reading and coordinate conversion functions.
- Required MATLAB functions include `rinexread`, `ecef2lla`, `timetable2table`, `rmmissing` and `writetable`.
- BDS dual-frequency pseudorange observations: C2I and C6I.
- BDS broadcast navigation data.

The availability of required functions depends on the MATLAB version and installed toolboxes.

## 4. File Structure

| File / Function     | Description                                                                                                         |
| ------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `main.m`            | Main program: data reading, observation processing, satellite positioning, SPP solution and plotting                |
| `approx_position.m` | Reads approximate receiver coordinates from the RINEX observation file header                                       |
| `GEO_STATE_tsv.m`   | Calculates BeiDou GEO satellite coordinates from broadcast ephemeris                                                |
| `MEO_STATE_tsv.m`   | Calculates BeiDou MEO satellite coordinates from broadcast ephemeris                                                |
| `least_squares.m`   | Performs weighted least squares positioning and residual-based satellite rejection                                  |
| `calculate_B_L.m`   | Constructs the design matrix and observation residual vector, including geometric range and atmospheric corrections |

*Note: The functions are currently implemented as local functions at the end of the main MATLAB script. The table describes their logical modules; they do not necessarily exist as separate files.*

## 5. Input Data

The program requires two RINEX files:

- Observation file (`.rnx`): BDS dual-frequency pseudorange observations, including C2I and C6I.
- Navigation file (`.rnx`): BDS broadcast ephemeris.

The input filenames are specified at the beginning of the main program.

## 6. Output

The program generates the following files:

- `AAAbeidou_coords.csv`: Calculated satellite coordinates and satellite clock corrections.
- `AAAresults.csv`: Estimated receiver coordinates and iteration counts for each epoch.
- `AAAremoved_satellites.csv`: Records of satellites removed due to large residuals.

It also generates a plot showing coordinate differences between the estimated positions and the approximate receiver coordinates.

## 7. Usage

1. Download or clone this repository.
2. Open the project in MATLAB.
3. Specify the RINEX observation and navigation filenames in the main program.
4. Ensure the required MATLAB functions and toolboxes are available.
5. Run the main script.

## 8. Notes

This project is intended for GNSS learning and research purposes. The implementation and positioning results may require further validation and improvement.

## 9. License

License to be determined.
