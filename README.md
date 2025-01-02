# Tesla Coil Simulation 

## Overview

This software simulates the behavior of a Tesla Coil system over a specified time span using C++ Object-Oriented Programming 
(OOP) approach. The interactions between its components, including the power supply, transformer, spark gap, capacitor, and coils, are modeled using mathematical principles and numerical methods. The results are saved to files for analysis, allowing an in-depth understanding of the system dynamics.

---

## Features

- Models interactions between Tesla Coil components.
- Simulates the charging and discharging behavior of capacitors.
- Calculates voltages, currents, and resistances dynamically.
- Exports detailed simulation data to CSV files for analysis.

---

## Simulation Details

- The simulation runs over a user-specified time duration with a fixed time step.
- Output files include:
  1. **All Data File** (`allData.csv`):
   - Comprehensive data, including all component parameters over time.
  2. **Voltage-Time Data File** (`voltTimeData.csv`):
   - Time-series data for voltages across key components.
---

## Notes

- The time step (`dt`) should be sufficiently small for accurate simulation results.
- Given the specific operational requirements of a Tesla Coil, the initial values for each component in the Tesla coil simulation have been carefully selected to ensure a functioning model; however, you're free to experiment with different values to explore various behavior and outcomes of the simulation.

---

Feel free to contribute or adapt this project based on your specific needs!
