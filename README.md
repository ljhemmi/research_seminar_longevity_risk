To run the code, download the full repository folder as a zip file from GitHub. 

Within this folder open the Model.Rproj file that will open a specific R-Studio project environment. 
Within the file explorer in R-Studio click on longevity_risk_pricing.R and run the desired code. 
Ensure that all relevant libraries are installed.

The code is structured into the following sub-sections:

  0.	Setup containing libraries and relevant data handling and simulation functions
  1.	Data import of mortality data
  2.	Manual Lee-Carter fitting
  3.	Exemplary Lee-Carter mortality forecasting of mortality index
  4.	Pricing of S-forward

Sections 2 and 3 are not required to perform the actual numerical calculations since all manual steps can also directly be performed using the authors’ functions “mx_to_LC_model()”, “LC_to_simulation_optimized()” and “pricing_calculation()”. 
In the code, sections 2 and 3 are primarily used to create data visualizations for which the manual approach is better suited. 
To recreate out final numerical results in Table 1 it is sufficient to run sections 0, 1 and 4.
