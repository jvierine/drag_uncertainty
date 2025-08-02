# The effect of atmospheric drag uncertainty on orbital elements

The script uses the SGP4 propagator to investigate how errors in $B^* = \rho_0 A \frac{C_d}{2m}$ affect position uncertainty. It doesn't matter if the uncertainty in the area-to-mass ratio ($A/m$) of the spacecraft, the drag coefficient ($C_d$), or the atmospheric density ($\rho_0$). They all have the same effect on orbit propagation.

The following figure shows how the position error accumulates as a functino of time when there is an error in the $B^*$ term:

<img width="595" height="438" alt="Screenshot 2025-08-01 at 10 46 32" src="https://github.com/user-attachments/assets/950ab753-f93f-4f29-95da-fb7e81e84890" />

Thermospheric density uncertainty typically several tens of percent and depends orbital parameters as well as geophysical forcing (Boniface and Bruinsma, 2021; Liying and Solomon 2012). 

This information can now be used to answer the question of how does the collision probability develop as a function of time due to position uncertainty increasing because of uncertainty in satellite drag:

<img width="587" height="430" alt="Screenshot 2025-08-02 at 10 37 54" src="https://github.com/user-attachments/assets/572798cf-a092-44d9-9be3-ff991cfef7be" />

Boniface, Claude, and Sean Bruinsma. "Uncertainty quantification of the DTM2020 thermosphere model." Journal of Space Weather and Space Climate 11 (2021): 53.

Qian, Liying, and Stanley C. Solomon. "Thermospheric density: An overview of temporal and spatial variations." Space science reviews 168.1 (2012): 147-173.

https://en.wikipedia.org/wiki/BSTAR
