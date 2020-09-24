# Stripointing
Simulation of Sun, Moon and planets observation time during the LSPE/Strip microwave telescope activity, based on the Stripeline.jl simulator (https://lspestrip.github.io/Stripeline.jl/latest/).

**Simulation**:
The function "Simulation" produces a .txt file for a celestial body and a horn in a certain time period specified by the user.

**Getting the results**
This is divided in two parts:
1) **Time to be discarded due to presence of Sun and Moon**. Sun and Moon emit whith high intensity in the microwave wavelenghts disturbing the observations for several reasons. The first task was to calculate the total time of observation one must discard due to Sun and/or Moon presence within a certain angular distance from one (or more than one) horn.
Use the function "Get_Results" to extract the results from a list of .txt files relative to different horns for Moon or Sun. 

2) **Planets observation time**. Solar system's planets are useful for Strip project in several ways, such as their fundamental role in the calibration procedure of the instrument. Use the function "Planets_Results" to extract the results from a list of .txt files relative to different horns for a certain planet.
