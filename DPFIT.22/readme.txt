This is my first attempt to improve a TLE using Doppler shift data.  Results look promising, but I would really like some feedback.  

OPERATION:

1.  Extract the contents of the DPFIT.ZIP file into a directory.

2.  Double click on the DPFIT.EXE file.  You will be looking at a TLE of sat 90027 and observations all listed in the default dp.obs file.

3.  The suggested transmit frequency has been computed as a rounded average of the reported observed frequencies.  If it is not already displayed, type in 256.375 as the most likely guess.  Pressing ENTER will accept the values for the transmit frequencies (6th column) in the input file.  This will allow mixing data lines using different frequencies.  If you press a+ENTER, the displayed average will be used as the transmit frequency for every observation instead of the values in the input file.  The default name for the input file is dp.obs.  A message appears stating that 88 observations have been found in the input file (dp.obs when you double click).

4.  The epoch of the TLE may be changed by pressing t+ENTER.  Press ENTER to accept all the prompt values.  These default values are the exact time of the ascending node just prior to the last observed time.

5.  The transmitter offset can be calculated by pressing o+ENTER and a+ENTER.  The offset is -16Hz.  Press ENTER again and the (F)it is displayed.  You will see observation number, station number, STA, azimuth, AZ, elevation, EL, aspect angle, ASP, observed doppler shift, OBS, computed doppler shift, CMP, and absolute value of percent difference scaled to transmitter frequency, %Diff.  Press ENTER four times to page through the observations.  Q+ENTER to quit this function.

6.  If the satellite orbit is unknown, it is best to find true anomaly.  If the orbit is being updated, skip this step.  

Press a+ENTER then s+ENTER to search for the optimum true anomaly.  Accept the default search limits, 360 to 0, and 20 steps by pressing ENTER 3 more times.  The rms sum is computed for each value of the true anomaly.  It looks like 18 is closest. Type in 18+ENTER and scroll through the (F)it output.  It appears that all the computed elevations, EL, are all well above the horizon, so we are on the right track.  Now select a+ENTER for an automatic search.  The final best value is 14.7598.  q+ENTER to quit this function. 

7.  The (S)earch function optimizes all six orbit elements in a TLE.  Press s+ENTER. Continue to press ENTER to continue optimizing as indicated by the decreasing rms sum.  When satisfied, press q+ENTER to quit.  Page through the new (F)it.  The doppler residuals are further reduced.

8.  To save these results, press w+ENTER.  Press u+ENTER to overwrite the old TLE with the updated TLE or press a+ENTER to append the new TLE to the bottom of the input/output file, dp.obs in this case.

9.  Press q+ENTER twice to quit the program.

Try opening DPFIT.EXE by running DPFIT.EXE ivan.obs.  The operation is the same, but calculations will take longer due to the larger number of observations.

DPFIT will read the observation format as per dp.obs.  This is the standard format for radio frequency observations.  Also note that STATIONS.IN needs to be in the same directory with DPFIT with all appropriate observer coordinates included.

Manual update of individual elements is possible with DPFIT in exactly the same manner as with SATFIT.  See the SATFIT page for more details of these operations, http://www.coastalbend.edu/acdem/math/sats/satfit.htm

10.  I have added the (U)Simulate function to model a satellite as seen by a single observer.  Press u+ENTER to begin a simulation.  Enter the observer site code and the start time, stop time, and step size of the simulation. (P)ass will only show times when the satellite is above the observers horizon. 

11. Satellites in Molniya orbits experience doppler maximums and minimums.  The files 90025.txt and 32378.txt are included to practice these capalilities.  The format for reporting times of observed Doppler min and max is as per the examples in 90025.txt.  DPFIT will recongnize the format type and be ready to optimize parameters based on the min/max oservations.

Another DPFIT format is shown in 32378.txt.  This format is for both Doppler Min/Max times and radio AOS/LOS times.  The optimization features only use the Min/Max times.  The AOS/LOS times are included to aid manual adjustment of the elements.


Scott Campbell
campbel7@hughes.net