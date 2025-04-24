This file contains my code created to calculate the Hubble-Constant using Gravitational wave events independently without any bias.
This is my Master's project.

The sky localisation of all the GW events is in the skymap folder.

distance.ipynb is a code that can find the distance of GW without relying on LIGO, by using Bayesian inference and the SNR of the event to get high and low priors.

statistical_redshift.ipynb calculates the statistical redshift of the events based on their sky localisation map and the galaxy catalogue GLADE+ or SDSS.
Note that the SDSS survey doesn't have distances on its own, so distances are calculated using H0(CMB) and H0(Ceph), and an average of the measurements could be done to reduce the bias.

GLADE+ galaxy catalogue can be downloaded at the following link (It is 6Gb+): GLADE.elte.hu 

SDSS catalogue can be downloaded here: https://www.skyserver.sdss.org/dr18/SearchTools/sql#
6 SQL commands need to be done to retrieve the whole galaxy dataset. 
SQL command is shown here: 
SELECT TOP 500000 ra, dec, z, zerr, SpecObjid
FROM SpecObj
WHERE class = 'GALAXY' 
AND z BETWEEN 0 AND 1
AND SpecObjid NOT BETWEEN 2.9949E+17 AND 1.80157E+18
ORDER BY SpecObjid;

Don't forget to change the specObjid to get different galaxies and not get the same.


Hubble_constant.ipynb calculate the H0 based on different criteria (GLADE+ <4000Mpc, SDSS(CMB, Ceph, Mean) <4000Mpc and GLADE+ <2500 Mpc with SDSS(Mean) 2500Mpc<d<4000 Mpc.
A cumulative mean for different events (up to 23) is plotted in a graph with previous measurements using CMB and Ceph to compare.  
