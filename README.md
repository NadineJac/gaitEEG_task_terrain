# gaitEEG_task_terrain
Analysis scripts for: Jacobsen N.S.J, Blum, S., Scanlon, J.E.M.,Witt,K. &amp; Debener S., 
Mobile EEG captures differences of walking over even and uneven terrain but not of single and dual-task gait 
(submitted)[Preprint link follows]
Please start by opening the master script naj_neurCorGait_master.m and download the required data and toolboxes mentioned in the header.
Then you can run the analysis from the master script. 

The spectral cleaning employed by specPCAdenoising.m (called by task_sPCA_sensor.m) is based on scrips by Martin Seeber. 
It performs a PCA of the gait ERSP and removes the first first principal component to attenuate muscle artifacts as describend in:
Seeber M, Scherer R, Wagner J, Solis-Escalante T, Müller-Putz GR. 
High and low gamma EEG oscillations in central sensorimotor areas are conversely modulated during the human gait cycle. 
Neuroimage. 2015 May 15;112:318–26. 
Available from: https://linkinghub.elsevier.com/retrieve/pii/S105381191500227X
