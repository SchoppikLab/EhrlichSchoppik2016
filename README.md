# EhrlichSchoppik2016
by David Ehrlich and David Schoppik
Department of Otolaryngology, Department of Neuroscience and Physiology, and the Neuroscience Institute
New York University School of Medicine, New York, NY 10016

PMID: 28111151
PMCID: PMC5421408
DOI: 10.1016/j.cub.2016.12.003

Please contact us with questions at ehrlichde@gmail.com or schoppik@gmail.com.

I. Description of contents and directory architecture.
II. Instructions to generate figures from preprocessed data.
III. Instructions to generate figures from raw data.
IV. Variable glossary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I. Description of contents and directory architecture.

The provided 'Raw_Behavior_Data' directory contains folders corresponding to groups: one for each developmental time-point (4, 7, 14, and 21 days post-fertilization), one for larvae with oil-filled swimbladders ('Control_swimbladders_05dpf') and one for their control siblings ('Mineral_oil_swimbladders_05dpf'). Within each group folder is a subfolder for each clutch (siblings), containing the raw data (tab delimited positions and orientations of fish with timestamps) from 48 consecutive hours of behavioral recording. To extract bouts from these groups for subsequent analysis, three similar analysis files are provided: 'Preprocess_By_Age', 'Preprocess_By_Clutch', and 'Preprocess_bouts_by_clutch_oil'. These functions have been run as described in Section III to generate the folder, 'Processed_Behavior_Data'.

'Preprocess_By_Age' iterates through each clutch subfolder to detect bouts, then analyzes all bouts from a given group collectively and saves data in the group folder.
'Preprocess_By_Clutch' iterates through each clutch subfolder, analyzes swim bouts for that clutch, and saves clutch-specific data in the subfolder for that clutch.
'Preprocess_bouts_by_clutch_oil' performs identical analysis but excludes bouts pointed vertically down. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II. To generate figures from preprocessed data:

Add Called_functions subfolder to Matlab path.

In the Processed_Behavior_Data directory run:
Make_Figure_1.m
Make_Figure_2.m
Make_Figure_3.m
Make_Figure_4.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
III. To generate the figures from raw data, execute the following steps in order:

Add EhrlichSchoppik2016 directory and Called_functions subfolder to Matlab path.

To compute swimming statistics and extract and save swim bout data for typically-developing larvae:
-In each developmental folder (04dpf, 07dpf, 14dpf, and 21dpf) run 'Preprocess_bouts_by_age.m' (this may take several minutes)
-Then, in each developmental folder (04dpf, 07dpf, 14dpf, and 21dpf) call 'Preprocess_bouts_by_clutch' (this may take several minutes)

To generate Figures 1 and 2, in the Raw_Behavior_Data directory run:
Make_Figure_1.m
Make_Figure_2.m

To generate Figure 4, begin by computing relative bout likelihood:
-In each developmental folder run 'Compute_bout_likelihood.m' (this may take tens of minutes).
-Then, in each developmental folder run 'Stochastic_swim_fit.m'
(this may take tens of minutes).
-Then, in the Raw_Behavior_Data directory run Make_Figure_4.m

To generate Figure 3, compute swimming statistics and extract and save swim bout data for larvae with oil-filled swimbladders and controls (this may take several minutes): 
-In folders 'Control_swimbladders_05dpf' and 'Mineral_oil_swimbladders_05dpf' run 'Preprocess_bouts_by_clutch_oil'
To generate figure 3:
-Then, in the Raw_Behavior_Data directory run Make_Figure_3.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IV. Variable glossary

Swimming variables are saved in individual group and clutch folders (see Section I) and are defined in preprocessing code where created or filled. Brief definitions of these variables follow.

	AVedges: Octiles of empirical, instantaneous angular velocities, used for calculating ProbAllObs and BoutProbability in Compute_bout_likelihood.m
	BodyAngles: Concatenated, pitch-axis postures (in deg) from all observed and non-excluded frames (0 is horizontal, positive is nose-up to horizontal).
	BoutProbability: number of observations of a given (posture, angular velocity) pair immediately before bouts, divided by the total observations of that pair when larvae were eligible to initiate a bout (excluding ongoing bouts and absolute refractory periods)
	coeffInit: initial values of parameters (alpha, zeta, beta) for fitting relative bout likelihood
	coeffs: estimates for parameters (alpha, zeta, beta, gamma) defining relative bout likelihood
	CovB: variance-covariance matrix for estimated relative bout likelihood coefficients, used to calculate confidence intervals
	GrabbedTimes: concatenated time-stamps from all observed and non-excluded frames
	HeadingMatchedAngles: observed pitch-axis postures corresponding to every instantaneous heading (translation direction)
	HeadingMatchedAngVels: observed pitch-axis angular velocities (smoothed by moving average with a span of 50 msec) corresponding to every instantaneous heading (translation direction)
	HeadingMatchedSpeeds: translation speeds corresponding to every instantaneous heading (translation direction)
	HeadingMatchedTimes: timestamps corresponding to every instantaneous heading (translation direction)
	PitchEdges: Octiles of empirical, instantaneous posture (in pitch axis), used for calculating ProbAllObs and BoutProbability in Compute_bout_likelihood.m
	ProbAllObs: Total observations of a given (posture, angular velocity) pair when larvae were eligible to initiate a bout (excluding ongoing bouts and absolute refractory periods).
	ProbPreBout: Number of observations of a given (posture, angular velocity) pair immediately before bouts.
	PropBoutAlignedAngVel: pitch-axis angular velocity during swim bouts for 750 msec preceding and 500 msec succeeding the instance of maximal translation speed
	PropBoutAlignedPitch: pitch-axis posture during swim bouts for 750 msec preceding and 500 msec succeeding the instance of maximal translation speed
	PropBoutAlignedSpeed: translation speed during swim bouts for 750 msec preceding and 500 msec succeeding the instance of maximal translation speed
	PropBoutAlignedTime: timestamps of observed swim bouts flanked with 750 msec of preceding and 500 msec of succeeding postural data, in hours
	PropBoutDisplacement: displacement across 600 msec temporally centered on each swim bout 
	PropBoutDuration: interpolated duration above 5 mm/sec translation speed for each swim bout
	PropBoutIEI: duration (in sec) between pairs of successive swim bouts
	PropBoutIEIalignedAngVel: cell array of pitch-axis angular velocities (smoothed by moving average with a span of 50 msec) between pairs of successive bouts, from 50 msec following to 100 msec preceding the instances of peak speed of the first and second bout, respectively
	PropBoutIEIalignedPitch: cell array of pitch-axis posture between pairs of successive bouts, from 50 msec following to 100 msec preceding the instances of peak speed of the first and second bout, respectively
	PropBoutIEIangAcc: mean pitch-axis angular acceleration during a window spanning from 300 msec following the instance of peak speed of one bout to 100 msec preceding the instance of peak speed for the next bout
	PropBoutIEIangVel: mean pitch-axis angular velocity during a window spanning from 300 msec following the instance of peak speed of one bout to 100 msec preceding the instance of peak speed for the next bout
	PropBoutIEIpitch: mean pitch-axis posture during a window spanning from 300 msec following the instance of peak speed of one bout to 100 msec preceding the instance of peak speed for the next bout
	PropBoutIEItime: timestamps of observed bout IEIs in hours
	PropBoutInitPitch: mean pitch-axis posture from 250 to 125 msec preceding instance of maximal speed during each swim bout
	PropBoutMaxAngVel: maximal angular speed while translation speed exceeds 5 mm/sec during each swim bout
	PropBoutMaxSpd: maximal translation speed during each swim bout
	PropBoutNetPitchChg: difference between mean pitch-axis postures from 250-125 msec preceding and following instance of maximal speed during each swim bout
	PropBoutPeakAngVel: largest magnitude angular velocity while translation speed exceeds 5 mm/sec during each swim bout
	PropBoutTime: timestamps of observed swim bouts in hours
	PsiCovB: variance of estimated coefficient of prior (psi) for bout initiation probability, used to calculate confidence intervals 
	PsiFit: estimates for parameter (psi) defining prior for bout probability
	PsiR: residuals from estimating coefficient of prior (psi) for bout initiation probability, used to calculate confidence intervals 
	R: residuals from estimating relative bout likelihood coefficients, used to calculate confidence intervals
