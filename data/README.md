# Intersection-point-height data

Contains intersection-point-height (zIP) data from unimpaired subjects from experiments by [dos Santos et al. (2017)](https://doi.org/10.7717/peerj.3626).

## Data file content

Each data file contains the following data:

    DataInfo  - subject information

    Frequency - frequency points at which the zIP is evaluated

    IP        - zIP evaluated at each frequency point (row), separately for each subject and each trial (column);
                For Bartloff2024_paretic and Bartloff2024_nonparetic, data from the corresponding limb is included;
                For dosSantos2017_old, zIP of each limb is included separately.

    IPDataAverage - zIP at each frequency point averaged over subjects and trials

    MeanHeight_m  - average height of all subjects in data set, in meters

    MeanMass_kg   - average mass of all subjects in data set, in kilograms

    NumSubjects   - number of subjects in data set

    Plane         - plane in which zIP was computed; 'sgt' = sagittal plane; 'frt' = frontal plane

    Pose          - pose of subjects, used to compute the inertial properties of the double-inverted-pendulum model;
                    'pose_I' = arms hanging by the side;
                    'pose_T' = arms out to the side, forming a T-pose in the frontal plane

    StandardDeviation - standard deviation of the average zIP data over subjects and trials, IPDataAverage

    SubjectType       - details on data group


## References:

dos Santos, D. A., Fukuchi, C. A., Fukuchi, R. K., & Duarte, M. (2017). A data set with kinematic and ground reaction forces of human balance. PeerJ, 2017(7), 1–17. https://doi.org/10.7717/peerj.3626