J/A+A/vol/page  Asteroseismic modelling of 36 β Cep stars  (Vanrespaille+, 2026)
================================================================================
Asteroseismic forward modelling of 36 β Cep pulsators and inferences on their 
internal differential rotation
     Vanrespaille M., Fritzewski D.J., Vanlaer V., Aerts C.
    <Astron. Astrophys. vol, page (year)>
    =2003A&A...VVV.ppppI     (SIMBAD/NED BibCode)
================================================================================
Keywords_ADC: Stars, ages ; Stars, OB ; Stars, variable ; Stars, normal
Keywords: Asteroseismology – Stars: oscillations (including pulsations) – 
	  Stars: massive – Stars: interiors – Stars: evolution –
	  Stars: rotation

Abstract: 
    Asteroseismic observations of the interior rotation of main-sequence stars 
    have shown that angular momentum transport is much more efficient than 
    expected. What transport mechanisms are responsible for this is still 
    unclear. The massive {beta} Cep pulsators are promising stars to constrain 
    these mechanisms as they are the only main-sequence pulsators to regularly 
    display differential rotation. However, fewer than ten {beta} Cep stars 
    have been asteroseismically modelled so far.

    We aim to expand the sample of asteroseismically forward modelled {beta} Cep 
    pulsators to maximally exploit their potential to observationally constrain
    angular momentum transport mechanisms. To that end, we seek measurements of 
    differential rotation.
    
    We searched for rotational splitting in a large {beta} Cep sample with 
    partial mode identifications. These are subjected to a novel forward 
    modelling approach that consistently takes second-order rotation effects 
    into account.

    36 {beta} Cep stars have successfully been modelled, constraining their 
    mass, age, rotation frequency, core overshoot and envelope mixing, among 
    other parameters. Herein, accounting for second-order effects of rotation 
    significantly affects the rotation estimates when the rotation rate exceeds 
    10% of the critical rotation rate. Like intermediate-mass main-sequence 
    stars, the rotation rate decreases as the stars evolve. Differential 
    rotation is constrained in 17 stars, with the rotation rate in at least 14 
    of them varying by over 10%.

    We show differential rotation is commonplace in {beta} Cep stars. Prevailing 
    features in the differential rotation suggests {beta} Cep typically have a 
    non-monotonic rotation profile.

Description:
    ObsInput.dat: Observational properties of beta Cephei stars observed with
    Gaia and TESS. It includes the effective temperature, luminosity and 
    quantities used to compute that luminosity, the identified signals 
    including their frequency and its uncertainty, amplitude, radial order, 
    degree, and azimuthal order. 
    
    Models.dat: The chosen parameters of the statistical models of each star 
    with a satisfactory fit, including their initial mass, envelope mixing 
    strength, core overshoot parameter, central hydrogen mass fraction, age, 
    optimal rotation frequency, radius, surface gravity, effective 
    temperature, luminosity, and convective core mass. 
    
    RotSplit.dat: The observed and best model's rotational splitting in each 
    identified multiplet used in the asteroseismic modelling. It includes the
    rotational splitting of each observed part of each multiplet as well as 
    the best model's counterpart. If sufficient modes are observed, the 
    asymmetry is also included. Note that there is one line per mode or 
    multiplet with a unique (n,l), not one line per star. 
    
    DiffRot.dat: The measurements of internal and surface rotation used to 
    study differential rotation in Section 6 of the accompanying paper. It 
    includes the Gaia surface rotation frequencies, if available, and the 
    optimal rotation associated with each rotationally split multiplet 
    and the relative radius where each multiplet is most sensitive.  
    
File Summary:
--------------------------------------------------------------------------------
  FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe.txt      80        .   This file
ObsInput.dat 10325       43   The 43 stars considered and their observational 
                               input for 43 stars considered
Models.dat     561       36   Essential parameters of asteroseismic forward 
                               models of 36 stars with a satisfactory fit
RotSplit.dat   428       81   Observed and modelled asymmetric rotational 
                               splitting of each multiplet
DiffRot.dat    594       36   Differential rotation constraints
--------------------------------------------------------------------------------

Byte-by-byte Description of file: ObsInput.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
    1-    9  I9     ---    TIC_ID             ? TESS object identifier
   11-   23  A13    ---    name               ? common name
   25-   43  I19    ---    Gaia_DR3_ID        ? Gaia DR3 identifier
   45-   63  F19.15 deg    ra                 [0.74/359.65]? right ascension
   65-   81  F17.13 deg    dec                [-71.63/68.28]? declination
   83-   99  F17.14 ---    Vmag               [1.97/13.23]? V-band apparent
                                             magnitude
        101  I1     ---    goodFit            ? Whether a satisfactory model was
                                             found
  103-  112  A10    ---    Teff_source        ? source of effective temperature
                                             measurement
  114-  129  F16.14 K      log_Teff           [4.3/4.55]? log10 of effective
                                             temperature
  131-  151  F21.19 K      log_Teff_err       [0.0/0.02]? uncertainty on log10
                                             of effective temperature
  153-  168  F16.14 solLum log_L              [3.72/5.25]? log10 luminosity
  170-  189  F20.18 solLum log_L_err          [0.01/0.04]? uncertainty on log10
                                             luminosity
  191-  208  F18.16 mas    parallax           [0.19/4.8]? Gaia parallax
  210-  229  F20.18 mas    parallax_err       [0.0/0.3]? uncertainty on Gaia
                                             parallax
  231-  248  F18.16 d-1    Gaia_frot_sini     ? ESP-HS projected surface
                                             rotation frequency
  250-  267  F18.16 d-1    Gaia_frot_sini_err [0.04/0.26]? uncertainty on ESP-HS
                                             projected surface rotation
                                             frequency
        269  I1     ---    Nm_IDs             [2/4]? number of identified radial
                                             modes and multiplets
  271-  272  I2     ---    Nf_IDs             [3/15]? number of identified
                                             frequencies
        274  I1     ---    i1                 [1]? index of multiplet this
                                             identified signal is part of
  276-  277  I2     ---    n1                 [-1/3]? radial order of identified
                                             signal
        279  I1     ---    l1                 ? degree of identified signal
  281-  282  I2     ---    m1                 ? azimuthal order of identified
                                             signal
  284-  301  F18.16 d-1    f1                 [3.99/9.74]? frequency of
                                             identified signal
  303-  307  F5.3   d-1    f_err1             [0.0/0.01]? uncertainty on
                                             frequency of identified signal
  309-  327  F19.16 ---    a1                 [0.0/40.8]? amplitude of
                                             identified signal
        329  I1     ---    i2                 [1/2]? index of multiplet this
                                             identified signal is part of
  331-  332  I2     ---    n2                 [-2/3]? radial order of identified
                                             signal
        334  I1     ---    l2                 [1/2]? degree of identified signal
  336-  337  I2     ---    m2                 [-1/1]? azimuthal order of
                                             identified signal
  339-  357  F19.16 d-1    f2                 [3.88/10.37]? frequency of
                                             identified signal
  359-  363  F5.3   d-1    f_err2             [0.0/0.01]? uncertainty on
                                             frequency of identified signal
  365-  383  F19.16 ---    a2                 [0.0/22.3]? amplitude of
                                             identified signal
        385  I1     ---    i3                 [1/2]? index of multiplet this
                                             identified signal is part of
  387-  388  I2     ---    n3                 [-2/3]? radial order of identified
                                             signal
        390  I1     ---    l3                 ? degree of identified signal
  392-  393  I2     ---    m3                 [-1/2]? azimuthal order of
                                             identified signal
  395-  413  F19.16 d-1    f3                 [3.97/11.04]? frequency of
                                             identified signal
  415-  419  F5.3   d-1    f_err3             [0.0/0.01]? uncertainty on
                                             frequency of identified signal
  421-  439  F19.16 ---    a3                 [0.0/27.1]? amplitude of
                                             identified signal
        441  I1     ---    i4                 [1/3]? index of multiplet this
                                             identified signal is part of
  443-  444  I2     ---    n4                 [-2/3]? radial order of identified
                                             signal
        446  I1     ---    l4                 ? degree of identified signal
  448-  449  I2     ---    m4                 [-2/2]? azimuthal order of
                                             identified signal
  451-  468  F18.16 d-1    f4                 [3.8/7.98]? frequency of
                                             identified signal
  470-  474  F5.3   d-1    f_err4             [0.0/0.01]? uncertainty on
                                             frequency of identified signal
  476-  494  F19.16 ---    a4                 [0.0/24.5]? amplitude of
                                             identified signal
        496  I1     ---    i5                 [1/3]? index of multiplet this
                                             identified signal is part of
  498-  499  I2     ---    n5                 [-2/2]? radial order of identified
                                             signal
        501  I1     ---    l5                 ? degree of identified signal
  503-  504  I2     ---    m5                 [-2/2]? azimuthal order of
                                             identified signal
  506-  523  F18.16 d-1    f5                 [3.25/7.65]? frequency of
                                             identified signal
  525-  529  F5.3   d-1    f_err5             [0.0/0.01]? uncertainty on
                                             frequency of identified signal
  531-  551  E21.16 ---    a5                 [0.0/12.7]? amplitude of
                                             identified signal
  553-  554  I2     ---    i6                 ? index of multiplet this
                                             identified signal is part of
  556-  557  I2     ---    n6                 ? radial order of identified
                                             signal
  559-  560  I2     ---    l6                 ? degree of identified signal
  562-  563  I2     ---    m6                 ? azimuthal order of identified
                                             signal
  565-  581  F17.15 d-1    f6                 ? frequency of identified signal
  583-  587  F5.3   d-1    f_err6             ? uncertainty on frequency of
                                             identified signal
  589-  602  E14.9  ---    a6                 ? amplitude of identified signal
  604-  605  I2     ---    i7                 ? index of multiplet this
                                             identified signal is part of
  607-  608  I2     ---    n7                 ? radial order of identified
                                             signal
  610-  611  I2     ---    l7                 ? degree of identified signal
  613-  614  I2     ---    m7                 ? azimuthal order of identified
                                             signal
  616-  633  F18.16 d-1    f7                 ? frequency of identified signal
  635-  639  F5.3   d-1    f_err7             ? uncertainty on frequency of
                                             identified signal
  641-  661  E21.16 ---    a7                 ? amplitude of identified signal
  663-  664  I2     ---    i8                 ? index of multiplet this
                                             identified signal is part of
  666-  667  I2     ---    n8                 ? radial order of identified
                                             signal
  669-  670  I2     ---    l8                 ? degree of identified signal
  672-  673  I2     ---    m8                 ? azimuthal order of identified
                                             signal
  675-  691  F17.15 d-1    f8                 ? frequency of identified signal
  693-  697  F5.3   d-1    f_err8             ? uncertainty on frequency of
                                             identified signal
  699-  712  E14.9  ---    a8                 ? amplitude of identified signal
  714-  715  I2     ---    i9                 ? index of multiplet this
                                             identified signal is part of
  717-  718  I2     ---    n9                 ? radial order of identified
                                             signal
  720-  721  I2     ---    l9                 ? degree of identified signal
  723-  724  I2     ---    m9                 ? azimuthal order of identified
                                             signal
  726-  743  F18.16 d-1    f9                 ? frequency of identified signal
  745-  749  F5.3   d-1    f_err9             ? uncertainty on frequency of
                                             identified signal
  751-  764  E14.9  ---    a9                 ? amplitude of identified signal
  766-  767  I2     ---    i10                ? index of multiplet this
                                             identified signal is part of
  769-  770  I2     ---    n10                ? radial order of identified
                                             signal
  772-  773  I2     ---    l10                ? degree of identified signal
  775-  776  I2     ---    m10                ? azimuthal order of identified
                                             signal
  778-  794  F17.15 d-1    f10                ? frequency of identified signal
  796-  800  F5.3   d-1    f_err10            ? uncertainty on frequency of
                                             identified signal
  802-  823  E22.17 ---    a10                ? amplitude of identified signal
  825-  826  I2     ---    i11                ? index of multiplet this
                                             identified signal is part of
  828-  829  I2     ---    n11                ? radial order of identified
                                             signal
  831-  832  I2     ---    l11                ? degree of identified signal
  834-  835  I2     ---    m11                ? azimuthal order of identified
                                             signal
  837-  853  F17.15 d-1    f11                ? frequency of identified signal
  855-  859  F5.3   d-1    f_err11            ? uncertainty on frequency of
                                             identified signal
  861-  878  F18.16 ---    a11                ? amplitude of identified signal
  880-  881  I2     ---    i12                ? index of multiplet this
                                             identified signal is part of
  883-  884  I2     ---    n12                ? radial order of identified
                                             signal
  886-  887  I2     ---    l12                ? degree of identified signal
  889-  890  I2     ---    m12                ? azimuthal order of identified
                                             signal
  892-  908  F17.15 d-1    f12                ? frequency of identified signal
  910-  914  F5.3   d-1    f_err12            ? uncertainty on frequency of
                                             identified signal
  916-  933  F18.16 ---    a12                ? amplitude of identified signal
  935-  936  I2     ---    i13                ? index of multiplet this
                                             identified signal is part of
  938-  939  I2     ---    n13                ? radial order of identified
                                             signal
  941-  942  I2     ---    l13                ? degree of identified signal
  944-  945  I2     ---    m13                ? azimuthal order of identified
                                             signal
  947-  964  F18.16 d-1    f13                ? frequency of identified signal
  966-  970  F5.3   d-1    f_err13            ? uncertainty on frequency of
                                             identified signal
  972-  984  E13.8  ---    a13                ? amplitude of identified signal
  986-  987  I2     ---    i14                ? index of multiplet this
                                             identified signal is part of
  989-  990  I2     ---    n14                ? radial order of identified
                                             signal
  992-  993  I2     ---    l14                ? degree of identified signal
  995-  996  I2     ---    m14                ? azimuthal order of identified
                                             signal
  998- 1005  F8.6   d-1    f14                ? frequency of identified signal
 1007- 1011  F5.3   d-1    f_err14            ? uncertainty on frequency of
                                             identified signal
 1013- 1026  F14.12 ---    a14                ? amplitude of identified signal
 1028- 1029  I2     ---    i15                ? index of multiplet this
                                             identified signal is part of
 1031- 1032  I2     ---    n15                ? radial order of identified
                                             signal
 1034- 1035  I2     ---    l15                ? degree of identified signal
 1037- 1038  I2     ---    m15                ? azimuthal order of identified
                                             signal
 1040- 1047  F8.6   d-1    f15                ? frequency of identified signal
 1049- 1053  F5.3   d-1    f_err15            ? uncertainty on frequency of
                                             identified signal
 1055- 1068  F14.12 ---    a15                ? amplitude of identified signal
 1070- 5309  A4240  d-1    all_f              ? all detected frequencies
 5311-10325  A5015  ---    all_a              ? all detected amplitudes
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9      ---    TIC_ID         ? TESS object identifier
 11- 23  A13     ---    name           ? common name
 25- 42  F18.15 solMass M              [8.03/29.86]? initial mass
 44- 63  F20.18 solMass M_err          [0.05/1.04]? uncertainty on initial mass
 65- 82  F18.16 cm2.s-1 logD           [1.0/6.0]? log10 of mixing strength at
                                     the base of the envelope
 84-101  F18.16 cm2.s-1 logD_err       [0.5/1.51]? uncertainty on log10 of
                                     mixing strength at the base of the envelope
103-123  F21.19  ---    fov            [0.0/0.04]? overshoot parameter
125-145  F21.19  ---    fov_err        [0.0/0.01]? uncertainty on overshoot
                                     parameter
147-165  F19.17  ---    Xc             [0.04/0.43]? central hydrogen mass
                                     fraction
167-187  F21.19  ---    Xc_err         [0.0/0.05]? uncertainty on central
                                     hydrogen mass fraction
189-208  F20.18  d-1    frot           [0.01/0.6]? rotation frequency measured
                                     when fitting all observed signals
210-230  E21.16  d-1    frot_err       [0.0/0.03]? uncertainty on rotation
                                     frequency measured when fitting all
                                     observed signals
232-252  F21.19  ---    frot_rel       [0.0/0.41]? rotation frequency relative
                                     to the critical rotation frequency measured
                                     when fitting all observed signals
254-275  E22.17  ---    frot_rel_err   [0.0/0.01]? uncertainty on rotation
                                     frequency relative to the critical rotation
                                     frequency measured when fitting all
                                     observed signals
277-294  F18.9   yr     age            [3282211.88/32511686.79]? age since ZAMS
296-314  F19.11  yr     age_err        [87491.68/1499569.31]? uncertainty on age
                                     since ZAMS
316-333  F18.16  K      log_Teff       [4.28/4.57]? log10 of effective
                                     temperature
335-355  F21.19  K      log_Teff_err   [0.0/0.02]? uncertainty of log10 of
                                     effective temperature
357-374  F18.16  solLum log_L          [3.66/5.21]? log10 of luminosity
376-395  F20.18  solLum log_L_err      [0.01/0.04]? uncertainty on log10 of
                                     luminosity
397-414  F18.16  cm.s-2 log_g          [3.56/3.95]? log10 of surface gravity
416-436  F21.19  cm.s-2 log_g_err      [0.0/0.05]? uncertainty on log10 of
                                     surface gravity
438-455  F18.16  solRad log_R          [0.69/1.07]? log10 of stellar radius
457-477  F21.19  solRad log_R_err      [0.0/0.02]? uncertainty of log10 of
                                     stellar radius
479-497  F19.16 solMass Mcc            [1.33/12.96]? convective core mass
499-519  F21.19 solMass Mcc_err        [0.0/0.32]? uncertainty on convective
                                     core mass
521-539  F19.17  ---    Mcc_over_M     [0.16/0.44]? convective core mass
                                     relative to total mass
541-561  F21.19  ---    Mcc_over_M_err [0.0/0.02]? uncertainty on convective
                                     core mass relative to total mass
--------------------------------------------------------------------------------

Byte-by-byte Description of file: RotSplit.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9     ---    TIC_ID           ? TESS object identifier
 11- 23  A13    ---    name             ? common name
 25- 26  I2     ---    n                [-2/4]? radial order
     28  I1     ---    l                ? degree
 30- 47  F18.16 d-1    f_zonal_obs      [3.5/9.14]? observed zonal frequency
 49- 66  F18.16 d-1    f_zonal_model    [3.5/9.14]? modelled zonal frequency
 68- 87  F20.18 d-1    Delta_f_1_obs    [0.0/0.69]? observed frequency
                                       difference between zonal and m=1 modes
 89-108  F20.18 d-1    Delta_f_-1_obs   [0.0/0.74]? observed frequency
                                       difference between zonal and m=1 modes
110-129  F20.18 d-1    Delta_f_2_obs    ? observed frequency difference between
                                       zonal and m=-1 modes
131-150  F20.18 d-1    Delta_f_-2_obs   ? observed frequency difference between
                                       zonal and m=2 modes
152-174  E23.16 ---    A_1_obs          [-0.14/0.21]? observed frequency
                                       difference between zonal and m=-2 modes
176-196  F21.18 ---    A_2_obs          ? observed asymmetry of m=1 and m=-1
                                       modes
198-217  F20.18 d-1    Delta_f_1_model  [0.0/0.6]? observed asymmetry of m=2 and
                                       m=-2 modes
219-238  F20.18 d-1    Delta_f_-1_model [0.0/0.84]? modelled frequency
                                       difference between zonal and m=1 modes
240-259  F20.18 d-1    Delta_f_2_model  ? modelled frequency difference between
                                       zonal and m=1 modes
261-280  F20.18 d-1    Delta_f_-2_model ? modelled frequency difference between
                                       zonal and m=-1 modes
282-303  F22.19 ---    A_1_model        [-0.2/0.68]? modelled frequency
                                       difference between zonal and m=2 modes
305-325  F21.19 ---    A_2_model        ? observed frequency difference between
                                       zonal and m=-2 modes
327-344  F18.16 d-1    frot_StORM       [0.01/0.6]? modelled asymmetry of m=1
                                       and m=-1 modes
346-367  E22.17 d-1    frot_StORM_err   [0.0/0.04]? modelled asymmetry of m=2
                                       and m=-2 modes
369-388  F20.18 d-1    frot_GYRE        [0.01/0.75]? rotation frequency
                                       optimised with StORM including second
                                       order effects
390-409  F20.18 ---    C_nl             ? uncertainty on rotation frequency
                                       optimised with StORM including second
                                       order effects
411-428  F18.16 ---    beta_nl          [0.48/1.0]? rotation frequency estimate
                                       from first order post-processing step
                                       using GYRE
--------------------------------------------------------------------------------

Byte-by-byte Description of file: DiffRot.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9     ---    TIC_ID                     ? TESS object identifier
 11- 23  A13    ---    name                       ? common name
 25- 42  F18.16 d-1    Gaia_frot_sini             ? ESP-HS projected surface
                                                 rotation frequency
 44- 61  F18.16 d-1    Gaia_frot_sini_err         [0.04/0.26]? uncertainty on
                                                 ESP-HS projected surface
                                                 rotation frequency
 63- 82  F20.18 d-1    frot                       [0.01/0.6]? rotation frequency
                                                 measured when fitting all
                                                 observed signals
 84-104  E21.16 d-1    frot_err                   [0.0/0.03]? uncertainty on
                                                 rotation frequency measured
                                                 when fitting all observed
                                                 signals
106-124  F19.17 ---    Xc                         [0.04/0.43]? rotation
                                                 frequency relative to the
                                                 critical rotation frequency
                                                 measured when fitting all
                                                 observed signals
126-146  F21.19 ---    Xc_err                     [0.0/0.05]? uncertainty on
                                                 rotation frequency relative to
                                                 the critical rotation frequency
                                                 measured when fitting all
                                                 observed signals
    148  I1     ---    N_multiplets               [1/4]? central hydrogen mass
                                                 fraction
150-167  F18.16 d-1    f_zonal_multiplet1         [3.88/9.14]? uncertainty on
                                                 central hydrogen mass fraction
169-170  I2     ---    n_multiplet1               [-2/3]? number of observed
                                                 rotationally split multiplets
    172  I1     ---    l_multiplet1               [1/2]? zonal frequency of
                                                 observed identified multiplet
174-191  F18.16 d-1    frot_multiplet1            [0.01/0.6]? radial order of
                                                 observed identified multiplet
193-213  E21.16 d-1    frot_err_multiplet1        [0.0/0.03]? degree of observed
                                                 identified multiplet
215-233  F19.17 ---    r_rel_KnlMode_multiplet1   [0.09/0.98]? rotation rate
                                                 optimised for individual
                                                 observed identified
                                                 multipletuncertainty on
                                                 rotation rate optimised for
                                                 individual observed identified
                                                 multiplet
235-253  F19.17 ---    r_rel_KnlMean_multiplet1   [0.3/0.8]? radius relative to
                                                 total radius of maximum of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
255-273  F19.17 ---    r_rel_KnlMedian_multiplet1 [0.14/0.89]? radius relative
                                                 to total radius of mean of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
275-291  F17.15 d-1    f_zonal_multiplet2         [3.5/7.45]? radius relative to
                                                 total radius of median of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
293-294  I2     ---    n_multiplet2               [-2/2]? zonal frequency of
                                                 observed identified multiplet
    296  I1     ---    l_multiplet2               [1/2]? radial order of
                                                 observed identified multiplet
298-315  F18.16 d-1    frot_multiplet2            [0.01/0.47]? degree of
                                                 observed identified multiplet
317-337  E21.16 d-1    frot_err_multiplet2        [0.0/0.04]? rotation rate
                                                 optimised for individual
                                                 observed identified
                                                 multipletuncertainty on
                                                 rotation rate optimised for
                                                 individual observed identified
                                                 multiplet
339-357  F19.17 ---    r_rel_KnlMode_multiplet2   [0.11/0.91]? radius relative
                                                 to total radius of maximum of
                                                 the sensitivity kernel of the
                                                 observed identified multiplet
359-377  F19.17 ---    r_rel_KnlMean_multiplet2   [0.2/0.75]? radius relative to
                                                 total radius of mean of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
379-397  F19.17 ---    r_rel_KnlMedian_multiplet2 [0.12/0.86]? radius relative
                                                 to total radius of median of
                                                 the sensitivity kernel of the
                                                 observed identified multiplet
399-415  F17.15 d-1    f_zonal_multiplet3         ? zonal frequency of observed
                                                 identified multiplet
417-418  I2     ---    n_multiplet3               ? radial order of observed
                                                 identified multiplet
420-421  I2     ---    l_multiplet3               ? degree of observed
                                                 identified multiplet
423-440  F18.16 d-1    frot_multiplet3            ? rotation rate optimised for
                                                 individual observed identified
                                                 multipletuncertainty on
                                                 rotation rate optimised for
                                                 individual observed identified
                                                 multiplet
442-463  E22.17 d-1    frot_err_multiplet3        ? radius relative to total
                                                 radius of maximum of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
465-483  F19.17 ---    r_rel_KnlMode_multiplet3   ? radius relative to total
                                                 radius of mean of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
485-502  F18.16 ---    r_rel_KnlMean_multiplet3   ? radius relative to total
                                                 radius of median of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
504-521  F18.16 ---    r_rel_KnlMedian_multiplet3 ? zonal frequency of observed
                                                 identified multiplet
523-530  F8.6   d-1    f_zonal_multiplet4         ? radial order of observed
                                                 identified multiplet
532-533  I2     ---    n_multiplet4               ? degree of observed
                                                 identified multiplet
535-536  I2     ---    l_multiplet4               ? rotation rate optimised for
                                                 individual observed identified
                                                 multipletuncertainty on
                                                 rotation rate optimised for
                                                 individual observed identified
                                                 multiplet
538-554  F17.15 d-1    frot_multiplet4            ? radius relative to total
                                                 radius of maximum of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
556-575  F20.18 d-1    frot_err_multiplet4        ? radius relative to total
                                                 radius of mean of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
577-594  F18.16 ---    r_rel_KnlMode_multiplet4   ? radius relative to total
                                                 radius of median of the
                                                 sensitivity kernel of the
                                                 observed identified multiplet
--------------------------------------------------------------------------------

================================================================================
(End)          Name [Institution]                                    01-Jan-2003
