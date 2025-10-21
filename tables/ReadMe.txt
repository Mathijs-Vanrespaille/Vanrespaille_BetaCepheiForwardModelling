J/A+A/vol/page  Asteroseismic modelling of 36 β Cep stars  (Vanrespaille+, 2026)
================================================================================
Asteroseismic forward modelling of 36 β Cep pulsators and inferences on their
internal differential rotation
     Vanrespaille M., Fritzewski D.J., Vanlaer V., Aerts C.
    <Astron. Astrophys. vol, page (year)>
    =2003A&A...VVV.ppppI     (SIMBAD/NED BibCode)
================================================================================
Keywords_ADC: Stars, ages ; Stars, OB ; Stars, variable ; Stars, normal
Keywords: Asteroseismology – Stars: oscillations – Stars: massive – 
          Stars: interiors – Stars: evolution – Stars: rotation

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

    Models1.dat: The chosen parameters of the statistical models of each star
    with a satisfactory fit, including their initial mass, envelope mixing
    strength, core overshoot parameter, central hydrogen mass fraction, age,
    optimal rotation frequency, radius, surface gravity, effective
    temperature, luminosity, and convective core mass. Herein, we used the
    of oscillation computations of StORM and included the rotational
    splitting in the optimisation of the model. These are the results
    presented in the paper.

    Models2.dat: The chosen parameters of the statistical models of each star
    with a satisfactory fit, including their initial mass, envelope mixing
    strength, core overshoot parameter, central hydrogen mass fraction, age,
    optimal rotation frequency, radius, surface gravity, effective
    temperature, luminosity, and convective core mass. Herein, we used the
    alternative method using the oscillation computations of StORM, but
    excluding the rotational splitting in the optimisation of the model.

    Models3.dat: The chosen parameters of the statistical models of each star
    with a satisfactory fit, including their initial mass, envelope mixing
    strength, core overshoot parameter, central hydrogen mass fraction, age,
    optimal rotation frequency, radius, surface gravity, effective
    temperature, luminosity, and convective core mass. Herein, we used the
    alternative method using the oscillation computations of GYRE,
    including the rotational splitting in the optimisation of the model.

    Models4.dat: The chosen parameters of the statistical models of each star
    with a satisfactory fit, including their initial mass, envelope mixing
    strength, core overshoot parameter, central hydrogen mass fraction, age,
    optimal rotation frequency, radius, surface gravity, effective
    temperature, luminosity, and convective core mass. Herein, we used the
    alternative method using the oscillation computations of GYRE, and
    excluding the rotational splitting in the optimisation of the model.

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
ObsInput.dat 10406       43   The 43 stars considered and their observational
                               input
Models1.dat    561       36   Essential parameters of asteroseismic forward
                               models of 36 stars with a satisfactory fit
                               from StORM including rotational splitting in
                               the optimisation
Models2.dat    561       36   Essential parameters of asteroseismic
                               models of 36 stars with a satisfactory fit
                               from StORM excluding rotational splitting in
                               the optimisation
Models3.dat    561       36   Essential parameters of asteroseismic
                               models of 36 stars with a satisfactory fit
                               from GYRE including rotational splitting in
                               the optimisation
Models4.dat    561       36   Essential parameters of asteroseismic
                               models of 36 stars with a satisfactory fit
                               from GYRE including rotational splitting in
                               the optimisation
RotSplit.dat   428       81   Observed and modelled asymmetric rotational
                               splitting of each multiplet
DiffRot.dat    594       36   Differential rotation constraints
--------------------------------------------------------------------------------

Byte-by-byte Description of file: ObsInput.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
    1-    9  I9     ---    TIC_ID           ? TESS object identifier
   11-   23  A13    ---    Name             ? common name
   25-   43  I19    ---    Gaia_DR3_ID      ? Gaia DR3 identifier
   45-   63  F19.15 deg    RAdeg            ? right ascension
   65-   81  F17.13 deg    DEdeg            ? declination
   83-   99  F17.14 ---    Vmag             [1.97/13.23]? V-band apparent
                                           magnitude
        101  I1     ---    goodFit          ? Whether a satisfactory model was
                                           found
  103-  126  A24    ---    Teff_source      ? source of effective temperature
                                           measurement
  128-  143  F16.14 K      log_Teff         [4.3/4.55]? log10 of effective
                                           temperature
  145-  165  F21.19 K      e_log_Teff       [0.0/0.02]? uncertainty on log10 of
                                           effective temperature
  167-  182  F16.14 solLum log_L            [3.72/5.25]? log10 luminosity
  184-  203  F20.18 solLum e_log_L          [0.01/0.04]? uncertainty on log10
                                           luminosity
  205-  222  F18.16 mas    plx              [0.19/4.8]? Gaia parallax
  224-  243  F20.18 mas    e_plx            [0.0/0.25]? uncertainty on Gaia
                                           parallax
  245-  254  F10.6  km.s-1 vsini            ? projected surface velocity
                                           measured by Gaia's ESP-HS pipeline
  256-  265  F10.7  km.s-1 e_vsini          [3.0/44.04]? uncertainty on
                                           projected surface velocity measured
                                           by Gaia's ESP-HS pipeline
  267-  284  F18.13 pc     distance         [274.93/5158.86]? Distance from
                                           Gaia's ESP-HS or GSPPHOT_OB pipelines
  286-  295  F10.8  ---    Ag               [0.13/4.48]? Extinction measured by
                                           Gaia's ESP-HS pipeline
  297-  314  F18.15 ---    Gmag             [3.25/12.91]? Mean G-band apparent
                                           magnitude
  316-  333  F18.16 d-1    Gaia_frot_sini   ? ESP-HS projected surface rotation
                                           frequency
  335-  352  F18.16 d-1    e_Gaia_frot_sini [0.04/0.26]? uncertainty on ESP-HS
                                           projected surface rotation frequency
        354  I1     ---    Nm_IDs           [2/4]? number of identified radial
                                           modes and multiplets
  356-  357  I2     ---    Nf_IDs           [3/15]? number of identified
                                           frequencies
        359  I1     ---    i1               [1]? index of multiplet this
                                           identified signal is part of
  361-  362  I2     ---    n1               [-1/3]? radial order of identified
                                           signal
        364  I1     ---    l1               ? degree of identified signal
  366-  367  I2     ---    m1               ? azimuthal order of identified
                                           signal
  369-  386  F18.16 d-1    f1               [3.99/9.74]? frequency of identified
                                           signal
  388-  392  F5.3   d-1    e_f1             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  394-  411  F18.16 ---    a1               [0.0/2.6]? amplitude of identified
                                           signal
        413  I1     ---    i2               [1/2]? index of multiplet this
                                           identified signal is part of
  415-  416  I2     ---    n2               [-2/3]? radial order of identified
                                           signal
        418  I1     ---    l2               [1/2]? degree of identified signal
  420-  421  I2     ---    m2               [-1/1]? azimuthal order of
                                           identified signal
  423-  441  F19.16 d-1    f2               [3.88/10.37]? frequency of
                                           identified signal
  443-  447  F5.3   d-1    e_f2             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  449-  466  F18.16 ---    a2               [0.0/0.1]? amplitude of identified
                                           signal
        468  I1     ---    i3               [1/2]? index of multiplet this
                                           identified signal is part of
  470-  471  I2     ---    n3               [-2/3]? radial order of identified
                                           signal
        473  I1     ---    l3               ? degree of identified signal
  475-  476  I2     ---    m3               [-1/2]? azimuthal order of
                                           identified signal
  478-  496  F19.16 d-1    f3               [3.97/11.04]? frequency of
                                           identified signal
  498-  502  F5.3   d-1    e_f3             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  504-  521  F18.16 ---    a3               [0.0/2.7]? amplitude of identified
                                           signal
        523  I1     ---    i4               [1/3]? index of multiplet this
                                           identified signal is part of
  525-  526  I2     ---    n4               [-2/3]? radial order of identified
                                           signal
        528  I1     ---    l4               ? degree of identified signal
  530-  531  I2     ---    m4               [-2/2]? azimuthal order of
                                           identified signal
  533-  550  F18.16 d-1    f4               [3.8/7.98]? frequency of identified
                                           signal
  552-  556  F5.3   d-1    e_f4             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  558-  575  F18.16 ---    a4               [0.0/0.03]? amplitude of identified
                                           signal
        577  I1     ---    i5               [1/3]? index of multiplet this
                                           identified signal is part of
  579-  580  I2     ---    n5               [-2/2]? radial order of identified
                                           signal
        582  I1     ---    l5               ? degree of identified signal
  584-  585  I2     ---    m5               [-2/2]? azimuthal order of
                                           identified signal
  587-  604  F18.16 d-1    f5               [3.25/7.65]? frequency of identified
                                           signal
  606-  610  F5.3   d-1    e_f5             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  612-  632  E21.16 ---    a5               [0.0/0.02]? amplitude of identified
                                           signal
  634-  635  I2     ---    i6               ? index of multiplet this identified
                                           signal is part of
  637-  638  I2     ---    n6               ? radial order of identified signal
  640-  641  I2     ---    l6               ? degree of identified signal
  643-  644  I2     ---    m6               ? azimuthal order of identified
                                           signal
  646-  662  F17.15 d-1    f6               ? frequency of identified signal
  664-  668  F5.3   d-1    e_f6             ? uncertainty on frequency of
                                           identified signal
  670-  683  E14.9  ---    a6               ? amplitude of identified signal
  685-  686  I2     ---    i7               ? index of multiplet this identified
                                           signal is part of
  688-  689  I2     ---    n7               ? radial order of identified signal
  691-  692  I2     ---    l7               ? degree of identified signal
  694-  695  I2     ---    m7               ? azimuthal order of identified
                                           signal
  697-  714  F18.16 d-1    f7               ? frequency of identified signal
  716-  720  F5.3   d-1    e_f7             ? uncertainty on frequency of
                                           identified signal
  722-  742  E21.16 ---    a7               ? amplitude of identified signal
  744-  745  I2     ---    i8               ? index of multiplet this identified
                                           signal is part of
  747-  748  I2     ---    n8               ? radial order of identified signal
  750-  751  I2     ---    l8               ? degree of identified signal
  753-  754  I2     ---    m8               ? azimuthal order of identified
                                           signal
  756-  772  F17.15 d-1    f8               ? frequency of identified signal
  774-  778  F5.3   d-1    e_f8             ? uncertainty on frequency of
                                           identified signal
  780-  793  E14.9  ---    a8               ? amplitude of identified signal
  795-  796  I2     ---    i9               ? index of multiplet this identified
                                           signal is part of
  798-  799  I2     ---    n9               ? radial order of identified signal
  801-  802  I2     ---    l9               ? degree of identified signal
  804-  805  I2     ---    m9               ? azimuthal order of identified
                                           signal
  807-  824  F18.16 d-1    f9               ? frequency of identified signal
  826-  830  F5.3   d-1    e_f9             ? uncertainty on frequency of
                                           identified signal
  832-  845  E14.9  ---    a9               ? amplitude of identified signal
  847-  848  I2     ---    i10              ? index of multiplet this identified
                                           signal is part of
  850-  851  I2     ---    n10              ? radial order of identified signal
  853-  854  I2     ---    l10              ? degree of identified signal
  856-  857  I2     ---    m10              ? azimuthal order of identified
                                           signal
  859-  875  F17.15 d-1    f10              ? frequency of identified signal
  877-  881  F5.3   d-1    e_f10            ? uncertainty on frequency of
                                           identified signal
  883-  904  E22.17 ---    a10              ? amplitude of identified signal
  906-  907  I2     ---    i11              ? index of multiplet this identified
                                           signal is part of
  909-  910  I2     ---    n11              ? radial order of identified signal
  912-  913  I2     ---    l11              ? degree of identified signal
  915-  916  I2     ---    m11              ? azimuthal order of identified
                                           signal
  918-  934  F17.15 d-1    f11              ? frequency of identified signal
  936-  940  F5.3   d-1    e_f11            ? uncertainty on frequency of
                                           identified signal
  942-  959  F18.16 ---    a11              ? amplitude of identified signal
  961-  962  I2     ---    i12              ? index of multiplet this identified
                                           signal is part of
  964-  965  I2     ---    n12              ? radial order of identified signal
  967-  968  I2     ---    l12              ? degree of identified signal
  970-  971  I2     ---    m12              ? azimuthal order of identified
                                           signal
  973-  989  F17.15 d-1    f12              ? frequency of identified signal
  991-  995  F5.3   d-1    e_f12            ? uncertainty on frequency of
                                           identified signal
  997- 1014  F18.16 ---    a12              ? amplitude of identified signal
 1016- 1017  I2     ---    i13              ? index of multiplet this identified
                                           signal is part of
 1019- 1020  I2     ---    n13              ? radial order of identified signal
 1022- 1023  I2     ---    l13              ? degree of identified signal
 1025- 1026  I2     ---    m13              ? azimuthal order of identified
                                           signal
 1028- 1045  F18.16 d-1    f13              ? frequency of identified signal
 1047- 1051  F5.3   d-1    e_f13            ? uncertainty on frequency of
                                           identified signal
 1053- 1065  E13.8  ---    a13              ? amplitude of identified signal
 1067- 1068  I2     ---    i14              ? index of multiplet this identified
                                           signal is part of
 1070- 1071  I2     ---    n14              ? radial order of identified signal
 1073- 1074  I2     ---    l14              ? degree of identified signal
 1076- 1077  I2     ---    m14              ? azimuthal order of identified
                                           signal
 1079- 1086  F8.6   d-1    f14              ? frequency of identified signal
 1088- 1092  F5.3   d-1    e_f14            ? uncertainty on frequency of
                                           identified signal
 1094- 1107  F14.12 ---    a14              ? amplitude of identified signal
 1109- 1110  I2     ---    i15              ? index of multiplet this identified
                                           signal is part of
 1112- 1113  I2     ---    n15              ? radial order of identified signal
 1115- 1116  I2     ---    l15              ? degree of identified signal
 1118- 1119  I2     ---    m15              ? azimuthal order of identified
                                           signal
 1121- 1128  F8.6   d-1    f15              ? frequency of identified signal
 1130- 1134  F5.3   d-1    e_f15            ? uncertainty on frequency of
                                           identified signal
 1136- 1149  F14.12 ---    a15              ? amplitude of identified signal
 1151- 5390  A4240  d-1    all_f            ? all detected frequencies
 5392-10406  A5015  ---    all_a            ? all detected amplitudes
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models1.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9      ---    TIC_ID           ? TESS object identifier
 11- 23  A13     ---    name             ? common name
     25  I1      ---    validation_star  ? whether this is one of the six
                                       validation stars
 27- 48  F22.16  ---    chi2             [0.78/94808.76]? the chi-square cost
                                       function
 50- 71  F22.16  ---    chi2_reduced     [0.78/94808.76]? the reduced chi-square
                                       cost function
 73- 93  E21.16  d-1    max_f_zonal_diff [0.0/0.11]? the absolute value of the
                                       largest difference between an observed
                                       and the best model's zonal frequency
 95-112  F18.15 solMass M                [8.03/29.86]? initial mass
114-133  F20.18 solMass e_M              [0.05/1.04]? uncertainty on initial
                                       mass
135-152  F18.16 cm2.s-1 logD             [1.0/6.0]? log10 of mixing strength at
                                       the base of the envelope
154-171  F18.16 cm2.s-1 e_logD           [0.5/1.51]? uncertainty on log10 of
                                       mixing strength at the base of the
                                       envelope
173-193  F21.19  ---    fov              [0.0/0.04]? overshoot parameter
195-215  F21.19  ---    e_fov            [0.0/0.01]? uncertainty on overshoot
                                       parameter
217-235  F19.17  ---    Xc               [0.04/0.43]? central hydrogen mass
                                       fraction
237-257  F21.19  ---    e_Xc             [0.0/0.05]? uncertainty on central
                                       hydrogen mass fraction
259-278  F20.18  d-1    frot             [0.01/0.6]? rotation frequency measured
                                       when fitting all observed signals
280-300  E21.16  d-1    e_frot           [0.0/0.03]? uncertainty on rotation
                                       frequency measured when fitting all
                                       observed signals
302-322  F21.19  ---    frot_rel         [0.0/0.41]? rotation frequency relative
                                       to the critical rotation frequency
                                       measured when fitting all observed
                                       signals
324-345  E22.17  ---    e_frot_rel       [0.0/0.01]? uncertainty on rotation
                                       frequency relative to the critical
                                       rotation frequency measured when fitting
                                       all observed signals
347-364  F18.9   yr     age              [3282211.88/32511686.79]? age since
                                       ZAMS
366-384  F19.11  yr     e_age            [87491.68/1499569.31]? uncertainty on
                                       age since ZAMS
386-403  F18.16  K      log_Teff         [4.28/4.57]? log10 of effective
                                       temperature
405-425  F21.19  K      e_log_Teff       [0.0/0.02]? uncertainty of log10 of
                                       effective temperature
427-444  F18.16  solLum log_L            [3.66/5.21]? log10 of luminosity
446-465  F20.18  solLum e_log_L          [0.01/0.04]? uncertainty on log10 of
                                       luminosity
467-484  F18.16  cm.s-2 log_g            [3.56/3.95]? log10 of surface gravity
486-506  F21.19  cm.s-2 e_log_g          [0.0/0.05]? uncertainty on log10 of
                                       surface gravity
508-525  F18.16  solRad log_R            [0.69/1.07]? log10 of stellar radius
527-547  F21.19  solRad e_log_R          [0.0/0.02]? uncertainty of log10 of
                                       stellar radius
549-567  F19.16 solMass Mcc              [1.33/12.96]? convective core mass
569-589  F21.19 solMass e_Mcc            [0.0/0.32]? uncertainty on convective
                                       core mass
591-609  F19.17  ---    Mcc_over_M       [0.16/0.44]? convective core mass
                                       relative to total mass
611-631  F21.19  ---    e_Mcc_over_M     [0.0/0.02]? uncertainty on convective
                                       core mass relative to total mass
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models2.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9      ---    TIC_ID           ? TESS object identifier
 11- 23  A13     ---    name             ? common name
     25  I1      ---    validation_star  ? whether this is one of the six
                                       validation stars
 27- 47  F21.16  ---    chi2             [0.04/9672.75]? the chi-square cost
                                       function
 49- 69  F21.16  ---    chi2_reduced     [0.04/9672.75]? the reduced chi-square
                                       cost function
 71- 90  E20.15  d-1    max_f_zonal_diff [0.0/0.1]? the absolute value of the
                                       largest difference between an observed
                                       and the best model's zonal frequency
 92-109  F18.15 solMass M                [7.86/27.51]? initial mass
111-130  F20.18 solMass e_M              [0.05/1.1]? uncertainty on initial mass
132-149  F18.16 cm2.s-1 logD             [1.0/6.0]? log10 of mixing strength at
                                       the base of the envelope
151-168  F18.16 cm2.s-1 e_logD           [0.5/1.51]? uncertainty on log10 of
                                       mixing strength at the base of the
                                       envelope
170-189  F20.18  ---    fov              [0.0/0.04]? overshoot parameter
191-211  F21.19  ---    e_fov            [0.0/0.02]? uncertainty on overshoot
                                       parameter
213-232  F20.18  ---    Xc               [0.04/0.42]? central hydrogen mass
                                       fraction
234-254  F21.19  ---    e_Xc             [0.0/0.05]? uncertainty on central
                                       hydrogen mass fraction
256-275  F20.18  d-1    frot             [0.01/0.59]? rotation frequency
                                       measured when fitting all observed
                                       signals
277-297  E21.16  d-1    e_frot           [0.0/0.23]? uncertainty on rotation
                                       frequency measured when fitting all
                                       observed signals
299-318  F20.18  ---    frot_rel         [0.0/0.41]? rotation frequency relative
                                       to the critical rotation frequency
                                       measured when fitting all observed
                                       signals
320-341  E22.17  ---    e_frot_rel       [0.0/0.04]? uncertainty on rotation
                                       frequency relative to the critical
                                       rotation frequency measured when fitting
                                       all observed signals
343-361  F19.10  yr     age              [3810812.81/33045606.93]? age since
                                       ZAMS
363-381  F19.11  yr     e_age            [88669.64/1653686.97]? uncertainty on
                                       age since ZAMS
383-400  F18.16  K      log_Teff         [4.27/4.56]? log10 of effective
                                       temperature
402-422  F21.19  K      e_log_Teff       [0.0/0.02]? uncertainty of log10 of
                                       effective temperature
424-441  F18.16  solLum log_L            [3.66/5.18]? log10 of luminosity
443-462  F20.18  solLum e_log_L          [0.01/0.04]? uncertainty on log10 of
                                       luminosity
464-481  F18.16  cm.s-2 log_g            [3.56/3.95]? log10 of surface gravity
483-503  F21.19  cm.s-2 e_log_g          [0.0/0.05]? uncertainty on log10 of
                                       surface gravity
505-522  F18.16  solRad log_R            [0.69/1.08]? log10 of stellar radius
524-544  F21.19  solRad e_log_R          [0.0/0.02]? uncertainty of log10 of
                                       stellar radius
546-564  F19.16 solMass Mcc              [1.32/12.24]? convective core mass
566-586  F21.19 solMass e_Mcc            [0.0/0.34]? uncertainty on convective
                                       core mass
588-606  F19.17  ---    Mcc_over_M       [0.16/0.45]? convective core mass
                                       relative to total mass
608-628  F21.19  ---    e_Mcc_over_M     [0.0/0.03]? uncertainty on convective
                                       core mass relative to total mass
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models3.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9      ---    TIC_ID           ? TESS object identifier
 11- 23  A13     ---    name             ? common name
     25  I1      ---    validation_star  ? whether this is one of the six
                                       validation stars
 27- 49  F23.16  ---    chi2             [0.78/188935.94]? the chi-square cost
                                       function
 51- 72  F22.16  ---    chi2_reduced     [0.78/94467.97]? the reduced chi-square
                                       cost function
 74- 95  E22.17  d-1    max_f_zonal_diff [0.0/0.18]? the absolute value of the
                                       largest difference between an observed
                                       and the best model's zonal frequency
 97-114  F18.15 solMass M                [8.03/25.12]? initial mass
116-135  F20.18 solMass e_M              [0.05/1.04]? uncertainty on initial
                                       mass
137-154  F18.16 cm2.s-1 logD             [1.0/6.0]? log10 of mixing strength at
                                       the base of the envelope
156-173  F18.16 cm2.s-1 e_logD           [0.5/1.51]? uncertainty on log10 of
                                       mixing strength at the base of the
                                       envelope
175-195  F21.19  ---    fov              [0.0/0.04]? overshoot parameter
197-217  F21.19  ---    e_fov            [0.0/0.01]? uncertainty on overshoot
                                       parameter
219-237  F19.17  ---    Xc               [0.04/0.4]? central hydrogen mass
                                       fraction
239-259  F21.19  ---    e_Xc             [0.0/0.02]? uncertainty on central
                                       hydrogen mass fraction
261-280  F20.18  d-1    frot             [0.01/0.74]? rotation frequency
                                       measured when fitting all observed
                                       signals
282-302  E21.16  d-1    e_frot           [0.0/0.02]? uncertainty on rotation
                                       frequency measured when fitting all
                                       observed signals
304-324  F21.19  ---    frot_rel         [0.0/0.59]? rotation frequency relative
                                       to the critical rotation frequency
                                       measured when fitting all observed
                                       signals
326-347  E22.17  ---    e_frot_rel       [0.0/0.01]? uncertainty on rotation
                                       frequency relative to the critical
                                       rotation frequency measured when fitting
                                       all observed signals
349-366  F18.9   yr     age              [4499226.1/32499990.34]? age since ZAMS
368-386  F19.11  yr     e_age            [101119.11/1499581.47]? uncertainty on
                                       age since ZAMS
388-405  F18.16  K      log_Teff         [4.27/4.54]? log10 of effective
                                       temperature
407-427  F21.19  K      e_log_Teff       [0.0/0.01]? uncertainty of log10 of
                                       effective temperature
429-446  F18.16  solLum log_L            [3.66/5.08]? log10 of luminosity
448-467  F20.18  solLum e_log_L          [0.01/0.04]? uncertainty on log10 of
                                       luminosity
469-486  F18.16  cm.s-2 log_g            [3.54/3.95]? log10 of surface gravity
488-508  F21.19  cm.s-2 e_log_g          [0.0/0.02]? uncertainty on log10 of
                                       surface gravity
510-527  F18.16  solRad log_R            [0.69/1.03]? log10 of stellar radius
529-549  F21.19  solRad e_log_R          [0.0/0.02]? uncertainty of log10 of
                                       stellar radius
551-568  F18.16 solMass Mcc              [1.32/9.73]? convective core mass
570-590  F21.19 solMass e_Mcc            [0.0/0.42]? uncertainty on convective
                                       core mass
592-610  F19.17  ---    Mcc_over_M       [0.16/0.39]? convective core mass
                                       relative to total mass
612-632  F21.19  ---    e_Mcc_over_M     [0.0/0.02]? uncertainty on convective
                                       core mass relative to total mass
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models4.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9      ---    TIC_ID           ? TESS object identifier
 11- 23  A13     ---    name             ? common name
     25  I1      ---    validation_star  ? whether this is one of the six
                                       validation stars
 27- 48  F22.16  ---    chi2             [0.31/30041.15]? the chi-square cost
                                       function
 50- 71  F22.16  ---    chi2_reduced     [0.31/30041.15]? the reduced chi-square
                                       cost function
 73- 94  E22.17  d-1    max_f_zonal_diff [0.0/0.18]? the absolute value of the
                                       largest difference between an observed
                                       and the best model's zonal frequency
 96-113  F18.15 solMass M                [8.03/25.12]? initial mass
115-134  F20.18 solMass e_M              [0.05/1.34]? uncertainty on initial
                                       mass
136-153  F18.16 cm2.s-1 logD             [1.0/6.0]? log10 of mixing strength at
                                       the base of the envelope
155-172  F18.16 cm2.s-1 e_logD           [0.5/1.51]? uncertainty on log10 of
                                       mixing strength at the base of the
                                       envelope
174-193  F20.18  ---    fov              [0.0/0.04]? overshoot parameter
195-215  F21.19  ---    e_fov            [0.0/0.02]? uncertainty on overshoot
                                       parameter
217-235  F19.17  ---    Xc               [0.04/0.4]? central hydrogen mass
                                       fraction
237-257  F21.19  ---    e_Xc             [0.0/0.03]? uncertainty on central
                                       hydrogen mass fraction
259-278  F20.18  d-1    frot             [0.01/0.74]? rotation frequency
                                       measured when fitting all observed
                                       signals
280-300  E21.16  d-1    e_frot           [0.0/0.02]? uncertainty on rotation
                                       frequency measured when fitting all
                                       observed signals
302-321  F20.18  ---    frot_rel         [0.0/0.59]? rotation frequency relative
                                       to the critical rotation frequency
                                       measured when fitting all observed
                                       signals
323-344  E22.17  ---    e_frot_rel       [0.0/0.01]? uncertainty on rotation
                                       frequency relative to the critical
                                       rotation frequency measured when fitting
                                       all observed signals
346-363  F18.9   yr     age              [4499226.1/32982479.01]? age since ZAMS
365-383  F19.11  yr     e_age            [101116.03/1500241.88]? uncertainty on
                                       age since ZAMS
385-402  F18.16  K      log_Teff         [4.27/4.54]? log10 of effective
                                       temperature
404-424  F21.19  K      e_log_Teff       [0.0/0.01]? uncertainty of log10 of
                                       effective temperature
426-443  F18.16  solLum log_L            [3.66/5.08]? log10 of luminosity
445-464  F20.18  solLum e_log_L          [0.01/0.04]? uncertainty on log10 of
                                       luminosity
466-483  F18.16  cm.s-2 log_g            [3.54/3.95]? log10 of surface gravity
485-505  F21.19  cm.s-2 e_log_g          [0.0/0.03]? uncertainty on log10 of
                                       surface gravity
507-524  F18.16  solRad log_R            [0.69/1.03]? log10 of stellar radius
526-546  F21.19  solRad e_log_R          [0.0/0.01]? uncertainty of log10 of
                                       stellar radius
548-565  F18.16 solMass Mcc              [1.32/9.73]? convective core mass
567-587  F21.19 solMass e_Mcc            [0.0/0.42]? uncertainty on convective
                                       core mass
589-607  F19.17  ---    Mcc_over_M       [0.16/0.39]? convective core mass
                                       relative to total mass
609-629  F21.19  ---    e_Mcc_over_M     [0.0/0.02]? uncertainty on convective
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
346-367  E22.17 d-1    e_frot_StORM     [0.0/0.04]? modelled asymmetry of m=2
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
 44- 61  F18.16 d-1    e_Gaia_frot_sini           [0.04/0.26]? uncertainty on
                                                 ESP-HS projected surface
                                                 rotation frequency
 63- 82  F20.18 d-1    frot                       [0.01/0.6]? rotation frequency
                                                 measured when fitting all
                                                 observed signals
 84-104  E21.16 d-1    e_frot                     [0.0/0.03]? uncertainty on
                                                 rotation frequency measured
                                                 when fitting all observed
                                                 signals
106-124  F19.17 ---    Xc                         [0.04/0.43]? rotation
                                                 frequency relative to the
                                                 critical rotation frequency
                                                 measured when fitting all
                                                 observed signals
126-146  F21.19 ---    e_Xc                       [0.0/0.05]? uncertainty on
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
