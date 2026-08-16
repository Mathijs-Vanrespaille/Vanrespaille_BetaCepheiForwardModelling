J/A+A/vol/page Asteroseismic modelling of 36 {beta} Cep stars (Vanrespaille+, 2026)
================================================================================
Asteroseismic forward modelling of 36 {beta} Cep pulsators and inferences on
their internal differential rotation
     Vanrespaille M., Fritzewski D.J., Vanlaer V., Aerts C.
    <Astron. Astrophys. vol, page (year)>
    =2003A&A...VVV.ppppI     (SIMBAD/NED BibCode)
================================================================================
Keywords_ADC: Stars, ages ; Stars, OB ; Stars, variable ; Stars, normal
Keywords: Asteroseismology – Stars: oscillations – Stars: massive –
          Stars: interiors – Stars: evolution – Stars: rotation

Abstract:
    Context. Asteroseismic observations of the interior rotation of main
    sequence stars have shown that angular momentum transport is much more
    efficient than expected. Which transport mechanisms are responsible for this
    is still unclear. Detections of radial differential rotation provide
    valuable constraints on these transport mechanisms. Differential rotation
    has been detected in several massive main sequence {beta} Cephei
    ({beta} Cep) pulsators, even though fewer than ten {beta} Cep stars have
    been asteroseismically modelled in detail so far. This suggests that
    differential rotation may be common in high-mass main sequence stars.

    Aims. We aim to expand the sample of asteroseismically forward modelled
    {beta} Cep pulsators and exploit their potential to constrain angular
    momentum transport mechanisms. To that end, we seek to place constraints on
    their rotation profiles.

    Methods. We searched for rotational splitting of non-radial modes in a large
    {beta} Cep sample with identified mode degrees. These were subjected to a
    novel forward modelling approach involving a 6-dimensional parameter space,
    which consistently accounts for second-order rotation effects using the
    state-of-the-art StORM oscillation code.

    Results. We successfully modelled 36 {beta} Cep stars and constrained
    crucial parameters such as their initial mass, internal rotation frequency,
    convective core mass, and age. Like in intermediate-mass stars, the internal
    rotation rate globally decreases in {beta} Cep stars as they evolve along
    the main sequence. Radial differential rotation is constrained in 17 {beta}
    Cep stars. We detect deviations from quasi-rigid rotation in 13 stars. Our
    constraints on eight of these 13 stars indicate the rotation rate decreases
    with radius while it increases in three other stars, and two display
    non-monotonic rotation profiles.

    Conclusions. We affirm that radial differential rotation is common in {beta}
    Cep stars. Moreover, our constrained rotation profiles suggest that the
    typical {beta} Cep rotation profile may be non-monotonic.

Description:
    ObsInput.dat: Observational properties of beta Cephei stars observed with
    Gaia and TESS. It includes the effective temperature, luminosity, the
    identified signals including their frequency and its uncertainty, amplitude,
    radial order, degree, and azimuthal order.

    Models1.dat: The statistical parameter estimates for each star with a
    satisfactory fit, including their initial mass, envelope mixing strength,
    core overshoot parameter, central hydrogen mass fraction, age, optimal
    rotation frequency, radius, surface gravity, effective temperature,
    luminosity, and convective core mass. Herein, we used the of oscillation
    computations of StORM and included the rotational splitting in our modelling
    step 3. These are the results presented in the main text of the paper.

    Models2.dat: The chosen parameters of the statistical models of each star
    with a satisfactory fit, including their initial mass, envelope mixing
    strength, core overshoot parameter, central hydrogen mass fraction, age,
    optimal rotation frequency, radius, surface gravity, effective
    temperature, luminosity, and convective core mass. Herein, we used the
    alternative method using the oscillation computations of StORM, but
    excluding the rotational splitting in modelling step 3.

    Models3.dat: The chosen parameters of the statistical models of each star
    with a satisfactory fit, including their initial mass, envelope mixing
    strength, core overshoot parameter, central hydrogen mass fraction, age,
    optimal rotation frequency, radius, surface gravity, effective
    temperature, luminosity, and convective core mass. Herein, we used the
    alternative method using the oscillation computations of GYRE,
    including the rotational splitting in modelling step 3.

    Models4.dat: The chosen parameters of the statistical models of each star
    with a satisfactory fit, including their initial mass, envelope mixing
    strength, core overshoot parameter, central hydrogen mass fraction, age,
    optimal rotation frequency, radius, surface gravity, effective
    temperature, luminosity, and convective core mass. Herein, we used the
    alternative method using the oscillation computations of GYRE, and
    excluding the rotational splitting in modelling step 3.

    RotSplit.dat: The observed and best model's rotational splitting in each
    identified multiplet used in the asteroseismic modelling. It includes the
    rotational splitting of each observed part of each multiplet as well as
    the best model's counterpart. If sufficient modes are observed, the
    asymmetry is also included. Note that there is one line per mode or
    multiplet with a unique (n,l), not one line per star.

    DiffRot.dat: The measurements of internal and surface rotation used to
    study differential rotation in Section 7 of the accompanying paper. It
    includes the Gaia surface rotation frequencies, if available, and the
    optimal rotation associated with each rotationally split multiplet
    and the relative radius where each multiplet is most sensitive.

File Summary:
--------------------------------------------------------------------------------
  FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe.txt      80        .   This file
ObsInput.dat 10294       43   The 43 stars considered and their observational
                               input
Models1.dat    646       36   Essential parameters of asteroseismic forward
                               models of 36 stars with a satisfactory fit
                               from StORM including rotational splitting in
                               the optimisation
Models2.dat    646       36   Essential parameters of asteroseismic
                               models of 36 stars with a satisfactory fit
                               from StORM excluding rotational splitting in
                               the optimisation
Models3.dat    641       36   Essential parameters of asteroseismic
                               models of 36 stars with a satisfactory fit
                               from GYRE including rotational splitting in
                               the optimisation
Models4.dat    645       36   Essential parameters of asteroseismic
                               models of 36 stars with a satisfactory fit
                               from GYRE including rotational splitting in
                               the optimisation
RotSplit.dat   451       81   Observed and modelled asymmetric rotational
                               splitting of each multiplet
DiffRot.dat    688       36   Differential rotation constraints
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
  145-  165  F21.19 K      e_log_Teff       [0.0/0.06]? uncertainty on log10 of
                                           effective temperature
  167-  182  F16.14 solLum log_L            [3.72/5.25]? log10 luminosity
  184-  203  F20.18 solLum e_log_L          [0.01/0.26]? uncertainty on log10
                                           luminosity
  205-  222  F18.16 d-1    Gaia_frot_sini   ? ESP-HS projected surface rotation
                                           frequency
  224-  241  F18.16 d-1    e_Gaia_frot_sini [0.04/0.26]? uncertainty on ESP-HS
                                           projected surface rotation frequency
        243  I1     ---    Ni_ID            [2/4]? number of identified radial
                                           modes and multiplets
  245-  246  I2     ---    Nm_ID            [3/15]? number of identified
                                           frequencies
        248  I1     ---    i1               [1]? index of multiplet this
                                           identified signal is part of
  250-  251  I2     ---    n1               [-1/3]? radial order of identified
                                           signal
        253  I1     ---    l1               ? degree of identified signal
  255-  256  I2     ---    m1               ? azimuthal order of identified
                                           signal
  258-  275  F18.16 d-1    f1               [3.99/9.74]? frequency of identified
                                           signal
  277-  281  F5.3   d-1    e_f1             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  283-  300  F18.16 ---    a1               [0.0/2.6]? amplitude of identified
                                           signal
        302  I1     ---    i2               [1/2]? index of multiplet this
                                           identified signal is part of
  304-  305  I2     ---    n2               [-2/3]? radial order of identified
                                           signal
        307  I1     ---    l2               [1/2]? degree of identified signal
  309-  310  I2     ---    m2               [-1/1]? azimuthal order of
                                           identified signal
  312-  330  F19.16 d-1    f2               [3.88/10.37]? frequency of
                                           identified signal
  332-  336  F5.3   d-1    e_f2             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  338-  355  F18.16 ---    a2               [0.0/0.1]? amplitude of identified
                                           signal
        357  I1     ---    i3               [1/2]? index of multiplet this
                                           identified signal is part of
  359-  360  I2     ---    n3               [-2/3]? radial order of identified
                                           signal
        362  I1     ---    l3               ? degree of identified signal
  364-  365  I2     ---    m3               [-1/2]? azimuthal order of
                                           identified signal
  367-  385  F19.16 d-1    f3               [3.97/11.04]? frequency of
                                           identified signal
  387-  391  F5.3   d-1    e_f3             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  393-  410  F18.16 ---    a3               [0.0/2.7]? amplitude of identified
                                           signal
        412  I1     ---    i4               [1/3]? index of multiplet this
                                           identified signal is part of
  414-  415  I2     ---    n4               [-2/3]? radial order of identified
                                           signal
        417  I1     ---    l4               ? degree of identified signal
  419-  420  I2     ---    m4               [-2/2]? azimuthal order of
                                           identified signal
  422-  439  F18.16 d-1    f4               [3.8/7.98]? frequency of identified
                                           signal
  441-  445  F5.3   d-1    e_f4             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  447-  464  F18.16 ---    a4               [0.0/0.03]? amplitude of identified
                                           signal
        466  I1     ---    i5               [1/3]? index of multiplet this
                                           identified signal is part of
  468-  469  I2     ---    n5               [-2/2]? radial order of identified
                                           signal
        471  I1     ---    l5               ? degree of identified signal
  473-  474  I2     ---    m5               [-2/2]? azimuthal order of
                                           identified signal
  476-  493  F18.16 d-1    f5               [3.25/7.65]? frequency of identified
                                           signal
  495-  499  F5.3   d-1    e_f5             [0.0/0.01]? uncertainty on frequency
                                           of identified signal
  501-  521  E21.16 ---    a5               [0.0/0.02]? amplitude of identified
                                           signal
  523-  524  I2     ---    i6               ? index of multiplet this identified
                                           signal is part of
  526-  527  I2     ---    n6               ? radial order of identified signal
  529-  530  I2     ---    l6               ? degree of identified signal
  532-  533  I2     ---    m6               ? azimuthal order of identified
                                           signal
  535-  551  F17.15 d-1    f6               ? frequency of identified signal
  553-  557  F5.3   d-1    e_f6             ? uncertainty on frequency of
                                           identified signal
  559-  572  E14.9  ---    a6               ? amplitude of identified signal
  574-  575  I2     ---    i7               ? index of multiplet this identified
                                           signal is part of
  577-  578  I2     ---    n7               ? radial order of identified signal
  580-  581  I2     ---    l7               ? degree of identified signal
  583-  584  I2     ---    m7               ? azimuthal order of identified
                                           signal
  586-  603  F18.16 d-1    f7               ? frequency of identified signal
  605-  609  F5.3   d-1    e_f7             ? uncertainty on frequency of
                                           identified signal
  611-  631  E21.16 ---    a7               ? amplitude of identified signal
  633-  634  I2     ---    i8               ? index of multiplet this identified
                                           signal is part of
  636-  637  I2     ---    n8               ? radial order of identified signal
  639-  640  I2     ---    l8               ? degree of identified signal
  642-  643  I2     ---    m8               ? azimuthal order of identified
                                           signal
  645-  661  F17.15 d-1    f8               ? frequency of identified signal
  663-  667  F5.3   d-1    e_f8             ? uncertainty on frequency of
                                           identified signal
  669-  682  E14.9  ---    a8               ? amplitude of identified signal
  684-  685  I2     ---    i9               ? index of multiplet this identified
                                           signal is part of
  687-  688  I2     ---    n9               ? radial order of identified signal
  690-  691  I2     ---    l9               ? degree of identified signal
  693-  694  I2     ---    m9               ? azimuthal order of identified
                                           signal
  696-  713  F18.16 d-1    f9               ? frequency of identified signal
  715-  719  F5.3   d-1    e_f9             ? uncertainty on frequency of
                                           identified signal
  721-  734  E14.9  ---    a9               ? amplitude of identified signal
  736-  737  I2     ---    i10              ? index of multiplet this identified
                                           signal is part of
  739-  740  I2     ---    n10              ? radial order of identified signal
  742-  743  I2     ---    l10              ? degree of identified signal
  745-  746  I2     ---    m10              ? azimuthal order of identified
                                           signal
  748-  764  F17.15 d-1    f10              ? frequency of identified signal
  766-  770  F5.3   d-1    e_f10            ? uncertainty on frequency of
                                           identified signal
  772-  793  E22.17 ---    a10              ? amplitude of identified signal
  795-  796  I2     ---    i11              ? index of multiplet this identified
                                           signal is part of
  798-  799  I2     ---    n11              ? radial order of identified signal
  801-  802  I2     ---    l11              ? degree of identified signal
  804-  805  I2     ---    m11              ? azimuthal order of identified
                                           signal
  807-  823  F17.15 d-1    f11              ? frequency of identified signal
  825-  829  F5.3   d-1    e_f11            ? uncertainty on frequency of
                                           identified signal
  831-  848  F18.16 ---    a11              ? amplitude of identified signal
  850-  851  I2     ---    i12              ? index of multiplet this identified
                                           signal is part of
  853-  854  I2     ---    n12              ? radial order of identified signal
  856-  857  I2     ---    l12              ? degree of identified signal
  859-  860  I2     ---    m12              ? azimuthal order of identified
                                           signal
  862-  878  F17.15 d-1    f12              ? frequency of identified signal
  880-  884  F5.3   d-1    e_f12            ? uncertainty on frequency of
                                           identified signal
  886-  903  F18.16 ---    a12              ? amplitude of identified signal
  905-  906  I2     ---    i13              ? index of multiplet this identified
                                           signal is part of
  908-  909  I2     ---    n13              ? radial order of identified signal
  911-  912  I2     ---    l13              ? degree of identified signal
  914-  915  I2     ---    m13              ? azimuthal order of identified
                                           signal
  917-  934  F18.16 d-1    f13              ? frequency of identified signal
  936-  940  F5.3   d-1    e_f13            ? uncertainty on frequency of
                                           identified signal
  942-  954  E13.8  ---    a13              ? amplitude of identified signal
  956-  957  I2     ---    i14              ? index of multiplet this identified
                                           signal is part of
  959-  960  I2     ---    n14              ? radial order of identified signal
  962-  963  I2     ---    l14              ? degree of identified signal
  965-  966  I2     ---    m14              ? azimuthal order of identified
                                           signal
  968-  975  F8.6   d-1    f14              ? frequency of identified signal
  977-  981  F5.3   d-1    e_f14            ? uncertainty on frequency of
                                           identified signal
  983-  996  F14.12 ---    a14              ? amplitude of identified signal
  998-  999  I2     ---    i15              ? index of multiplet this identified
                                           signal is part of
 1001- 1002  I2     ---    n15              ? radial order of identified signal
 1004- 1005  I2     ---    l15              ? degree of identified signal
 1007- 1008  I2     ---    m15              ? azimuthal order of identified
                                           signal
 1010- 1017  F8.6   d-1    f15              ? frequency of identified signal
 1019- 1023  F5.3   d-1    e_f15            ? uncertainty on frequency of
                                           identified signal
 1025- 1038  F14.12 ---    a15              ? amplitude of identified signal
 1040- 5279  A4240  d-1    all_f            ? all detected frequencies
 5281-10295  A5015  ---    all_a            ? all detected amplitudes
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models1.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9        ---    TIC_ID           ? TESS object identifier
 11- 23  A13       ---    name             ? common name
     25  I1        ---    validation_star  ? whether this is one of the six
                                     validation stars
 27- 48  F22.16    ---    chi2             [0.01/41956.99]? the reduced chi-
                                     square cost function
 50- 71  F22.20    d-1    max_f_zonal_diff [0.0/0.12]? the absolute value of the
                                     largest difference between an observed
                                     and the best model's zonal frequency
 73- 90  F18.15   solMass M                [7.74/29.62]? current mass
 92-110  F19.17   solMass e_M              [0.13/3.8]? uncertainty on current
                                     mass
112-129  F18.16   cm2.s-1 logD             [1.2/6.01]? log10 of mixing strength
                                     at the base of the envelope
131-149  F19.16   cm2.s-1 e_logD           [0.5/19.39]? uncertainty on log10 of
                                     mixing strength at the base of the
                                     envelope
151-171  F21.19    ---    fov              [0.0/0.04]? overshoot parameter
173-192  F20.18    ---    e_fov            [0.0/0.12]? uncertainty on overshoot
                                     parameter
194-213  F20.18    ---    Xc               [0.04/0.44]? central hydrogen mass
                                     fraction
215-235  F21.19    ---    e_Xc             [0.0/0.12]? uncertainty on central
                                     hydrogen mass fraction
237-255  F19.17    ---    X'c              [0.06/0.64]? normalised central
                                     hydrogen mass fraction
257-276  F20.18    ---    e_X'c            [0.0/0.16]? uncertainty on normalised
                                     central hydrogen mass fraction
278-297  F20.18    d-1    frot             [0.01/0.59]? rotation frequency
                                     measured when fitting all observed
                                     signals
299-320  E22.17    d-1    e_frot           [0.0/0.31]? uncertainty on rotation
                                     frequency measured when fitting all
                                     observed signals
322-344  E23.17    d-1    b_frot           [-0.07/0.11]? bias on rotation
                                     frequency measured when fitting all
                                     observed signals
346-365  F20.18    ---    frot_rel         [0.0/0.41]? rotation frequency
                                     relative to the critical rotation
                                     frequency measured when fitting all
                                     observed signals
367-388  E22.17    ---    e_frot_rel       [0.0/0.05]? uncertainty on rotation
                                     frequency relative to the critical
                                     rotation frequency measured when fitting
                                     all observed signals
390-391  F2.0      ---    b_frot_rel       ? bias on rotation frequency relative
                                     to the critical rotation frequency
                                     measured when fitting all observed
                                     signals
393-412  F20.16    km.s-1 vrot             [3.73/299.38]? surface rotation speed
414-435  F22.18    km.s-1 e_vrot           [0.01/167.72]? uncertainty on surface
                                     rotation speed
437-455  F19.16 1000000yr age              [2.8/35.54]? age since ZAMS
457-476  F20.17 1000000yr e_age            [0.12/10.33]? uncertainty on age
                                     since ZAMS
478-495  F18.16    K      log_Teff         [4.28/4.56]? log10 of effective
                                     temperature
497-517  F21.19    K      e_log_Teff       [0.0/0.06]? uncertainty of log10 of
                                     effective temperature
519-536  F18.16    solLum log_L            [3.68/5.2]? log10 of luminosity
538-557  F20.18    solLum e_log_L          [0.01/0.26]? uncertainty on log10 of
                                     luminosity
559-576  F18.16    cm.s-2 log_g            [3.57/3.95]? log10 of surface gravity
578-598  F21.19    cm.s-2 e_log_g          [0.0/0.09]? uncertainty on log10 of
                                     surface gravity
600-617  F18.16    solRad log_R            [0.69/1.07]? log10 of stellar radius
619-639  F21.19    solRad e_log_R          [0.0/0.04]? uncertainty of log10 of
                                     stellar radius
641-659  F19.16   solMass Mcc              [1.33/13.37]? convective core mass
661-680  F20.18   solMass e_Mcc            [0.01/1.83]? uncertainty on
                                     convective core mass
682-700  F19.17    ---    Mcc_over_M       [0.16/0.46]? convective core mass
                                     relative to total mass
702-722  F21.19    ---    e_Mcc_over_M     [0.0/0.16]? uncertainty on convective
                                     core mass relative to total mass
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models2.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9        ---    TIC_ID           ? TESS object identifier
 11- 23  A13       ---    name             ? common name
     25  I1        ---    validation_star  ? whether this is one of the six
                                     validation stars
 27- 47  F21.16    ---    chi2             [0.0/4631.44]? the reduced chi-square
                                     cost function
 49- 69  E21.16    d-1    max_f_zonal_diff [0.0/0.1]? the absolute value of the
                                     largest difference between an observed
                                     and the best model's zonal frequency
 71- 88  F18.15   solMass M                [7.85/27.99]? current mass
 90-108  F19.17   solMass e_M              [0.13/1.06]? uncertainty on current
                                     mass
110-127  F18.16   cm2.s-1 logD             [1.2/6.01]? log10 of mixing strength
                                     at the base of the envelope
129-146  F18.16   cm2.s-1 e_logD           [0.5/1.46]? uncertainty on log10 of
                                     mixing strength at the base of the
                                     envelope
148-167  F20.18    ---    fov              [0.0/0.04]? overshoot parameter
169-189  F21.19    ---    e_fov            [0.0/0.01]? uncertainty on overshoot
                                     parameter
191-209  F19.17    ---    Xc               [0.04/0.43]? central hydrogen mass
                                     fraction
211-231  F21.19    ---    e_Xc             [0.0/0.04]? uncertainty on central
                                     hydrogen mass fraction
233-251  F19.17    ---    X'c              [0.05/0.6]? normalised central
                                     hydrogen mass fraction
253-273  F21.19    ---    e_X'c            [0.0/0.05]? uncertainty on normalised
                                     central hydrogen mass fraction
275-294  F20.18    d-1    frot             [0.01/0.59]? rotation frequency
                                     measured when fitting all observed
                                     signals
296-316  E21.16    d-1    e_frot           [0.0/0.08]? uncertainty on rotation
                                     frequency measured when fitting all
                                     observed signals
318-340  F23.20    d-1    b_frot           [-0.61/0.14]? bias on rotation
                                     frequency measured when fitting all
                                     observed signals
342-361  F20.18    ---    frot_rel         [0.0/0.4]? rotation frequency
                                     relative to the critical rotation
                                     frequency measured when fitting all
                                     observed signals
363-384  E22.17    ---    e_frot_rel       [0.0/0.01]? uncertainty on rotation
                                     frequency relative to the critical
                                     rotation frequency measured when fitting
                                     all observed signals
386-387  F2.0      ---    b_frot_rel       ? bias on rotation frequency relative
                                     to the critical rotation frequency
                                     measured when fitting all observed
                                     signals
389-408  F20.16    km.s-1 vrot             [3.73/293.8]? surface rotation speed
410-430  F21.18    km.s-1 e_vrot           [0.01/36.5]? uncertainty on surface
                                     rotation speed
432-450  F19.16 1000000yr age              [3.62/31.36]? age since ZAMS
452-470  F19.17 1000000yr e_age            [0.22/3.08]? uncertainty on age since
                                     ZAMS
472-489  F18.16    K      log_Teff         [4.28/4.56]? log10 of effective
                                     temperature
491-511  F21.19    K      e_log_Teff       [0.0/0.02]? uncertainty of log10 of
                                     effective temperature
513-530  F18.16    solLum log_L            [3.67/5.17]? log10 of luminosity
532-552  F21.19    solLum e_log_L          [0.0/0.08]? uncertainty on log10 of
                                     luminosity
554-571  F18.16    cm.s-2 log_g            [3.57/3.95]? log10 of surface gravity
573-593  F21.19    cm.s-2 e_log_g          [0.0/0.02]? uncertainty on log10 of
                                     surface gravity
595-612  F18.16    solRad log_R            [0.69/1.07]? log10 of stellar radius
614-634  F21.19    solRad e_log_R          [0.0/0.01]? uncertainty of log10 of
                                     stellar radius
636-654  F19.16   solMass Mcc              [1.32/12.75]? convective core mass
656-676  F21.19   solMass e_Mcc            [0.0/0.77]? uncertainty on convective
                                     core mass
678-696  F19.17    ---    Mcc_over_M       [0.16/0.46]? convective core mass
                                     relative to total mass
698-718  F21.19    ---    e_Mcc_over_M     [0.0/0.02]? uncertainty on convective
                                     core mass relative to total mass
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models3.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9        ---    TIC_ID           ? TESS object identifier
 11- 23  A13       ---    name             ? common name
     25  I1        ---    validation_star  ? whether this is one of the six
                                     validation stars
 27- 48  F22.16    ---    chi2             [0.01/89151.27]? the reduced chi-
                                     square cost function
 50- 71  E22.17    d-1    max_f_zonal_diff [0.0/0.14]? the absolute value of the
                                     largest difference between an observed
                                     and the best model's zonal frequency
 73- 90  F18.15   solMass M                [7.74/28.63]? current mass
 92-111  F20.17   solMass e_M              [0.13/11.19]? uncertainty on current
                                     mass
113-130  F18.16   cm2.s-1 logD             [1.04/6.01]? log10 of mixing strength
                                     at the base of the envelope
132-150  F19.16   cm2.s-1 e_logD           [0.5/22.38]? uncertainty on log10 of
                                     mixing strength at the base of the
                                     envelope
152-172  F21.19    ---    fov              [0.0/0.04]? overshoot parameter
174-194  F21.19    ---    e_fov            [0.0/0.16]? uncertainty on overshoot
                                     parameter
196-214  F19.17    ---    Xc               [0.03/0.4]? central hydrogen mass
                                     fraction
216-236  F21.19    ---    e_Xc             [0.0/0.58]? uncertainty on central
                                     hydrogen mass fraction
238-256  F19.17    ---    X'c              [0.04/0.56]? normalised central
                                     hydrogen mass fraction
258-278  F21.19    ---    e_X'c            [0.0/0.86]? uncertainty on normalised
                                     central hydrogen mass fraction
280-299  F20.18    d-1    frot             [0.01/0.74]? rotation frequency
                                     measured when fitting all observed
                                     signals
301-322  E22.17    d-1    e_frot           [0.0/0.24]? uncertainty on rotation
                                     frequency measured when fitting all
                                     observed signals
324-326  F3.1      d-1    b_frot           ? bias on rotation frequency measured
                                     when fitting all observed signals
328-347  F20.18    ---    frot_rel         [0.0/0.59]? rotation frequency
                                     relative to the critical rotation
                                     frequency measured when fitting all
                                     observed signals
349-370  E22.17    ---    e_frot_rel       [0.0/0.03]? uncertainty on rotation
                                     frequency relative to the critical
                                     rotation frequency measured when fitting
                                     all observed signals
372-374  F3.1      ---    b_frot_rel       ? bias on rotation frequency relative
                                     to the critical rotation frequency
                                     measured when fitting all observed
                                     signals
376-395  F20.16    km.s-1 vrot             [3.73/351.71]? surface rotation speed
397-416  F20.17    km.s-1 e_vrot           [0.01/85.44]? uncertainty on surface
                                     rotation speed
418-436  F19.16 1000000yr age              [3.98/35.79]? age since ZAMS
438-457  F20.17 1000000yr e_age            [0.09/11.74]? uncertainty on age
                                     since ZAMS
459-476  F18.16    K      log_Teff         [4.28/4.56]? log10 of effective
                                     temperature
478-498  F21.19    K      e_log_Teff       [0.0/0.17]? uncertainty of log10 of
                                     effective temperature
500-517  F18.16    solLum log_L            [3.68/5.21]? log10 of luminosity
519-538  F20.18    solLum e_log_L          [0.0/0.82]? uncertainty on log10 of
                                     luminosity
540-557  F18.16    cm.s-2 log_g            [3.56/3.95]? log10 of surface gravity
559-579  F21.19    cm.s-2 e_log_g          [0.0/0.19]? uncertainty on log10 of
                                     surface gravity
581-598  F18.16    solRad log_R            [0.69/1.03]? log10 of stellar radius
600-620  F21.19    solRad e_log_R          [0.0/0.09]? uncertainty of log10 of
                                     stellar radius
622-640  F19.16   solMass Mcc              [1.32/12.36]? convective core mass
642-662  F21.18   solMass e_Mcc            [0.01/10.03]? uncertainty on
                                     convective core mass
664-682  F19.17    ---    Mcc_over_M       [0.16/0.44]? convective core mass
                                     relative to total mass
684-704  F21.19    ---    e_Mcc_over_M     [0.0/0.49]? uncertainty on convective
                                     core mass relative to total mass
--------------------------------------------------------------------------------

Byte-by-byte Description of file: Models4.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  9  I9        ---    TIC_ID           ? TESS object identifier
 11- 23  A13       ---    name             ? common name
     25  I1        ---    validation_star  ? whether this is one of the six
                                       validation stars
 27- 47  F21.16    ---    chi2             [0.01/9691.67]? the reduced chi-
                                       square cost function
 49- 70  E22.17    d-1    max_f_zonal_diff [0.0/0.14]? the absolute value of the
                                       largest difference between an observed
                                       and the best model's zonal frequency
 72- 89  F18.15   solMass M                [7.74/28.63]? current mass
 91-109  F19.17   solMass e_M              [0.13/0.71]? uncertainty on current
                                       mass
111-128  F18.16   cm2.s-1 logD             [1.22/6.0]? log10 of mixing strength
                                       at the base of the envelope
130-147  F18.16   cm2.s-1 e_logD           [0.5/4.31]? uncertainty on log10 of
                                       mixing strength at the base of the
                                       envelope
149-168  F20.18    ---    fov              [0.0/0.04]? overshoot parameter
170-190  F21.19    ---    e_fov            [0.0/0.01]? uncertainty on overshoot
                                       parameter
192-210  F19.17    ---    Xc               [0.03/0.4]? central hydrogen mass
                                       fraction
212-232  F21.19    ---    e_Xc             [0.0/0.04]? uncertainty on central
                                       hydrogen mass fraction
234-252  F19.17    ---    X'c              [0.04/0.56]? normalised central
                                       hydrogen mass fraction
254-274  F21.19    ---    e_X'c            [0.0/0.05]? uncertainty on normalised
                                       central hydrogen mass fraction
276-295  F20.18    d-1    frot             [0.01/0.74]? rotation frequency
                                       measured when fitting all observed
                                       signals
297-318  E22.17    d-1    e_frot           [0.0/0.04]? uncertainty on rotation
                                       frequency measured when fitting all
                                       observed signals
320-322  F3.1      d-1    b_frot           ? bias on rotation frequency measured
                                       when fitting all observed signals
324-343  F20.18    ---    frot_rel         [0.0/0.59]? rotation frequency
                                       relative to the critical rotation
                                       frequency measured when fitting all
                                       observed signals
345-366  E22.17    ---    e_frot_rel       [0.0/0.01]? uncertainty on rotation
                                       frequency relative to the critical
                                       rotation frequency measured when fitting
                                       all observed signals
368-370  F3.1      ---    b_frot_rel       ? bias on rotation frequency relative
                                       to the critical rotation frequency
                                       measured when fitting all observed
                                       signals
372-391  F20.16    km.s-1 vrot             [3.73/351.71]? surface rotation speed
393-413  F21.18    km.s-1 e_vrot           [0.01/13.77]? uncertainty on surface
                                       rotation speed
415-432  F18.15 1000000yr age              [3.98/35.79]? age since ZAMS
434-452  F19.17 1000000yr e_age            [0.1/3.28]? uncertainty on age since
                                       ZAMS
454-471  F18.16    K      log_Teff         [4.28/4.56]? log10 of effective
                                       temperature
473-493  F21.19    K      e_log_Teff       [0.0/0.02]? uncertainty of log10 of
                                       effective temperature
495-512  F18.16    solLum log_L            [3.68/5.21]? log10 of luminosity
514-533  F20.18    solLum e_log_L          [0.0/0.1]? uncertainty on log10 of
                                       luminosity
535-552  F18.16    cm.s-2 log_g            [3.56/3.95]? log10 of surface gravity
554-574  F21.19    cm.s-2 e_log_g          [0.0/0.03]? uncertainty on log10 of
                                       surface gravity
576-593  F18.16    solRad log_R            [0.69/1.03]? log10 of stellar radius
595-615  F21.19    solRad e_log_R          [0.0/0.01]? uncertainty of log10 of
                                       stellar radius
617-635  F19.16   solMass Mcc              [1.32/12.36]? convective core mass
637-656  F20.18   solMass e_Mcc            [0.01/0.33]? uncertainty on
                                       convective core mass
658-676  F19.17    ---    Mcc_over_M       [0.16/0.44]? convective core mass
                                       relative to total mass
678-698  F21.19    ---    e_Mcc_over_M     [0.0/0.02]? uncertainty on convective
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
                                       difference between zonal and m=-1 modes
110-129  F20.18 d-1    Delta_f_2_obs    ? observed frequency difference between
                                       zonal and m=2 modes
131-150  F20.18 d-1    Delta_f_-2_obs   ? observed frequency difference between
                                       zonal and m=-2 modes
152-174  E23.16 ---    A_1_obs          [-0.14/0.21]? observed asymmetry of m=1
                                       and m=-1 modes
176-196  F21.18 ---    A_2_obs          ? observed asymmetry of m=2 and m=-2
                                       modes
198-217  F20.18 d-1    Delta_f_1_model  [0.0/0.64]? modelled frequency
                                       difference between zonal and m=1 modes
219-238  F20.18 d-1    Delta_f_-1_model [0.0/0.83]? modelled frequency
                                       difference between zonal and m=-1 modes
240-259  F20.18 d-1    Delta_f_2_model  ? modelled frequency difference between
                                       zonal and m=2 modes
261-280  F20.18 d-1    Delta_f_-2_model ? observed frequency difference between
                                       zonal and m=-2 modes
282-303  F22.19 ---    A_1_model        [-0.2/0.51]? modelled asymmetry of m=1
                                       and m=-1 modes
305-325  F21.19 ---    A_2_model        ? modelled asymmetry of m=2 and m=-2
                                       modes
327-344  F18.16 d-1    frot_StORM       [0.01/0.59]? rotation frequency
                                       optimised with StORM including second
                                       order effects
346-367  E22.17 d-1    e_frot_StORM     [0.0/0.31]? uncertainty on rotation
                                       frequency optimised with StORM including
                                       second order effects
369-391  E23.17 d-1    e_frot_StORM     [-0.14/0.22]? bias on rotation frequency
                                       optimised with StORM including second
                                       order effects
393-412  F20.18 d-1    frot_GYRE        [0.01/0.77]? rotation frequency estimate
                                       from first order post-processing step
                                       using GYRE
414-433  F20.18 ---    C_nl             ? Ledoux constant
435-453  F19.17 ---    beta_nl          [0.49/1.0]? 1 - Ledoux constant
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
 63- 82  F20.18 d-1    frot                       [0.01/0.59]? rotation
                                                 frequency measured when fitting
                                                 all observed signals
 84-105  E22.17 d-1    e_l_frot                   [0.0/0.16]? lower uncertainty
                                                 on rotation frequency measured
                                                 when fitting all observed
                                                 signals
107-129  E23.17 d-1    e_u_frot                   [-0.07/0.11]? upper
                                                 uncertainty on rotation
                                                 frequency measured when fitting
                                                 all observed signals
131-149  F19.17 ---    X'c                        [0.06/0.64]? rotation
                                                 frequency relative to the
                                                 critical rotation frequency
                                                 measured when fitting all
                                                 observed signals
151-171  F21.19 ---    e_X'c                      [0.0/0.21]? uncertainty on
                                                 rotation frequency relative to
                                                 the critical rotation frequency
                                                 measured when fitting all
                                                 observed signals
    173  I1     ---    N_multiplets               [1/4]? normalised central
                                                 hydrogen mass fraction
175-192  F18.16 d-1    f_zonal_multiplet1         [3.88/9.14]? uncertainty on
                                                 normalised central hydrogen
                                                 mass fraction
194-195  I2     ---    n_multiplet1               [-2/3]? number of observed
                                                 rotationally split multiplets
    197  I1     ---    l_multiplet1               [1/2]? zonal frequency of the
                                                 first observed identified
                                                 multiplet
199-216  F18.16 d-1    frot_multiplet1            [0.01/0.59]? radial order of
                                                 the first observed identified
                                                 multiplet
218-238  E21.16 d-1    e_frot_multiplet1          [0.0/0.12]? degree of the
                                                 first observed identified
                                                 multiplet
240-261  E22.17 d-1    b_frot_multiplet1          [-0.14/0.22]? rotation rate
                                                 optimised for the first
                                                 observed identified
                                                 multipletlower uncertainty on
                                                 rotation rate optimised for the
                                                 first observed identified
                                                 multiplet
263-280  F18.16 ---    r_rel_KnlMode_multiplet1   [0.09/0.98]? upper uncertainty
                                                 on rotation rate optimised for
                                                 the first observed identified
                                                 multiplet
282-298  F17.15 ---    r_rel_KnlMean_multiplet1   [0.3/0.77]? radius relative to
                                                 total radius of maximum of the
                                                 sensitivity kernel of the first
                                                 observed identified multiplet
300-316  F17.15 ---    r_rel_KnlMedian_multiplet1 [0.14/0.86]? radius relative
                                                 to total radius of mean of
                                                 sensitivity kernel of the first
                                                 observed identified multiplet
318-334  F17.15 d-1    f_zonal_multiplet2         [3.5/7.45]? radius relative to
                                                 total radius of median of the
                                                 sensitivity kernel of the first
                                                 observed identified multiplet
336-337  I2     ---    n_multiplet2               [-2/2]? zonal frequency of the
                                                 second observed identified
                                                 multiplet
    339  I1     ---    l_multiplet2               [1/2]? radial order of the
                                                 second observed identified
                                                 multiplet
341-358  F18.16 d-1    frot_multiplet2            [0.01/0.47]? degree of the
                                                 second observed identified
                                                 multiplet
360-379  E20.15 d-1    e_frot_multiplet2          [0.0/0.08]? rotation rate
                                                 optimised for the second
                                                 observed identified
                                                 multipletlower uncertainty on
                                                 rotation rate optimised for the
                                                 second observed identified
                                                 multiplet
381-403  E23.16 d-1    b_frot_multiplet2          [-0.01/0.09]? upper
                                                 uncertainty on rotation rate
                                                 optimised for the second
                                                 observed identified multiplet
405-421  F17.15 ---    r_rel_KnlMode_multiplet2   [0.11/0.91]? radius relative
                                                 to total radius of maximum of
                                                 the sensitivity kernel of the
                                                 second observed identified
                                                 multiplet
423-439  F17.15 ---    r_rel_KnlMean_multiplet2   [0.21/0.75]? radius relative
                                                 to total radius of mean of
                                                 sensitivity kernel of the
                                                 second observed identified
                                                 multiplet
441-457  F17.15 ---    r_rel_KnlMedian_multiplet2 [0.12/0.86]? radius relative
                                                 to total radius of median of
                                                 the sensitivity kernel of the
                                                 second observed identified
                                                 multiplet
459-475  F17.15 d-1    f_zonal_multiplet3         ? zonal frequency of the third
                                                 observed identified multiplet
477-478  I2     ---    n_multiplet3               ? radial order of the third
                                                 observed identified multiplet
480-481  I2     ---    l_multiplet3               ? degree of the third observed
                                                 identified multiplet
483-500  F18.16 d-1    frot_multiplet3            ? rotation rate optimised for
                                                 the third observed identified
                                                 multipletlower uncertainty on
                                                 rotation rate optimised for the
                                                 third observed identified
                                                 multiplet
502-522  F21.19 d-1    e_frot_multiplet3          ? upper uncertainty on
                                                 rotation rate optimised for the
                                                 third observed identified
                                                 multiplet
524-544  E21.16 d-1    b_frot_multiplet3          ? radius relative to total
                                                 radius of maximum of the
                                                 sensitivity kernel of the third
                                                 observed identified multiplet
546-562  F17.15 ---    r_rel_KnlMode_multiplet3   ? radius relative to total
                                                 radius of mean of sensitivity
                                                 kernel of the third observed
                                                 identified multiplet
564-580  F17.15 ---    r_rel_KnlMean_multiplet3   ? radius relative to total
                                                 radius of median of the
                                                 sensitivity kernel of the third
                                                 observed identified multiplet
582-598  F17.15 ---    r_rel_KnlMedian_multiplet3 ? zonal frequency of the
                                                 fourth observed identified
                                                 multiplet
600-607  F8.6   d-1    f_zonal_multiplet4         ? radial order of the fourth
                                                 observed identified multiplet
609-610  I2     ---    n_multiplet4               ? degree of the fourth
                                                 observed identified multiplet
612-613  I2     ---    l_multiplet4               ? rotation rate optimised for
                                                 the fourth observed identified
                                                 multipletlower uncertainty on
                                                 rotation rate optimised for the
                                                 fourth observed identified
                                                 multiplet
615-632  F18.16 d-1    frot_multiplet4            ? upper uncertainty on
                                                 rotation rate optimised for the
                                                 fourth observed identified
                                                 multiplet
634-654  F21.19 d-1    e_frot_multiplet4          ? radius relative to total
                                                 radius of maximum of the
                                                 sensitivity kernel of the
                                                 fourth observed identified
                                                 multiplet
656-677  F22.19 d-1    b_frot_multiplet4          ? radius relative to total
                                                 radius of mean of sensitivity
                                                 kernel of the fourth observed
                                                 identified multiplet
679-695  F17.15 ---    r_rel_KnlMode_multiplet4   ? radius relative to total
                                                 radius of median of the
                                                 sensitivity kernel of the
                                                 fourth observed identified
                                                 multiplet
--------------------------------------------------------------------------------

================================================================================
(End)          Name [Institution]                                    01-Jan-2003
