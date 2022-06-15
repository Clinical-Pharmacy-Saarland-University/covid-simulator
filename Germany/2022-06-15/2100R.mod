;; 1. Based on: 2100
;; 2. Description: hosp+d+ICU+F_LD_ICU+Vent+LD_KH_L(O); cp Inf+Hosp(BaWue); F RR_Hosp+ICU+dout; F_Hosp(BL); hill ICU+d+hosp; PQ new
;; x1. Author: user

;; 3. Label: LD_B_L(MT), hosp(MT,NT,A,D,V), dout(MT,PQ), dKH(PQ), ICU(MT,A,V,D)
;; Kommentar: Für die Simulationen zukünftiger Impfungen muss die Anzahl an Impfungen auf die Erstimpfung und mit 6 Wochen Zeitverzug auf die Zweitimpfungen addiert werden

$SIZES NO=-10000 PD=-5000 LVR=50 PC=32 PG=-500

$PROBLEM COVID-19 cases

$INPUT ID Week Year STATE=DROP DROP TIME TIME2 DV OBS CMT AMT EVID MDV NUMBER CASE0 DIVIFLAG SCHOOLS STORECLOSURE CURFEW RESTRAININGORDER CARNIVAL FLAGZEROOBS FLAGRKI FLAGMP DOSING2 MAXIMUMICU Zweit05to14 Zweit15to59 Zweit60 ZweitGesamt Booster05to14 Booster15to59 Booster60 BoosterGesamt EWZ_0 EWZ_5 EWZ_15 EWZ_60 EWZ_Gesamt T0 A00toA04_M A00toA04_W A05toA14_M A05toA14_W A15toA34_M A15toA34_W A35toA59_M A35toA59_W A60toA79_M A60toA79_W A80_M A80_W Positivenquote AnzahlTestungen PREDIC fixR0_49 fixCP_49 fixR0_48 fixCP_48 fixR0_47 fixCP_47 fixR0_46 fixCP_46 fixR0_45 fixCP_45 fixR0_44 fixCP_44 fixR0_43 fixCP_43 fixR0_0 fixR0_1 fixR0_2 fixR0_3 fixR0_4 fixR0_5 fixR0_6 fixR0_7 fixR0_8 fixR0_9 fixR0_10 fixR0_11 fixR0_12 fixR0_13 fixR0_14 fixR0_15 fixR0_16 fixR0_17 fixR0_18 fixR0_19 fixR0_20 fixR0_21 fixR0_22 fixR0_23 fixR0_24 fixR0_25 fixR0_26 fixR0_27 fixR0_28 fixR0_29 fixR0_30 fixR0_31 fixR0_32 fixR0_33 fixR0_34 fixR0_35 fixR0_36 fixR0_37 fixR0_38 fixR0_39 fixR0_40 fixR0_41 fixR0_42 fixCP_4 fixCP_5 fixCP_6 fixCP_7 fixCP_8 fixCP_9 fixCP_10 fixCP_11 fixCP_12 fixCP_13 fixCP_14 fixCP_15 fixCP_16 fixCP_17 fixCP_18 fixCP_19 fixCP_20 fixCP_21 fixCP_22 fixCP_23 fixCP_24 fixCP_25 fixCP_26 fixCP_27 fixCP_28 fixCP_29 fixCP_30 fixCP_31 fixCP_32 fixCP_33 fixCP_34 fixCP_35 fixCP_36 fixCP_37 fixCP_38 fixCP_39 fixCP_40 fixCP_41 fixCP_42 fixF2 fixR0FREE MDVLAST
$DATA ..\DATASET\COVID_ds_228_Impf_VOC_Betten_V2_R0_2001_41.csv IGNORE=@ IGNORE(ID.EQ.17) IGNORE(DOSING2.EQ.1) IGNORE(DIVIFLAG.EQ.2) IGNORE(CMT.EQ.8) IGNORE(CMT.EQ.11) IGNORE(CMT.EQ.12) IGNORE(CMT.EQ.9) IGNORE(FLAGZEROOBS.EQ.1) IGNORE(FLAGRKI.EQ.1) IGNORE(PREDIC.EQ.1)

$SUBROUTINES ADVAN6 TOL=4
$MODEL NCOMPS=31

$PK
fBER = 1
IF(ID.EQ.3) fBER = THETA(8)
fHE = 1
IF(ID.EQ.7) fHE = 1.15
fBR = 1
IF(ID.EQ.5) fBR = 1.3
fHH = 1
fHH_KH = 1
IF(ID.EQ.6) fHH = 1.2
IF(ID.EQ.6.AND.TIME.GE.569) fHH_KH = 1.2 ;ab 13.07.2021 werden in HH nur noch stationär gesamt gemeldet (Patienten mit und ohne Wohnsitz in Hamburg)

R0START = fixR0_0

COV1 = SCHOOLS
E1 = fixR0_1
COV2 = CARNIVAL
E2 = fixR0_2
COV3 = MAX(CURFEW,RESTRAININGORDER)
E3 = fixR0_3
MTDIFF = 1
MTIME(1) = fixCP_4
COV4 = MPAST(1)
E4 = fixR0_4
MTIME(2) = fixCP_5
COV5 = MPAST(2)
E5 = fixR0_5
MTIME(3) = fixCP_6
COV6 = MPAST(3)
E6 = fixR0_6
MTIME(4) = fixCP_7
COV7 = MPAST(4)
E7 = fixR0_7
MTIME(7) = fixCP_8
COV8 = MPAST(7)
E8 = fixR0_8
MTIME(8) = fixCP_9
COV9 = MPAST(8)
E9 = fixR0_9
MTIME(9) = fixCP_10
COV10 = MPAST(9)
E10 = fixR0_10
MTIME(11) = fixCP_11
COV11 = MPAST(11)
E11 =fixR0_11
MTIME(12) = fixCP_12
COV12 = MPAST(12)
E12 = fixR0_12
MTIME(13) = fixCP_13
COV13 = MPAST(13)
E13 = fixR0_13
MTIME(14) = fixCP_14
COV14 = MPAST(14)
E14 = fixR0_14
MTIME(16) = fixCP_15
COV15 = MPAST(16)
E15= fixR0_15
MTIME(17) = fixCP_16
COV16= MPAST(17)
E16= fixR0_16
MTIME(18) = fixCP_17
COV17= MPAST(18)
E17= fixR0_17
MTIME(19) = fixCP_18
COV18= MPAST(19)
E18= fixR0_18
MTIME(20) = fixCP_19
COV19= MPAST(20)
E19= fixR0_19
MTIME(21) = fixCP_20
COV20 = MPAST(21)
E20= fixR0_20
MTIME(22) = fixCP_21
COV21 = MPAST(22)
E21 = fixR0_21
MTIME(23) = fixCP_22
COV22 = MPAST(23)
E22 = fixR0_22
MTIME(24) = fixCP_23
COV23 = MPAST(24)
E23 = fixR0_23
MTIME(25) = fixCP_24
COV24 = MPAST(25)
E24 = fixR0_24
MTIME(26) = fixCP_25
COV25 = MPAST(26)
E25= fixR0_25
MTIME(26) = fixCP_25
COV25 = MPAST(26)
E25= fixR0_25
MTIME(27) = fixCP_26
COV26 = MPAST(27)
E26= fixR0_26
MTIME(28) = fixCP_27
COV27 = MPAST(28)
E27= fixR0_27
MTIME(29) = fixCP_28
COV28 = MPAST(29)
E28= fixR0_28
MTIME(30) = fixCP_29
COV29 = MPAST(30)
E29= fixR0_29
MTIME(31) = fixCP_30
COV30 = MPAST(31)
E30= fixR0_30
MTIME(32) = fixCP_31
COV31 = MPAST(32)
E31= fixR0_31
MTIME(33) = fixCP_32
COV32 = MPAST(33)
E32= fixR0_32
MTIME(34) = fixCP_33
COV33 = MPAST(34)
E33= fixR0_33
MTIME(35) = fixCP_34
COV34 = MPAST(35)
E34= fixR0_34
MTIME(36) = fixCP_35
COV35 = MPAST(36)
E35= fixR0_35
MTIME(37) = fixCP_36
COV36 = MPAST(37)
E36= fixR0_36
MTIME(38) = fixCP_37
COV37 = MPAST(38)
E37= fixR0_37
MTIME(39) = fixCP_38
COV38 = MPAST(39)
E38= fixR0_38
MTIME(40) = fixCP_39
COV39 = MPAST(40)
E39= fixR0_39
MTIME(41) = fixCP_40
COV40 = MPAST(41)
E40 = fixR0_40
MTIME(42) = fixCP_41
COV41 = MPAST(42)
E41 = fixR0_41
MTIME(43) = fixCP_42
COV42 = MPAST(43)
E42 = fixR0_42
MTIME(44) = fixCP_43
COV43 = MPAST(44)
E43 = fixR0_43
MTIME(45) = fixCP_44
COV44 = MPAST(45)
E44 = fixR0_44
MTIME(47) = fixCP_45
COV45 = MPAST(47)
E45 = fixR0_45
MTIME(48) = fixCP_46
COV46 = MPAST(48)
E46 = fixR0_46
MTIME(49) = fixCP_47
COV47 = MPAST(49)
E47 = fixR0_47
MTIME(50) = fixCP_48
COV48 = MPAST(50)
E48 = fixR0_48
MTIME(51) = fixCP_49
COV49 = MPAST(51)
E49 = fixR0_49

R0FREE_OLD = R0START*(1-COV1)+E2*COV2+E1*(COV1-COV3)+E3*(COV3-COV4)+E4*(COV4-COV5)+E5*(COV5-COV6)+E6*(COV6-COV7)+E7*(COV7-COV8)+E8*(COV8-COV9)+E9*(COV9-COV10)+E10*(COV10-COV11)+E11*(COV11-COV12)+E12*(COV12-COV13)+E13*(COV13-COV14)+E14*(COV14-COV15)+E15*(COV15-COV16)+E16*(COV16-COV17)+E17*(COV17-COV18)+E18*(COV18-COV19)+E19*(COV19-COV20)+E20*(COV20-COV21)+E21*(COV21-COV22)+E22*(COV22-COV23)+E23*(COV23-COV24)+E24*(COV24-COV25)+E25*(COV25-COV26)+E26*(COV26-COV27)+E27*(COV27-COV28)+E28*(COV28-COV29)+E29*(COV29-COV30)+E30*(COV30-COV31)+E31*(COV31-COV32)+E32*(COV32-COV33)+E33*(COV33-COV34)+E34*(COV34-COV35)+E35*(COV35-COV36)+E36*(COV36-COV37)+E37*(COV37-COV38)+E38*(COV38-COV39)+E39*(COV39-COV40)+E40*(COV40-COV41)+E41*(COV41-COV42)+E42*(COV42-COV43)+E43*(COV43-COV44)+E44*(COV44-COV45)+E45*(COV45-COV46)+E46*(COV46-COV47)+E47*(COV47-COV48)+E48*(COV48-COV49)+E49*COV49

F2 = fixF2
ALAG2 = 9

isolation = 100
infectious = 7

gamma = 1/(infectious)

;Parameter VOC Alpha
StarttagAlpha= 347
kAlpha=0.072
initAlpha=0.002
FractionInfectious_A = 0.08 ;Fraction infectious after vaccination for Alpha+Wildtype
DSMAlpha = TIME-StarttagAlpha
IF(TIME.LT.StarttagAlpha) DSMAlpha = 0

;Parameter VOC Delta
StarttagDelta = 347+154
kDelta= 0.13
initDelta= 0.002
FaktorDelta = 1.97 ;Factor Delta vs WT
FractionInfectious_D = 0.15 ;Fraction infectious after vaccination for Delta
FractionInfectious_BD = 0.15 ;Fraction infectious after booster for Delta
DSMDelta = TIME-StarttagDelta
IF(TIME.LT.StarttagDelta) DSMDelta = 0

;Parameter VOC Omicron
StarttagOmicron = 347+T0
MTIME(10) = StarttagOmicron
kOmicron= 0.21
initOmicron= 0.002
FaktorOmicron = 1.97*1.62 ;Factor Omicron vs WT
FractionInfectious_O = 1 ;Fraction infectious after vaccination for Omicron
FractionInfectious_BO = 0.25 ;Fraction infectious after booster for Omicron
DSMOmicron = (TIME-StarttagOmicron)*MPAST(10)

;Parameter VOC Omicron BA.2
StarttagBA2 = 347+T0
kBA2= 0.0689
initBA2= 0.002
FaktorBA2 = 1.97*1.62*1.5 ;Factor BA2 vs WT
;Parameter VOC Omicron BA.4_5
StarttagBA4_5 = 347+499
kBA4_5= 0.105
initBA4_5= 0.002
FaktorBA4_5 = 1.97*1.62*1.5*1.52   ;Factor BA4_5 vs WT
DSMBA4_5 = (TIME-StarttagBA4_5)*MPAST(37)

AnteilAlpha = 0.92/(1+((1-initAlpha)/initAlpha)*EXP(-kAlpha*DSMAlpha))
AnteilDelta= 0.99/(1+((1-initDelta)/initDelta)*EXP(-kDelta*DSMDelta))
AnteilOmicron= 0.99/(1+((1-initOmicron)/initOmicron)*EXP(-kOmicron*DSMOmicron))*MPAST(10)
AnteilBA4_5 = 0.99/(1+((1-initBA4_5)/initBA4_5)*EXP(-kBA4_5*DSMBA4_5))*MPAST(37)
RR_Inf = FractionInfectious_A*(1-AnteilDelta) + FractionInfectious_D*(AnteilDelta-AnteilOmicron) + FractionInfectious_O*AnteilOmicron
FractionInfectious = FractionInfectious_A*(1-AnteilDelta) + FractionInfectious_D*(AnteilDelta-AnteilOmicron) + FractionInfectious_O*AnteilOmicron
FractionInfectious_Booster = FractionInfectious_A*(1-AnteilDelta) + FractionInfectious_BD*(AnteilDelta-AnteilOmicron) + FractionInfectious_BO*AnteilOmicron

A00toA04_M_COV = A00toA04_M;/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A00toA04_W_COV = A00toA04_W;/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A05toA14_M_COV = A05toA14_M;/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A05toA14_W_COV = A05toA14_W;/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A15toA34_M_COV = A15toA34_M;*(1-Factor15)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A15toA34_W_COV = A15toA34_W;*(1-Factor15)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A35toA59_M_COV = A35toA59_M;*(1-Factor35)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A35toA59_W_COV = A35toA59_W;*(1-Factor35)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A60toA79_M_COV = A60toA79_M;*(1-Factor60)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A60toA79_W_COV = A60toA79_W;*(1-Factor60)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A80_M_COV = A80_M;*(1-Factor80)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A80_W_COV = A80_W;*(1-Factor80)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))

PQ= Positivenquote
IF(Positivenquote.LT.0.1) PQ=0.1

MTIME(6) = THETA(31)
F_INF = THETA(32)
RR_Inf = FractionInfectious_A*(1-AnteilDelta) + FractionInfectious_D*(AnteilDelta-AnteilOmicron)*F_INF**MPAST(6) + FractionInfectious_O*AnteilOmicron
RR_Inf_Booster = FractionInfectious_A*(1-AnteilDelta) + FractionInfectious_BD*(AnteilDelta-AnteilOmicron) + FractionInfectious_BO*AnteilOmicron; 0.08 ;DOI: 10.1056/NEJMoa2101765
F_Hosp_new = THETA(33)
RR_Hosp_WT_60 = 0.12;(1-0.748)
RR_Hosp_WT_00 = 0.08;(1-0.782)
RR_Hosp_D_60 = 0.12 /(1+THETA(28))   ;0.38
RR_Hosp_D_00 = 0.08 /(1+THETA(28))   ;0.38
RR_Hosp_O = THETA(38)
RR_Hosp_60 = (RR_Hosp_WT_60*(1-AnteilDelta) + RR_Hosp_D_60*(AnteilDelta-AnteilOmicron)*F_Hosp_new**MPAST(6) + RR_Hosp_O*AnteilOmicron)
RR_Hosp_00 = (RR_Hosp_WT_00*(1-AnteilDelta) + RR_Hosp_D_00*(AnteilDelta-AnteilOmicron)*F_Hosp_new**MPAST(6) + RR_Hosp_O*AnteilOmicron)
RR_Hosp_60_Booster = (RR_Hosp_WT_60*(1-AnteilDelta) + RR_Hosp_D_60*(AnteilDelta-AnteilOmicron) + RR_Hosp_O*AnteilOmicron)
RR_Hosp_00_Booster = (RR_Hosp_WT_00*(1-AnteilDelta) + RR_Hosp_D_00*(AnteilDelta-AnteilOmicron) + RR_Hosp_O*AnteilOmicron)

FKH = THETA(5)
hill2 = THETA(22)
FKHnTEST = THETA(11)    
hosptrans = THETA(10)

F_Hosp_all = 1
IF(ID.EQ.1) F_Hosp_all = THETA(39)**MPAST(6)
IF(ID.EQ.2) F_Hosp_all = THETA(39)**MPAST(6)
IF(ID.EQ.3) F_Hosp_all = THETA(39)**MPAST(6)
IF(ID.EQ.4) F_Hosp_all = THETA(39)**MPAST(6)
IF(ID.EQ.13) F_Hosp_all = THETA(39)**MPAST(6)
IF(ID.EQ.14) F_Hosp_all = THETA(39)**MPAST(6)
IF(ID.EQ.16) F_Hosp_all = THETA(39)**MPAST(6)
IF(ID.EQ.12) F_Hosp_all = THETA(40)**MPAST(6)
IF(ID.EQ.17) F_Hosp_all = THETA(41)**MPAST(6)

MTIME(46) = 600
hill_2022 = THETA(50)
FHosp_2022 = 1-THETA(49)*(MPAST(46)*(TIME-600))**hill_2022/((MPAST(46)*(TIME-600))**hill_2022+(THETA(46)-600)**hill_2022)
;slider for THETA(51) (see next line): change in hospitalization rate for BA4/5 vs BA2
hospfactor = F_Hosp_all*FHosp_2022*(1-FKH*TIME**hill2/(THETA(4)**hill2+TIME**hill2))*EXP(-FKHnTEST*AnzahlTestungen/1000000)*(1+(AnteilAlpha-AnteilDelta*0.92/0.99)*THETA(20))*(1+AnteilDelta*THETA(28))*(1+AnteilOmicron*THETA(34))*(1+AnteilBA4_5*THETA(51))

;Hosprates decrease with increasing vaccination according to the estimated rate of cases that is vaccinated: 
;RR_Inf*Geimpft80/Agegroup80/(1-Geimpft80/Agegroup80*(1-RR_Inf))
Hosp0 = 0.0735*THETA(1)*hospfactor
Hosp5 = 0.0178*THETA(1)*hospfactor
Hosp15= 0.0305*THETA(1)*hospfactor
Hosp35= 0.148*THETA(1)*hospfactor
Hosp60= 0.578*THETA(1)*hospfactor
IF(Hosp60.GT.1) Hosp60 = 1
Hosp80= THETA(1)*hospfactor
IF(Hosp80.GT.1) Hosp80 = 1

Ungeimpft = 1 - ZweitGesamt - BoosterGesamt
AnteilUngeimpft05 = 1 - Zweit05to14 - Booster05to14
AnteilUngeimpft15 = 1 - Zweit15to59 - Booster15to59
AnteilUngeimpft60 = 1 - Zweit60 - Booster60
All05 = AnteilUngeimpft05+Zweit05to14*RR_Inf+Booster05to14*RR_Inf_Booster
All15 = AnteilUngeimpft15+Zweit15to59*RR_Inf+Booster15to59*RR_Inf_Booster
All60 = AnteilUngeimpft60+Zweit60*RR_Inf+Booster60*RR_Inf_Booster

H0 = Hosp0*(A00toA04_M_COV+A00toA04_W_COV*0.817)
H5 = Hosp5*(A05toA14_M_COV+A05toA14_W_COV)        /All05*(AnteilUngeimpft05+RR_Hosp_00*RR_Inf*Zweit05to14+RR_Hosp_00_Booster*RR_Inf_Booster*Booster05to14)
H15 = Hosp15*(A15toA34_M_COV+A15toA34_W_COV*1.55) /All15*(AnteilUngeimpft15+RR_Hosp_00*RR_Inf*Zweit15to59+RR_Hosp_00_Booster*RR_Inf_Booster*Booster15to59)
H35 = Hosp35*(A35toA59_M_COV+A35toA59_W_COV*0.597)/All15*(AnteilUngeimpft15+RR_Hosp_00*RR_Inf*Zweit15to59+RR_Hosp_00_Booster*RR_Inf_Booster*Booster15to59)
H60 = Hosp60*(A60toA79_M_COV+A60toA79_W_COV*0.708)/All60*(AnteilUngeimpft60+RR_Hosp_60*RR_Inf*Zweit60+RR_Hosp_60_Booster*RR_Inf_Booster*Booster60)
H80 = Hosp80*(A80_M_COV+A80_W_COV*0.588)          /All60*(AnteilUngeimpft60+RR_Hosp_60*RR_Inf*Zweit60+RR_Hosp_60_Booster*RR_Inf_Booster*Booster60)

hosp_TV = H0+H5+H15+H35+H60+H80
hosp = hosp_TV
isolated = (1-hosp)*isolation
hospitalized = hosp*isolation

POP1 = 0
POP2 = 0
IF(ID.EQ.2) POP1 = 1
IF(ID.EQ.3) POP1 = 1
IF(ID.EQ.5) POP1 = 1
IF(ID.EQ.6) POP1 = 1
IF(ID.EQ.7) POP1 = 1
IF(ID.EQ.10) POP1 = 1
IF(ID.EQ.12) POP1 = 1
IF(ID.EQ.17) POP2 = 1

FICU_2022 = 1-THETA(47)*(MPAST(46)*(TIME-600))**hill_2022/((MPAST(46)*(TIME-600))**hill_2022+(THETA(46)-600)**hill_2022)
CPs_ICU = (THETA(6)+(POP1*THETA(13)+(1-POP1)*THETA(24)+POP2*THETA(25))*TIME**hill2/(THETA(4)**hill2+TIME**hill2))*FICU_2022
;slider for THETA(52) (see next line): change in ICU rate for BA4/5 vs BA2
FICU= CPs_ICU*(1+(AnteilAlpha-AnteilDelta*0.92/0.99)*THETA(21))*(1+AnteilDelta*THETA(30))*(1+AnteilOmicron*THETA(42))*(1+AnteilBA4_5*THETA(52))

I0 = 0.439*FICU
IF(I0.GT.1) I0 = 1
I5 = 0.363*FICU
IF(I5.GT.1) I5 = 1
I15 = 0.41*FICU
IF(I15.GT.1) I15 = 1
I35 = 0.66*FICU
IF(I35.GT.1) I35 = 1
I60 = FICU
IF(I60.GT.1) I60 = 1
I80 = 0.65*FICU
IF(I80.GT.1) I80 = 1

RR_ICU_00 = 0.05*THETA(37)
RR_ICU_60 = 0.07*THETA(37)

ICU0 = Hosp0*(A00toA04_M_COV+A00toA04_W_COV*0.817*0.255)*I0/hosp_TV
ICU5 = Hosp5*(A05toA14_M_COV+A05toA14_W_COV*0.299)*I5/hosp_TV/All05*(AnteilUngeimpft05+RR_ICU_00*RR_Hosp_00*RR_Inf*Zweit05to14+RR_ICU_00*RR_Hosp_00_Booster*RR_Inf_Booster*Booster05to14)
ICU15 = Hosp15*(A15toA34_M_COV+A15toA34_W_COV*1.55*0.435)*I15/hosp_TV/All15*(AnteilUngeimpft15+RR_ICU_00*RR_Hosp_00*RR_Inf*Zweit15to59+RR_ICU_00*RR_Hosp_00_Booster*RR_Inf_Booster*Booster15to59)
ICU35 = Hosp35*(A35toA59_M_COV+A35toA59_W_COV*0.597*0.600)*I35/hosp_TV/All15*(AnteilUngeimpft15+RR_ICU_00*RR_Hosp_00*RR_Inf*Zweit15to59+RR_ICU_00*RR_Hosp_00_Booster*RR_Inf_Booster*Booster15to59)
ICU60 = Hosp60*(A60toA79_M_COV+A60toA79_W_COV*0.708*0.652)*I60/hosp_TV/All60*(AnteilUngeimpft60+RR_ICU_60*RR_Hosp_60*RR_Inf*Zweit60+RR_ICU_60*RR_Hosp_60_Booster*RR_Inf_Booster*Booster60)
ICU80 = Hosp80*(A80_M_COV+A80_W_COV*0.588*0.608)*I80/hosp_TV/All60*(AnteilUngeimpft60+RR_ICU_60*RR_Hosp_60*RR_Inf*Zweit60+RR_ICU_60*RR_Hosp_60_Booster*RR_Inf_Booster*Booster60)

revFRACTHOSP =(ICU0+ICU5+ICU15+ICU35+ICU60+ICU80)
FRACTHOSP = 1-revFRACTHOSP
toKH =  hosptrans*FRACTHOSP

; FRACTION ICU (non ventilated)

F_VENT =  1+AnteilOmicron*THETA(44)
B0 = ICU0/revFRACTHOSP*0.308*F_VENT
B5 = ICU5/revFRACTHOSP*0.444*F_VENT
B15 = ICU15/revFRACTHOSP*0.463*F_VENT
B35 = ICU35/revFRACTHOSP*0.600*F_VENT
B60 = ICU60/revFRACTHOSP*0.719*F_VENT
B80 = ICU80/revFRACTHOSP*0.666*F_VENT
FRACTBEAT = (B0+B5+B15+B35+B60+B80)
FRACTICU = 1-FRACTBEAT
toICU = hosptrans*(1-FRACTHOSP)*FRACTICU
toBEAT = hosptrans*(1-FRACTHOSP)*(1-FRACTICU)

KHu60_M = (1-0.66*FICU)*(Hosp35*A35toA59_M_COV/All15*(AnteilUngeimpft15+RR_Hosp_00*RR_Inf*Zweit15to59+RR_Hosp_00_Booster*RR_Inf_Booster*Booster15to59))/hosp_TV/FRACTHOSP
KHue60_M = (1-FICU)*(Hosp60*A60toA79_M_COV/All60*(AnteilUngeimpft60+RR_Hosp_60*RR_Inf*Zweit60+RR_Hosp_60_Booster*RR_Inf_Booster*Booster60))/hosp_TV/FRACTHOSP
KHue80_M = (1-0.65*FICU)*(Hosp80*A80_M_COV/All60*(AnteilUngeimpft60+RR_Hosp_60*RR_Inf*Zweit60+RR_Hosp_60_Booster*RR_Inf_Booster*Booster60))/hosp_TV/FRACTHOSP
KHu60_W = (1-0.600*0.66*FICU)*(Hosp35*A35toA59_W_COV*0.597/All15*(AnteilUngeimpft15+RR_Hosp_00*RR_Inf*Zweit15to59+RR_Hosp_00_Booster*RR_Inf_Booster*Booster15to59))/hosp_TV/FRACTHOSP
KHue60_W = (1-0.652*FICU)*(Hosp60*A60toA79_W_COV*0.708/All60*(AnteilUngeimpft60+RR_Hosp_60*RR_Inf*Zweit60+RR_Hosp_60_Booster*RR_Inf_Booster*Booster60))/hosp_TV/FRACTHOSP
KHue80_W = (1-0.608*0.65*FICU)*(Hosp80*A80_W_COV*0.588/All60*(AnteilUngeimpft60+RR_Hosp_60*RR_Inf*Zweit60+RR_Hosp_60_Booster*RR_Inf_Booster*Booster60))/hosp_TV/FRACTHOSP
ICUue35 = ICU35*(1-0.6*F_VENT)/revFRACTHOSP/FRACTICU
ICUue60 =  ICU60*(1-0.719*F_VENT)/revFRACTHOSP/FRACTICU
ICUue80 =  ICU80*(1-0.666*F_VENT)/revFRACTHOSP/FRACTICU
Bue05 = (B5)/(FRACTBEAT+0.001)
Bue15 = (B15)/(FRACTBEAT+0.001)
Bue35 = (B35)/(FRACTBEAT+0.001)
Bue60 = (B60)/(FRACTBEAT+0.001)
Bue80 = (B80)/(FRACTBEAT+0.001)

PQdeath = (THETA(7)-THETA(7)*THETA(12)*EXP(-PQ*THETA(9)))
PQdeathICU = (THETA(7)-THETA(7)*THETA(12)*EXP(-3.12*THETA(9)))
FDeath_2022 = 1-THETA(48)*(MPAST(46)*(TIME-600))**hill_2022/((MPAST(46)*(TIME-600))**hill_2022+(THETA(46)-600)**hill_2022)
;slider for THETA(53) (see next 4 lines and line 416): change in fatality rate for BA4/5 vs BA2
KHd = (0.0125*KHu60_M+0.0125*KHu60_W+0.147*KHue60_M+0.101*KHue60_W+0.414*KHue80_M+0.334*KHue80_W)*(1+(AnteilAlpha-AnteilDelta*0.92/0.99)*THETA(19))*(1+AnteilDelta*THETA(19))*PQdeath*(1+AnteilOmicron*THETA(35))*(1+AnteilBA4_5*THETA(53))*FDeath_2022
KHd80 = 0.414*PQdeath*(1+(AnteilAlpha-AnteilDelta*0.92/0.99)*THETA(19))*(1+AnteilDelta*THETA(19))*(1+AnteilOmicron*THETA(35))*(1+AnteilBA4_5*THETA(53))*FDeath_2022
ICUd = (0.0453*ICUue35+0.193*ICUue60+0.477*ICUue80) *(1+(AnteilAlpha-AnteilDelta*0.92/0.99)*THETA(19))*(1+AnteilDelta*THETA(19))*PQdeath*(1+AnteilOmicron*THETA(35))*(1+AnteilBA4_5*THETA(53))*FDeath_2022
BEATd = (0.25*Bue05+0.18*Bue15+0.372*Bue35+0.653*Bue60+0.841*Bue80) *(1+(AnteilAlpha-AnteilDelta*0.92/0.99)*THETA(19))*(1+AnteilDelta*THETA(19))*PQdeathICU*(1+AnteilOmicron*THETA(35))*(1+AnteilBA4_5*THETA(53))*FDeath_2022

;MTIME(10) = THETA(18)
hill3 = THETA(23)
FDOUT = THETA(17)*TIME**hill3/(TIME**hill3+THETA(18)**hill3)
hill4 = THETA(29)
FDOUT2 = THETA(26)*THETA(17)*TIME**hill4/(TIME**hill4+THETA(27)**hill4)
RR_Death = 0.16*THETA(36) ;DOI: 10.1056/NEJMoa2101765

out80M = A80_M_COV*((1-Hosp80)*AnteilUngeimpft60+(1-Hosp80*RR_Inf*RR_Hosp_60*RR_Death)*Zweit60+(1-Hosp80*RR_Inf_Booster*RR_Hosp_60_Booster*RR_Death)*Booster60)/All60
out80W = A80_W_COV*((1-Hosp80*0.588)*AnteilUngeimpft60+(1-Hosp80*0.588*RR_Inf*RR_Hosp_60*RR_Death)*Zweit60+(1-Hosp80*0.588*RR_Inf_Booster*RR_Hosp_60_Booster*RR_Death)*Booster60)/All60
out60M = A60toA79_M_COV*((1-Hosp60)*AnteilUngeimpft60+(1-Hosp60*RR_Inf*RR_Hosp_60*RR_Death)*Zweit60+(1-Hosp60*RR_Inf_Booster*RR_Hosp_60_Booster*RR_Death)*Booster60)/All60
out60W = A60toA79_W_COV*((1-Hosp60*0.708)*AnteilUngeimpft60+(1-Hosp60*0.708*RR_Inf*RR_Hosp_60*RR_Death)*Zweit60+(1-Hosp60*0.708*RR_Inf_Booster*RR_Hosp_60_Booster*RR_Death)*Booster60)/All60

OUTd = (THETA(16)+FDOUT-FDOUT2)*PQdeath*(1+(AnteilAlpha-AnteilDelta*0.92/0.99)*THETA(19))*(1+AnteilDelta*THETA(19))*(1+AnteilOmicron*THETA(35))*(1+AnteilBA4_5*THETA(53))*FDeath_2022
death = (out80M+out80W*0.334/0.414+(out60M+out60W*0.101/0.147)*1.22/10.1)*OUTd   ;Verhältnis rate M/W wie Sterberate Normalstation; Verhältnis Ü60/Ü80 = 1.22/10.1 (=Saarland)

;----- Metakis
IIV1 = EXP(ETA(1))
hill1 = THETA(14)
F_LD_BEAT_L = (THETA(2)-(THETA(2)*THETA(3))*TIME**hill1/(TIME**hill1+THETA(15)**hill1))
F_LD_OMICRON = (1+AnteilOmicron*THETA(43))
F_LD_ICU = (1+AnteilOmicron*THETA(45))
;;;;Metakis-Information;;;;
LD_BEAT_D = 15.5  ;Liegedauer Beatmet verstorben
LD_BEAT_L = 28.6*F_LD_BEAT_L ;Liegedauer Beatmet entlassen
LD_ICU_D = 20 ;Liegedauer ICU verstorben
LD_ICU_L = 20.4 ;Liegedauer ICU entlassen
LD_KH_D = 10.6 ;Liegedauer KH verstorben
LD_KH_L = 11.5*F_LD_OMICRON ;Liegedauer KH entlassen
P_BEAT_AUFENTHALT_D = 0.63*F_LD_ICU;Beatmung % Aufenthalt verstorben
P_BEAT_AUFENTHALT_L = 0.28*F_LD_ICU;Beatmung % Aufenthalt entlassen
P_BEAT_on_ICU_AUFENTHALT_D = 0.68*F_LD_ICU;beatmet ICU % Aufenthalt verstorben
P_BEAT_on_ICU_AUFENTHALT_L = 0.43*F_LD_ICU;beatmet ICU % Aufenthalt entlassen
P_ICU_AUFENTHALT_D = 0.44*F_LD_ICU;unbeatmet ICU % Aufenthalt verstorben
P_ICU_AUFENTHALT_L = 0.29*F_LD_ICU;unbeatmet ICU % Aufenthalt entlassen

LDdeathKH = LD_KH_D
LDdeathICU = LD_ICU_D
LDdeathBEAT = LD_BEAT_D
LDrecoveryKH = LD_KH_L
LDrecoveryICU = LD_ICU_L
LDrecoveryBEAT = LD_BEAT_L
deathKH = 2/(LD_KH_D)
deathICU = 2/(LD_ICU_D)
deathBEAT = 2/(LD_BEAT_D)
recovery_rateKH = 2/(LD_KH_L)
recovery_rateICU = 2/(LD_ICU_L)
recovery_rateBEAT = 2/(LD_BEAT_L)

MTIME(5) = 70
REC_TIME = MPAST(5)
recovery = 14
recovery_rate = 2/recovery
recovery_hosp = 14
recovery_ratehosp = 2/recovery_hosp
death_rate = 2/LD_KH_D

KHl = 1-KHd
ICUl = 1-ICUd
BEATl = 1-BEATd

Q_0 = 0
;IF(ID.EQ.2) Q_0=14;ein paar mehr f�r die Bayern
A_0(1) = NUMBER

$DES
DaysSinceAlpha = T-StarttagAlpha
IF(T.LT.StarttagAlpha) DaysSinceAlpha = 0
Alpha = 0.92/(1+((1-initAlpha)/initAlpha)*EXP(-kAlpha*DaysSinceAlpha))
IF(T.LT.StarttagAlpha) Alpha = 0

DaysSinceDelta = T-StarttagDelta
IF(T.LT.StarttagDelta) DaysSinceDelta = 0
Delta = 0.99/(1+((1-initDelta)/initDelta)*EXP(-kDelta*DaysSinceDelta))
IF(T.LT.StarttagDelta) Delta = 0

DaysSinceOmicron = (T-StarttagOmicron)*MPAST(10)
Omicron = 0.99/(1+((1-initOmicron)/initOmicron)*EXP(-kOmicron*DaysSinceOmicron))*MPAST(10)

DaysSinceBA2 = (T-StarttagBA2)*MPAST(10)
BA2 = 0.99/(1+((1-initBA2)/initBA2)*EXP(-kBA2*DaysSinceBA2))*MPAST(10)

DaysSinceBA4_5 = (T-StarttagBA4_5)*MPAST(37)
BA4_5 = 0.99/(1+((1-initBA4_5)/initBA4_5)*EXP(-kBA4_5*DaysSinceBA4_5))*MPAST(37)

R0FREE = (R0FREE_OLD*(1-(Alpha-Delta*0.92/0.99)-Delta)+R0FREE_OLD*1.35*(Alpha-Delta*0.92/0.99)) + R0FREE_OLD*FaktorDelta*(Delta-Omicron) + R0FREE_OLD*FaktorOmicron*Omicron*(1-BA2) + R0FREE_OLD*FaktorBA2*Omicron*(BA2-BA4_5) + R0FREE_OLD*FaktorBA4_5*Omicron*BA4_5
beta = R0FREE*gamma
FractionInfectious_DES = FractionInfectious_A*(1-Delta) + FractionInfectious_D*(Delta-Omicron) + FractionInfectious_O*Omicron
FractionInfectious_Booster_DES = FractionInfectious_A*(1-Delta) + FractionInfectious_BD*(Delta-Omicron) + FractionInfectious_BO*Omicron

;Susceptibles
DADT(1) = - beta*A(1)/NUMBER*A(2)*(Ungeimpft+ZweitGesamt*FractionInfectious_DES+BoosterGesamt*FractionInfectious_Booster_DES)
;Infected
DADT(2) = beta*(A(1)+Omicron*A(30))/NUMBER*A(2)*(Ungeimpft+ZweitGesamt*FractionInfectious_DES+BoosterGesamt*FractionInfectious_Booster_DES) - gamma*A(2) 
;Cumultaive cases
DADT(3) = gamma*A(2) - isolated*A(3) - hospitalized*A(3)
DADT(4) = isolated*A(3)*(1-death) - recovery_rate*A(4) 	            ;R1 
DADT(5) = KHl*toKH*A(26) - recovery_rateKH*A(5) 						;KH L 1
DADT(6) = recovery_rateKH*A(5) - recovery_rateKH*A(6)				;KH L 2
DADT(7) = KHd*toKH*A(26)- deathKH*A(7) 								;KH D 1
DADT(8) = deathKH*A(7) - deathKH*A(8) 								;KH D 1
DADT(9) = ICUl*toICU*A(26) - recovery_rateICU*A(9) 					;ICU L 1
DADT(10) = recovery_rateICU*A(9) - recovery_rateICU*A(10) 			;ICU L 2
DADT(11) = ICUd*toICU*A(26)- deathICU*A(11) 						;ICU D 1
DADT(12) = deathICU*A(11)- deathICU*A(12) 							;ICU D 2
DADT(13) = BEATl*toBEAT*A(26) - recovery_rateBEAT*A(13) 				;BEAT L 1
DADT(14) = recovery_rateBEAT*A(13) - recovery_rateBEAT*A(14) 		;BEAT L 2
DADT(15) = BEATd*toBEAT*A(26) - deathBEAT*A(15) 						;BEAT D 1
DADT(16) = deathBEAT*A(15)- deathBEAT*A(16) 						;BEAT D 2
DADT(17) = recovery_rate*A(20) +  recovery_ratehosp*A(24)           ;RECOVERED
DADT(18) = death_rate*A(29) +deathKH*A(8) + deathICU*A(12) + deathBEAT*A(16);DEAD
DADT(19) = gamma*A(2)                       ;total cases
DADT(20) = recovery_rate*A(4) - recovery_rate*A(20)   ;R2
DADT(21) = (toICU + toBEAT + toKH)*A(26)             ;KH kumulativ
DADT(22) = (toICU + toBEAT)*A(26)                    ;ICU kumulativ
DADT(23) = recovery_rateKH*A(6) + recovery_rateICU*A(10) + recovery_rateBEAT*A(14) - recovery_ratehosp*A(23) ;RKH1
DADT(24) = recovery_ratehosp*A(23) - recovery_ratehosp*A(24) ;RKH2
DADT(25) = hospitalized*A(3) - hosptrans*A(25);DEAD KH
DADT(26) = hosptrans*A(25) - hosptrans*A(26)
DADT(27) = death_rate*A(29) +deathKH*A(8) + deathICU*A(12) + deathBEAT*A(16);DEAD
DADT(28) = isolated*A(3)*death - death_rate*A(28)
DADT(29) = death_rate*A(28) - death_rate*A(29)
;Cases without Omicron --> susceptible to Omicron
DADT(30) = (1-Omicron)*gamma*A(2) - beta*Omicron*A(30)/NUMBER*A(2)*(Ungeimpft+ZweitGesamt*FractionInfectious_DES+BoosterGesamt*FractionInfectious_Booster_DES)
;Cases with Omicron --> immune after Omicron
DADT(31) = Omicron*gamma*A(2)



$ERROR
DaysSinceAlpha_Err = TIME-StarttagAlpha
IF(TIME.LT.StarttagAlpha) DaysSinceAlpha_Err = 0
Alpha_Err = 0.92/(1+((1-initAlpha)/initAlpha)*EXP(-kAlpha*DaysSinceAlpha_Err))
IF(TIME.LT.StarttagAlpha) Alpha_Err = 0
DaysSinceDelta_Err = TIME-StarttagDelta
IF(TIME.LT.StarttagDelta) DaysSinceDelta_Err = 0
Delta_Err = 0.99/(1+((1-initDelta)/initDelta)*EXP(-kDelta*DaysSinceDelta_Err))
IF(TIME.LT.StarttagDelta) Delta_Err = 0
DaysSinceOmicron_Err = (TIME-StarttagOmicron)*MPAST(10)
Omicron_Err = 0.99/(1+((1-initOmicron)/initOmicron)*EXP(-kOmicron*DaysSinceOmicron_Err))*MPAST(10)
DaysSinceBA2_Err = (TIME-StarttagBA2)*MPAST(10)
BA2_Err = 0.99/(1+((1-initBA2)/initBA2)*EXP(-kBA2*DaysSinceBA2_Err))*MPAST(10)
DaysSinceBA4_5_Err = TIME-StarttagBA4_5
BA4_5_Err = 0.99/(1+((1-initBA4_5)/initBA4_5)*EXP(-kBA4_5*DaysSinceBA4_5_Err))
R0 = R0FREE_OLD*(1-(Alpha_Err-Delta_Err*0.92/0.99)-Delta_Err)+R0FREE_OLD*1.35*(Alpha_Err-Delta_Err*0.92/0.99) + R0FREE_OLD*FaktorDelta*(Delta_Err-Omicron_Err) + R0FREE_OLD*FaktorOmicron*Omicron_Err*(1-BA2_Err) + R0FREE_OLD*FaktorBA2*Omicron_Err*(BA2_Err-BA4_5_Err) + R0FREE_OLD*FaktorBA4_5*Omicron_Err*BA4_5_Err
SUSCEPTIBLES = A(1)
FractionInfectious_Err = FractionInfectious_A*(1-Delta_Err) + FractionInfectious_D*(Delta_Err-Omicron_Err) + FractionInfectious_O*Omicron_Err
FractionInfectious_Booster_Err = FractionInfectious_A*(1-Delta_Err) + FractionInfectious_BD*(Delta_Err-Omicron_Err) + FractionInfectious_BO*Omicron_Err
Rt = R0*(SUSCEPTIBLES+A(30)*Omicron_Err)/NUMBER*(Ungeimpft+ZweitGesamt*FractionInfectious_Err+BoosterGesamt*FractionInfectious_Booster_Err)

CASES = A(19)+Q_0
BEATMET = P_BEAT_AUFENTHALT_L*(A(13)+A(14))+P_BEAT_AUFENTHALT_D*(A(15)+A(16))
BEATMET_ABS = BEATMET
BEATMET_ICU = P_BEAT_on_ICU_AUFENTHALT_L*(A(13)+A(14))+P_BEAT_on_ICU_AUFENTHALT_D*(A(15)+A(16))
ICUgesamt = P_ICU_AUFENTHALT_L*(A(9)+A(10))+ P_ICU_AUFENTHALT_D*(A(11)+A(12))+BEATMET_ICU
ICU_ABS = ICUgesamt
KHgesamt = A(5)+A(6)+A(7)+A(8)+A(9)+A(10)+A(11)+A(12)+A(13)+A(14)+A(15)+A(16)
DAILY_DEAD = A(18)
DAILY_KH = A(21)
DEAD = A(27)
ICUkumulativ = A(22)
RECOVERED = CASE0+Q_0+A(17)*REC_TIME
IPRED = CASES
IF(CMT.EQ.4) IPRED = RECOVERED 
IF(CMT.EQ.5) IPRED= ICU_ABS*fBER*fBR*fHH
IF(CMT.EQ.6) IPRED = DEAD
IF(CMT.EQ.7) IPRED = KHgesamt*fBER*fHE*fHH_KH
IF(CMT.EQ.10) IPRED = BEATMET_ABS*fBER*fBR*fHH
IF(CMT.EQ.12) IPRED = ICUkumulativ
IF(CMT.EQ.18) IPRED = DAILY_DEAD
IF(CMT.EQ.21) IPRED = DAILY_KH
W = IPRED
DEL = 0
IF(IPRED.EQ.0) DEL = 0.001
IRES	= DV -  IPRED
IWRES   = IRES/(W+DEL)
Y  = IPRED *EXP(EPS(1)) + EPS(2)
IF(CMT.EQ.4) Y = IPRED + W*EPS(15) + EPS(16)
IF(CMT.EQ.5) Y = IPRED + W*EPS(3) + EPS(4)
IF(CMT.EQ.6) Y = IPRED + W*EPS(5)  + EPS(6)
IF(CMT.EQ.7) Y = IPRED + W*EPS(7) + EPS(8)
IF(CMT.EQ.10) Y = IPRED + W*EPS(9)  + EPS(10)
IF(CMT.EQ.18) Y = IPRED + W*EPS(11) + EPS(12)
IF(CMT.EQ.21) Y = IPRED + W*EPS(13) + EPS(14)

$THETA
(0.8, 0.999,1) FIX ;1 hosp80
(1, 1.69,2) FIX ;2 ftoBEAT1
(0, 0.407,1) FIX ;3 factor ftoBEAT2
(50, 280,400) FIX ;4 MTKH
(0.1, 0.458,1) FIX ;5 fKH
(0, 0.478,1) FIX ;6 toICU 
(0.5, 1.23,1.7) FIX ;7 fdeath
(0.3, 1.52,2) FIX ;8 fBER
(0, 0.0697) FIX ;9 PQ death
(100) FIX ; 10 hosptrans
(0, 0.409) FIX ;11 FKHnTEST
(0, 0.509,1) FIX ;12 fmin PQ Death
(0, 0.171,1) FIX ;13 +toICU POP1
(100) FIX ;14 hill1 vent
(50, 225,350) FIX ;15 MTICU
(0) FIX ;16 dout
(0, 0.27,1) FIX ;17 dout2
(0, 348,450) FIX ;18 MTIME FDOUT
(0) FIX ;19 factor D Alpha
(0, 0.124,1) FIX ;20 factor H Alpha
(0, 0.474,1) FIX ;21 factor ICU Alpha
(100) FIX ;22 hill2 hosp+ICU
(0, 20.3,100) FIX ;23 hill3 dout
(0) FIX ;24 +toICU 1-POP1
(-0.2, 0.0976,1) FIX ;25 +toICU Germany
(0, 0.864,1) FIX ; 26 Anteil -FDOUT
(360, 443,560) FIX ;27 MTIME FDOUT2
(-0.3, 0.284,1) FIX ;28 factor H  Delta
(100) FIX ;29 hill4
(0, 0.546,1) FIX ;30 factor ICU Delta 
(550, 646) FIX ;31 MTIME Inf+Hosp für ungeboostert
(0, 6) FIX ;32 F_INF für ungeboostert
(1, 3.78,6) FIX ;33 F_Hosp für ungeboostert
(-0.9, -0.589,0) FIX ;34 F_hosp_Omicron Anteilige Senkung Hospitalisierungsrate durch Omicron
(-0.9, -0.577,0) FIX ;35 F_death_Omicron Anteilige Senkung Sterberate durch Omicron
(0) FIX ;36 F_RR_dout
(0, 14) FIX ;37 F_RR_ICU
(0, 0.534,1) FIX ;38 F_RR_Hosp_O
(0, 0.888,2) FIX ;39 F_Hosp
(0, 1.08,2) FIX ;40 F_Hosp
(0, 1,2) FIX ;41 F_Hosp
(-0.9, -0.493,0) FIX ;42 F_ICU_Omicron
(0) FIX ;43 F_LD_Omicron
(-0.9, -0.302,0) FIX ;44 F_Vent_Omicron
(0, 0.309,1) FIX ; 45 F_LD_ICU_Om
(600, 792,830) ;46 MT ICU
(0, 0.785,1) ;47 F_ICU2022
(0, 0.539,1) ;48 F_Death2022
(0, 0.308,1) ;49 F_Hosp2022
(0, 8.3,50) ;50 hill_2022
(-1, 0,10) FIX ;51 F_Hosp_BA45 Zugewinn Hospitalisierungsrate durch BA4/5 
(-1, 0,10) FIX ;52 F_ICU_BA45 Zugewinn ICUrate durch BA4/5 
(-1, 0,10) FIX ;53 F_Death_BA45 Zugewinn Sterberate durch BA4/5 

$OMEGA
 0 FIX  ;1 hosptrans

$SIGMA 
 0.000027 ;prop c
 1240 ;add c
 0.0842 ;prop i
 3 FIX  ;add i
 0.00424 ;prop d
 390 ;add d
 0.116 ;prop k
 10 FIX  ;add k
 0.09 ;prop v
 2 FIX  ;add v
 1.03 ;prop dd
 1 FIX  ;add dd
 7.41 ;prop dk
 3.2 ;add dk
 0.1   ;prop r
 10   ;add r

$EST METHOD=1 INTER MAXEVAL=0 NOABORT NSIG=2 SIGL=6 PRINT=1 POSTHOC
;$EST METHOD=SAEM INTERACTION PRINT=100 NBURN=1000 NITER=1000 ISAMPLE=1000
;$EST METHOD=IMP EONLY=1 ISAMPLE=1000 NITER=5 PRINT=1 MAPITER = 0
$COV UNCOND

$TABLE ID TIME DV OBS CMT EVID MDV IPRED Rt R0FREE R0FREE_OLD hosp_TV hospfactor hosp FRACTHOSP FRACTICU KHd ICUd BEATd OUTd KHd80 death ONEHEADER NOPRINT FILE=sdtab2100R
