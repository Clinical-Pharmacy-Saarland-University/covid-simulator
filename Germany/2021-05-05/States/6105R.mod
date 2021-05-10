;; 1. Based on: 6105G
;; 2. Description: clean, Impfungen(-87%hosp,Impfwillige,Korrektur Erstimpfungen), MT LD BEAT L, hosp MT+n(Tests), MT PQ dout+dKH max-min, hosp+ICU(Mutante), to ICU groups each(+Ba), each hill, real lowering hosprate
;; x1. Author: dings

;; 3. Label:
;; Kommentar: Für die Simulationen zukünftiger Impfungen muss die Anzahl an Impfungen auf die Erstimpfung und mit 6 Wochen Zeitverzug auf die Zweitimpfungen addiert werden

$SIZES NO=-10000 PD=-5000 LVR=50 PC=31

$PROBLEM COVID-19 cases

$INPUT ID DROP STATE=DROP Week Year TIME TIME2 DV LNDV LNDVCR OBS CMT AMT EVID MDV NUMBER CASE0 ICUDATASET DIVIFLAG SCHOOLS STORECLOSURE CURFEW RESTRAININGORDER CARNIVAL FLAGZEROOBS FLAGRKI FLAGMP DOSING2 MAXIMUMICU FLAG20200714 MDV20200714 FLAG20200804 RKIDEATHS MPCASES MPCASESRECOVRKIDEATH Zweitimpfungen secondDosesCumulative A00toA04_M A00toA04_W A05toA14_M A05toA14_W A15toA34_M A15toA34_W A35toA59_M A35toA59_W A60toA79_M A60toA79_W A80_M A80_W Positivenquote AnzahlTestungen PREDIC SIMzweit SIMzweitCUM initialDosesCumulative fixR0_0 fixR0_1 fixR0_2 fixR0_3 fixR0_4 fixR0_5 fixR0_6 fixR0_7 fixR0_8 fixR0_9 fixR0_10 fixR0_11 fixR0_12 fixR0_13 fixR0_14 fixR0_15 fixR0_16 fixR0_17 fixR0_18 fixR0_19 fixR0_20 fixR0_21 fixR0_22 fixR0_23 fixR0_24 fixR0_25 fixR0_26 fixR0_27 fixR0_28 fixCP_4 fixCP_5 fixCP_6 fixCP_7 fixCP_8 fixCP_9 fixCP_10 fixCP_11 fixCP_12 fixCP_13 fixCP_14 fixCP_15 fixCP_16 fixCP_17 fixCP_18 fixCP_19 fixCP_20 fixCP_21 fixCP_22 fixCP_23 fixCP_24 fixCP_25 fixCP_26 fixCP_27 fixCP_28 fixF2 MDVLAST
$DATA ..\DATASET\COVID_ds_166_Impf_AGE_SEX_TESTS_BL_DAY_D_KHNEW_PRED_IMPF_R0_6000_GERMANY33.csv IGNORE=@ IGNORE(DOSING2.EQ.1) IGNORE(DIVIFLAG.EQ.2) IGNORE(CMT.EQ.8) IGNORE(CMT.EQ.11) IGNORE(CMT.EQ.12) IGNORE(CMT.EQ.9) IGNORE(FLAGZEROOBS.EQ.1) IGNORE(ID.EQ.17) IGNORE(FLAGRKI.EQ.1) IGNORE(PREDIC.EQ.1)

$SUBROUTINES ADVAN6 TOL=6
$MODEL NCOMPS=30

$PK
IIV1 = EXP(ETA(1))
hill1 = THETA(14)
F_LD_BEAT_L = THETA(2)-(THETA(2)*THETA(3))*TIME**hill1/(TIME**hill1+THETA(15)**hill1)
;;;;Metakis-Information;;;;;
LD_BEAT_D = 15.5  ;Liegedauer Beatmet verstorben
LD_BEAT_L = 28.6*F_LD_BEAT_L ;Liegedauer Beatmet entlassen
LD_ICU_D = 20 ;Liegedauer ICU verstorben
LD_ICU_L = 20.4 ;Liegedauer ICU entlassen
LD_KH_D = 10.6 ;Liegedauer KH verstorben
LD_KH_L = 11.5 ;Liegedauer KH entlassen
P_BEAT_AUFENTHALT_D = 0.63;Beatmung % Aufenthalt verstorben
P_BEAT_AUFENTHALT_L = 0.28;Beatmung % Aufenthalt entlassen
P_BEAT_on_ICU_AUFENTHALT_D = 0.68;beatmet ICU % Aufenthalt verstorben
P_BEAT_on_ICU_AUFENTHALT_L = 0.43;beatmet ICU % Aufenthalt entlassen
P_ICU_AUFENTHALT_D = 0.44;unbeatmet ICU % Aufenthalt verstorben
P_ICU_AUFENTHALT_L = 0.29;unbeatmet ICU % Aufenthalt entlassen

fBER = 1
IF(ID.EQ.3) fBER = THETA(8)
fHE = 1
IF(ID.EQ.7) fHE = 1.15
fBR = 1
IF(ID.EQ.5) fBR = 1.3
fHH = 1
IF(ID.EQ.6) fHH = 1.2

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

R0FREE_OLD = R0START*(1-COV1)+E2*COV2+E1*(COV1-COV3)+E3*(COV3-COV4)+E4*(COV4-COV5)+E5*(COV5-COV6)+E6*(COV6-COV7)+E7*(COV7-COV8)+E8*(COV8-COV9)+E9*(COV9-COV10)+E10*(COV10-COV11)+E11*(COV11-COV12)+E12*(COV12-COV13)+E13*(COV13-COV14)+E14*(COV14-COV15)+E15*(COV15-COV16)+E16*(COV16-COV17)+E17*(COV17-COV18)+E18*(COV18-COV19)+E19*(COV19-COV20)+E20*(COV20-COV21)+E21*(COV21-COV22)+E22*(COV22-COV23)+E23*(COV23-COV24)+E24*(COV24-COV25)+E25*(COV25-COV26)+E26*(COV26-COV27)+E27*(COV27-COV28)+E28*COV28

F2 = fixF2
ALAG2 = 9

MTIME(6) = 1

isolation = 100
infectious = 7

gamma = 1/(infectious)
FractionInfectious = 0.08
FractionUninfetious = 1-FractionInfectious

TglImfpungabsim = SIMzweit
SimulatedImpfungen = SIMzweitCUM
SIMtglImpfungen = TglImfpungabsim
FractionAge = 0.547 ;=1-Anteil Berufe
FractionImpfwillig = 0.956
Fract15 = 0.230
Fract35 = 0.348
Fract60 = 0.21
Fract80 = 0.068

Agegroup80 = NUMBER*Fract80
Agegroup60 = NUMBER*Fract60
Agegroup35 = NUMBER*Fract35
Agegroup15 = NUMBER*Fract15

RelativeGeimpft80 = (secondDosesCumulative  + SimulatedImpfungen)*FractionAge/Agegroup80
IF(RelativeGeimpft80.GT.FractionImpfwillig) RelativeGeimpft80 = FractionImpfwillig
RelativeGeimpft60 = ((secondDosesCumulative  + SimulatedImpfungen)*FractionAge - Agegroup80*FractionImpfwillig)/Agegroup60
IF(RelativeGeimpft60.LT.0) RelativeGeimpft60 = 0
IF(RelativeGeimpft60.GT.FractionImpfwillig) RelativeGeimpft60 = FractionImpfwillig
RelativeGeimpft35 = ((secondDosesCumulative + SimulatedImpfungen)*(1-FractionAge))/Agegroup35
IF(RelativeGeimpft60.GE.FractionImpfwillig) RelativeGeimpft35 =  (secondDosesCumulative + SimulatedImpfungen - (Agegroup80 + Agegroup60)*FractionImpfwillig)/Agegroup35
IF(RelativeGeimpft35.GT.FractionImpfwillig) RelativeGeimpft35 = FractionImpfwillig
RelativeGeimpft15 = 0
IF(RelativeGeimpft35.GE.FractionImpfwillig) RelativeGeimpft15 = (secondDosesCumulative + SimulatedImpfungen - (Agegroup80 + Agegroup60 + Agegroup35)*FractionImpfwillig)/Agegroup15
IF(RelativeGeimpft15.GT.FractionImpfwillig) RelativeGeimpft15 = FractionImpfwillig
RelativeGeimpft = (secondDosesCumulative + SimulatedImpfungen)/NUMBER
IF(RelativeGeimpft.GT.FractionImpfwillig) RelativeGeimpft = FractionImpfwillig

RelativeErstgeimpft80 = initialDosesCumulative*FractionAge/Agegroup80
IF(RelativeErstgeimpft80.GT.FractionImpfwillig) RelativeErstgeimpft80 = FractionImpfwillig
RelativeErstgeimpft60 = (initialDosesCumulative*FractionAge - Agegroup80*FractionImpfwillig)/Agegroup60
IF(RelativeErstgeimpft60.LT.0) RelativeErstgeimpft60 = 0
IF(RelativeErstgeimpft60.GT.FractionImpfwillig) RelativeErstgeimpft60 = FractionImpfwillig
RelativeErstgeimpft35 = (initialDosesCumulative*(FractionImpfwillig-FractionAge))/Agegroup35
IF(RelativeErstgeimpft60.GE.FractionImpfwillig) RelativeErstgeimpft35 =  (initialDosesCumulative - (Agegroup80 + Agegroup60)*FractionImpfwillig)/Agegroup35
IF(RelativeErstgeimpft35.GT.FractionImpfwillig) RelativeErstgeimpft35 = FractionImpfwillig
RelativeErstgeimpft15 = 0
IF(RelativeErstgeimpft35.GE.FractionImpfwillig) RelativeErstgeimpft15 = (initialDosesCumulative - (Agegroup80 + Agegroup60 + Agegroup35)*FractionImpfwillig)/Agegroup15
IF(RelativeErstgeimpft15.GT.FractionImpfwillig) RelativeErstgeimpft15 = FractionImpfwillig
RelativeErstgeimpft = initialDosesCumulative/NUMBER
IF(RelativeErstgeimpft.GT.FractionImpfwillig) RelativeErstgeimpft = FractionImpfwillig

Geimpft80 = (secondDosesCumulative + SimulatedImpfungen)*FractionAge
Neugeimpft80 = SimulatedImpfungen*FractionAge
IF(RelativeGeimpft80.GE.FractionImpfwillig) Geimpft80 = Agegroup80*FractionImpfwillig
Geimpft35 = (secondDosesCumulative + SimulatedImpfungen)*(1-FractionAge)
Neugeimpft35 =  SimulatedImpfungen*(1-FractionAge)
IF(RelativeGeimpft60.GE.FractionImpfwillig) Geimpft35 = ((secondDosesCumulative + SimulatedImpfungen)*(1-FractionAge) + (secondDosesCumulative + SimulatedImpfungen)*FractionAge - Agegroup80*FractionImpfwillig - Agegroup60*FractionImpfwillig)
IF(RelativeGeimpft35.GT.FractionImpfwillig) Geimpft35 = Agegroup35*FractionImpfwillig
Geimpft60 = 0
Neugeimpft60 = 0
IF(RelativeGeimpft80.GE.FractionImpfwillig) Geimpft60 = ((secondDosesCumulative + SimulatedImpfungen)*FractionAge - Agegroup80*FractionImpfwillig)
IF(RelativeGeimpft80.GE.FractionImpfwillig) Neugeimpft60 = ((secondDosesCumulative + SimulatedImpfungen)*FractionAge - Agegroup80*FractionImpfwillig)
IF(RelativeGeimpft80.GE.FractionImpfwillig.AND.RelativeGeimpft60.GE.FractionImpfwillig) Geimpft60 = Agegroup60*FractionImpfwillig
;IF(RelativeGeimpft80.GE.FractionImpfwillig.AND.RelativeGeimpft60.GE.FractionImpfwillig) Neugeimpft60 = Agegroup60*FractionImpfwillig
Geimpft15 = 0
Neugeimpft15 = 0
IF(RelativeGeimpft35.GE.FractionImpfwillig) Geimpft15 = (secondDosesCumulative+SimulatedImpfungen) - Agegroup80*FractionImpfwillig - Agegroup60*FractionImpfwillig - Agegroup35*FractionImpfwillig
IF(RelativeGeimpft35.GE.FractionImpfwillig) Neugeimpft15 = (secondDosesCumulative+SimulatedImpfungen) - Agegroup80*FractionImpfwillig - Agegroup60*FractionImpfwillig - Agegroup35*FractionImpfwillig
IF(RelativeGeimpft15.GE.FractionImpfwillig) Geimpft15 = Agegroup15*FractionImpfwillig
IF(RelativeGeimpft15.GE.FractionImpfwillig) Neugeimpft15 = Agegroup15*FractionImpfwillig

Vacc15 = 0
IF(RelativeGeimpft80.GE.FractionImpfwillig.AND.RelativeGeimpft60.GE.FractionImpfwillig.AND.RelativeGeimpft35.GE.FractionImpfwillig) Vacc15 = 1
Vacc35 = (1-FractionAge)
IF(RelativeGeimpft80.GE.FractionImpfwillig.AND.RelativeGeimpft60.GE.FractionImpfwillig) Vacc35 = 1
Vacc60 = 0
IF(RelativeGeimpft80.GE.FractionImpfwillig) Vacc60 = FractionAge
Vacc80 = FractionAge

RR_Inf = 0.08 

;Factor dafür wie viele Zweitimpfungen in der Zukunft dazukommen
Factor80= (1-RR_Inf)*Neugeimpft80/(Agegroup80-secondDosesCumulative*FractionAge)
IF(RelativeGeimpft80.GE.FractionImpfwillig) Factor80 = (1-RR_Inf)*(FractionImpfwillig*Agegroup80-secondDosesCumulative*FractionAge)/(Agegroup80-secondDosesCumulative*FractionAge)
;Wenn mit SCHON GESCHEHENEN IMPFUNGEN die Altersgruppe durchgeimpft ist muss der Faktor auf 0 begrenzt werden.
IF(RelativeGeimpft80.GE.FractionImpfwillig.AND.FractionImpfwillig*Agegroup80.LE.secondDosesCumulative*FractionAge) Factor80 = 0
Factor60= (1-RR_Inf)*Neugeimpft60/Agegroup60
;Wenn Ü80 schon durchgeimpft sind müssen die Impfungen in die Ü60 
IF(FractionImpfwillig*Agegroup80.LE.secondDosesCumulative*FractionAge) Factor60 = (1-RR_Inf)*Neugeimpft60/(Agegroup60-(secondDosesCumulative*FractionAge-FractionImpfwillig*Agegroup80))
IF(RelativeGeimpft60.GE.FractionImpfwillig) Factor60 = (1-RR_Inf)*FractionImpfwillig*Neugeimpft60/Agegroup60
IF(RelativeGeimpft60.GE.FractionImpfwillig.AND.FractionImpfwillig*Agegroup80.LE.secondDosesCumulative*FractionAge) Factor60 = (1-RR_Inf)*(Agegroup60*FractionImpfwillig-(secondDosesCumulative*FractionAge-FractionImpfwillig*Agegroup80))/(Agegroup60-(secondDosesCumulative*FractionAge-FractionImpfwillig*Agegroup80))
;IF(Factor60.LE.0) Factor60 = 0
Factor35= (1-RR_Inf)*Neugeimpft35/(Agegroup35-secondDosesCumulative*(1-FractionAge))
IF(RelativeGeimpft35.GE.FractionImpfwillig) Factor35 = (1-RR_Inf)*FractionImpfwillig*Agegroup35/(Agegroup35-secondDosesCumulative*(1-FractionAge))
;IF(Factor60.LE.0) Factor35 = (1-RR_Inf)*FractionImpfwillig*Agegroup35/(Agegroup35-(secondDosesCumulative-Agegroup80-Agegroup60))
Factor15= (1-RR_Inf)*Neugeimpft15/Agegroup15
IF(RelativeGeimpft15.GE.FractionImpfwillig) Factor15= (1-RR_Inf)*FractionImpfwillig*Agegroup15/Agegroup15

A00toA04_M_COV = A00toA04_M/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A00toA04_W_COV = A00toA04_W/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A05toA14_M_COV = A05toA14_M/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A05toA14_W_COV = A05toA14_W/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A15toA34_M_COV = A15toA34_M*(1-Factor15)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A15toA34_W_COV = A15toA34_W*(1-Factor15)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A35toA59_M_COV = A35toA59_M*(1-Factor35)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A35toA59_W_COV = A35toA59_W*(1-Factor35)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A60toA79_M_COV = A60toA79_M*(1-Factor60)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A60toA79_W_COV = A60toA79_W*(1-Factor60)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A80_M_COV = A80_M*(1-Factor80)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))
A80_W_COV = A80_W*(1-Factor80)/(1-Factor80*(A80_M+A80_W)-Factor60*(A60toA79_M+A60toA79_W)-Factor35*(A35toA59_M+A35toA59_W)-Factor15*(A15toA34_M+A15toA34_W))

PQ= Positivenquote
IF(Positivenquote.LT.0.1) PQ=0.1

RR_Inf = 0.08 ;DOI: 10.1056/NEJMoa2101765
RR_Hosp = 0.13 ;DOI: 10.1056/NEJMoa2101765

;Hosprates decrease with increasing vaccination according to the estimated rate of cases that is vaccinated: 
;RR_Inf*Geimpft80/Agegroup80/(1-Geimpft80/Agegroup80*(1-RR_Inf))
Hosp80= THETA(1)*(1-(1-RR_Hosp)*RR_Inf*Geimpft80/Agegroup80/(1-Geimpft80/Agegroup80*(1-RR_Inf)))
Hosp60= 0.578*THETA(1)*(1-(1-RR_Hosp)*RR_Inf*Geimpft60/Agegroup60/(1-Geimpft60/Agegroup60*(1-RR_Inf)))
Hosp0 = 0.0735*THETA(1)
Hosp5 = 0.0178*THETA(1)
Hosp15= 0.0305*THETA(1)*(1-(1-RR_Hosp)*RR_Inf*Geimpft15/Agegroup15/(1-Geimpft15/Agegroup15*(1-RR_Inf)))
Hosp35= 0.148*THETA(1)*(1-(1-RR_Hosp)*RR_Inf*Geimpft35/Agegroup35/(1-Geimpft35/Agegroup35*(1-RR_Inf)))

H0 = Hosp0*(A00toA04_M_COV+A00toA04_W_COV*0.817)
H5 = Hosp5*(A05toA14_M_COV+A05toA14_W_COV)
H15 = Hosp15*(A15toA34_M_COV+A15toA34_W_COV*1.55)
H35 = Hosp35*(A35toA59_M_COV+A35toA59_W_COV*0.597)
H60 = Hosp60*(A60toA79_M_COV+A60toA79_W_COV*0.708)
H80 = Hosp80*(A80_M_COV+A80_W_COV*0.588)

Starttag = 347
k=0.072
initMutante=0.002
DSM = TIME-Starttag
IF(TIME.LT.Starttag) DSM = 0
AnteilMutante = 1/(1+((1-initMutante)/initMutante)*EXP(-k*DSM))

MTIME(15) = THETA(10)
FKH = THETA(5)
hill2 = THETA(22)
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
FICU= (THETA(6)+(POP1*THETA(13)+(1-POP1)*THETA(24)+POP2*THETA(25))*TIME**hill2/(THETA(4)**hill2+TIME**hill2))*(1+AnteilMutante*THETA(21))
FKHnTEST = THETA(11)    

hosptrans = THETA(10)
hosp_TV = H0+H5+H15+H35+H60+H80
hospfactor = (1-FKH*TIME**hill2/(THETA(4)**hill2+TIME**hill2))*EXP(-FKHnTEST*AnzahlTestungen/1000000)*(1+AnteilMutante*THETA(20))
hosp = hosp_TV*hospfactor
isolated = (1-hosp)*isolation
hospitalized = hosp*isolation

ICU0 = Hosp0*(A00toA04_M_COV+A00toA04_W_COV*0.817*0.255)*0.439*FICU/hosp_TV
ICU5 = Hosp5*(A05toA14_M_COV+A05toA14_W_COV*0.299)*0.363*FICU/hosp_TV
ICU15 = Hosp15*(A15toA34_M_COV+A15toA34_W_COV*1.55*0.435)*0.41*FICU/hosp_TV
ICU35 = Hosp35*(A35toA59_M_COV+A35toA59_W_COV*0.597*0.600)*0.66*FICU/hosp_TV
ICU60 = Hosp60*(A60toA79_M_COV+A60toA79_W_COV*0.708*0.652)*FICU/hosp_TV
ICU80 = Hosp80*(A80_M_COV+A80_W_COV*0.588*0.608)*0.65*FICU/hosp_TV

revFRACTHOSP =(ICU0+ICU5+ICU15+ICU35+ICU60+ICU80)
FRACTHOSP = 1-revFRACTHOSP
toKH =  hosptrans*FRACTHOSP

; FRACTION ICU (non ventilated)
B0 = ICU0/revFRACTHOSP*0.308
B5 = ICU5/revFRACTHOSP*0.444
B15 = ICU15/revFRACTHOSP*0.463
B35 = ICU35/revFRACTHOSP*0.600
B60 = ICU60/revFRACTHOSP*0.719
B80 = ICU80/revFRACTHOSP*0.666
FRACTBEAT = (B0+B5+B15+B35+B60+B80)
FRACTICU = 1-FRACTBEAT
toICU = hosptrans*(1-FRACTHOSP)*FRACTICU
toBEAT = hosptrans*(1-FRACTHOSP)*(1-FRACTICU)

KHu60_M = (1-0.66*FICU)*(Hosp35*A35toA59_M_COV)/hosp_TV/FRACTHOSP
KHue60_M = (1-FICU)*(Hosp60*A60toA79_M_COV)/hosp_TV/FRACTHOSP
KHue80_M = (1-0.65*FICU)*(Hosp80*A80_M_COV)/hosp_TV/FRACTHOSP
KHu60_W = (1-0.600*0.66*FICU)*(Hosp35*A35toA59_W_COV*0.597)/hosp_TV/FRACTHOSP
KHue60_W = (1-0.652*FICU)*(Hosp60*A60toA79_W_COV*0.708)/hosp_TV/FRACTHOSP
KHue80_W = (1-0.608*0.65*FICU)*(Hosp80*A80_W_COV*0.588)/hosp_TV/FRACTHOSP
ICUue35 = ICU35*(1-0.6)/revFRACTHOSP/FRACTICU
ICUue60 =  ICU60*(1-0.719)/revFRACTHOSP/FRACTICU
ICUue80 =  ICU80*(1-0.666)/revFRACTHOSP/FRACTICU
Bue05 = (B5)/(FRACTBEAT+0.001)
Bue15 = (B15)/(FRACTBEAT+0.001)
Bue35 = (B35)/(FRACTBEAT+0.001)
Bue60 = (B60)/(FRACTBEAT+0.001)
Bue80 = (B80)/(FRACTBEAT+0.001)

PQdeath = (THETA(7)-THETA(7)*THETA(12)*EXP(-PQ*THETA(9)))
PQdeathICU = (THETA(7)-THETA(7)*THETA(12)*EXP(-3.12*THETA(9)))
KHd = (0.0125*KHu60_M+0.0125*KHu60_W+0.147*KHue60_M+0.101*KHue60_W+0.414*KHue80_M+0.334*KHue80_W)*(1+AnteilMutante*THETA(19))*PQdeath
ICUd = (0.0453*ICUue35+0.193*ICUue60+0.477*ICUue80) *(1+AnteilMutante*THETA(19))*PQdeathICU
BEATd = (0.25*Bue05+0.18*Bue15+0.372*Bue35+0.653*Bue60+0.841*Bue80) *(1+AnteilMutante*THETA(19))*PQdeathICU

MTIME(10) = THETA(18)
hill3 = THETA(23)
FDOUT = THETA(17)*TIME**hill3/(TIME**hill3+THETA(18)**hill3)
FDOUT2 = THETA(26)*THETA(17)*TIME**hill3/(TIME**hill3+THETA(27)**hill3)
RR_Death = 0.16 ;DOI: 10.1056/NEJMoa2101765
out80M = A80_M_COV*(1-Hosp80*hospfactor)*(1-(1-RR_Death)*RR_Inf*Geimpft80/Agegroup80)/(1-Geimpft80/Agegroup80*(1-RR_Inf))
out80W = A80_W_COV*(1-Hosp80*0.588*hospfactor)*(1-(1-RR_Death)*RR_Inf*Geimpft80/Agegroup80)/(1-Geimpft80/Agegroup80*(1-RR_Inf))
out60M = A60toA79_M_COV*(1-Hosp60*hospfactor)*(1-(1-RR_Death)*RR_Inf*Geimpft60/Agegroup60)/(1-Geimpft60/Agegroup60*(1-RR_Inf))
out60W = A60toA79_W_COV*(1-Hosp60*0.708*hospfactor)*(1-(1-RR_Death)*RR_Inf*Geimpft60/Agegroup60)/(1-Geimpft60/Agegroup60*(1-RR_Inf))
OUTd = (THETA(16)+FDOUT-FDOUT2)*PQdeath*(1+AnteilMutante*THETA(19))
death = (out80M+out80W*0.334/0.414+(out60M+out60W*0.101/0.147)*1.22/10.1)/(1-hosp)*OUTd   ;Verhältnis rate M/W wie Sterberate Normalstation; Verhältnis Ü60/Ü80 = 1.22/10.1 (=Saarland)
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
DaysSinceMutante = T-Starttag
IF(T.LT.Starttag) DaysSinceMutante = 0
Mutante = 1/(1+((1-initMutante)/initMutante)*EXP(-k*DaysSinceMutante))
IF(T.LT.Starttag) Mutante = 0
R0FREE = R0FREE_OLD*(1-Mutante)+R0FREE_OLD*1.35*Mutante
beta = R0FREE*gamma
DADT(1) = -beta*A(1)/NUMBER*A(2) - (Zweitimpfungen+SIMtglImpfungen)*FractionUninfetious
DADT(2) = beta*A(1)/NUMBER*A(2) -  gamma*A(2) 
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


$ERROR
DaysSinceMutante_Err = TIME-Starttag
IF(TIME.LT.Starttag) DaysSinceMutante_Err = 0
Mutante_Err = 1/(1+((1-initMutante)/initMutante)*EXP(-k*DaysSinceMutante_Err))
IF(TIME.LT.Starttag) Mutante_Err = 0
R0FREE_Err = R0FREE_OLD*(1-Mutante_Err)+R0FREE_OLD*1.35*Mutante_Err
Rt = R0FREE_Err*A(1)/NUMBER
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
IF(CMT.EQ.7) IPRED = KHgesamt*fBER*fHE
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
(0, 0.616,1) FIX ;3 factor ftoBEAT2
(50, 280,400) FIX ;4 MTKH
(0.1, 0.478,1) FIX ;5 fKH
(0, 0.478,1) FIX ;6 toICU 
(0.5, 1.05,1.7) FIX ;7 fdeath
(0.3, 1.52,2) FIX ;8 fBER
(0, 0.132) FIX ;9 PQ death
(100) FIX ; 10 hosptrans
(0, 0.334) FIX ;11 FKHnTEST
(0, 0.484,1) FIX ;12 fmin PQ Death
(0, 0.255,1) FIX ;13 +toICU POP1
(100) FIX ;14 hill1 vent
(50, 225,350) FIX ;15 MTICU
(0) FIX ;16 dout
(0, 0.222,1) FIX ;17 dout2
(0, 348,450) FIX ;18 MTIME FDOUT
(0) FIX ;19 factor D Mutante
(0, 0.0836,1) FIX ;20 factor H Mutante
(0, 0.43,1) FIX ;21 factor ICU Mutante
(100) FIX ;22 hill2 hosp+ICU
(0, 29.6,100) FIX ;23 hill3 dout
(-0.2, 0.0588,1) FIX ;24 +toICU 1-POP1
(-0.2, 0.137,1) FIX;25 +toICU Germany
(1) FIX ; 26 Anteil -FDOUT
(360, 446,500) FIX ;27 MTIME FDOUT2

$OMEGA 
 0 FIX  ;1 hosptrans

$SIGMA 
 0.000011 ;prop c
 216000 ;add c
 0.0178 ;prop i
 3 FIX  ;add i
 0.000166 ;prop d
 38200 ;add d
 0.157 ;prop k
 10 FIX  ;add k
 0.0389 ;prop v
 2 FIX  ;add v
 0.234 ;prop dd
 1 FIX  ;add dd
 6.62 ;prop dk
 2.12 ;add dk
 0.1   ;prop r
 10   ;add r

$EST METHOD=1 INTER MAXEVAL=100 NOABORT NSIG=3 SIGL=6 PRINT=1 POSTHOC
;$EST METHOD=SAEM INTERACTION PRINT=100 NBURN=1000 NITER=1000 ISAMPLE=1000
;$EST METHOD=IMP EONLY=1 ISAMPLE=1000 NITER=5 PRINT=1 MAPITER = 0
;$COV UNCOND

$TABLE ID TIME DV OBS CMT EVID MDV IPRED DIVIFLAG Rt Hosp80 A80_W_COV A80_M_COV RelativeGeimpft80 RelativeGeimpft60 RelativeGeimpft Factor80 Factor60 Factor35 Factor15 secondDosesCumulative Agegroup80 R0FREE R0FREE_OLD hosp FRACTHOSP FRACTICU KHd ICUd BEATd OUTd death ONEHEADER NOPRINT FILE=sdtab6105R
