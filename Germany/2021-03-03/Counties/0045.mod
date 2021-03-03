;; 1. Based on: 0044
;; 2. Description: FIX RKI Modell (Ges.Dtl), open IIVs, 24 Effects RKI (fix 19), Mutante ab 347
;; x1. Author: dings

;; 3. Label:

$SIZES NO=10000 PD=-100 LVR=50

$PROBLEM COVID-19 cases

$INPUT ID TIME DROP Landkreis=DROP IdLandkreis Bundesland=DROP ID_BL DV AMT MDV CMT Week OBS EWZ DV_CASE0 SCHOOLS STORECLOSURE CURFEW RESTRAININGORDER CARNIVAL day=DROP MDVWKND=DROP fixR0_0 fixR0_1 fixR0_2 fixR0_3 fixR0_4 fixR0_5 fixR0_6 fixR0_7 fixR0_8 fixR0_9 fixR0_10 fixR0_11 fixR0_12 fixR0_13 fixR0_14 fixR0_15 fixR0_16 fixR0_17 fixR0_18 fixR0_19 fixR0_20=DROP fixR0_21=DROP fixCP_4 fixCP_5 fixCP_6 fixCP_7 fixCP_8 fixCP_9 fixCP_10 fixCP_11 fixCP_12 fixCP_13 fixCP_14 fixCP_15 fixCP_16 fixCP_17 fixCP_18 fixCP_19 fixCP_20=DROP fixCP_21=DROP fixF2 PREDIC
$DATA ..\..\DATASET\Landkreise_V4_2021_03_03_R0_0044.csv IGNORE=@ IGNORE(PREDIC.EQ.1) IGNORE(CMT.GT.3)

$SUBROUTINES ADVAN13 TOL=6
$MODEL NCOMPS=3

$PK
NUMBER = EWZ
Q_0 = 0
;;;;Metakis-Information;;;;;
P_KH = 0.801 ;Prozent Normalstation
P_ICU = 0.252 ;Prozent ICU vs BEAT
P_BEAT_D = 0.514;Sterberate Beatmet
P_ICU_D = 0.204 ;Sterberate ICU
P_KH_D = 0.155 ;Sterberate KH
LD_BEAT_D = 13 ;Liegedauer Beatmet verstorben
LD_BEAT_L = 35;Liegedauer Beatmet entlassen
LD_ICU_D = 13;Liegedauer ICU verstorben
LD_ICU_L = 19;Liegedauer ICU entlassen
LD_KH_D = 8.8;Liegedauer KH verstorben
LD_KH_L = 12;Liegedauer KH entlassen
P_BEAT_AUFENTHALT_D = 0.76;Beatmung % Aufenthalt verstorben
P_BEAT_AUFENTHALT_L = 0.49;Beatmung % Aufenthalt entlassen
P_BEAT_on_ICU_AUFENTHALT_D = 0.88;beatmet ICU % Aufenthalt verstorben
P_BEAT_on_ICU_AUFENTHALT_L =  0.75;beatmet ICU % Aufenthalt entlassen
P_ICU_AUFENTHALT_D = 0.57;unbeatmet ICU % Aufenthalt verstorben
P_ICU_AUFENTHALT_L = 0.34;unbeatmet ICU % Aufenthalt entlassen

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
MTIME(5) = fixCP_8
COV8 = MPAST(5)
E8 = fixR0_8
MTIME(7) = fixCP_9
COV9 = MPAST(7)
E9 = fixR0_9
MTIME(8) = fixCP_10
COV10 = MPAST(8)
E10 = fixR0_10
MTIME(9) = fixCP_11
COV11 = MPAST(9)
E11 =fixR0_11
MTIME(10) = fixCP_12
COV12 = MPAST(10)
E12 = fixR0_12
MTIME(11) = fixCP_13
COV13 = MPAST(11)
E13 = fixR0_13
MTIME(12) = fixCP_14
COV14 = MPAST(12)
E14 = fixR0_14
MTIME(13) = fixCP_15
COV15 = MPAST(13)
E15 =  fixR0_15
MTIME(14) = fixCP_16
COV16 = MPAST(14)
E16 = fixR0_16

MTIME(15) = fixCP_17
COV17 = MPAST(15)
E17 = fixR0_17

MTIME(16) = fixCP_18
COV18 = MPAST(16)
E18 = fixR0_18

MTIME(17) = fixCP_19
COV19 = MPAST(17)
E19 = fixR0_19

MTIME(18) = MTIME(17)+THETA(7)
COV20 = MPAST(18)
E20 = THETA(8)*EXP(ETA(4))

MTIME(19) = MTIME(18)+THETA(9)
COV21 = MPAST(19)
E21 = THETA(10)*EXP(ETA(5))
MTIME(20) = MTIME(19)+THETA(11)
COV22 = MPAST(20)
E22 = THETA(12)*EXP(ETA(6))

MTIME(21) = MTIME(20)+THETA(13)
COV23 = MPAST(21)
E23 = THETA(14)*EXP(ETA(7))
MTIME(22) = MTIME(21)+THETA(15)
COV24 = MPAST(22)
E24 = E23+ THETA(16)*EXP(ETA(8))

R0FREE_OLD = R0START*(1-COV1)+E2*COV2+E1*(COV1-COV3)+E3*(COV3-COV4)+E4*(COV4-COV5)+E5*(COV5-COV6)+E6*(COV6-COV7)+E7*(COV7-COV8)+E8*(COV8-COV9)+E9*(COV9-COV10)+E10*(COV10-COV11)+E11*(COV11-COV12)+E12*(COV12-COV13)+E13*(COV13-COV14)+E14*(COV14-COV15)+E15*(COV15-COV16)+E16*(COV16-COV17)+E17*(COV17-COV18)+E18*(COV18-COV19)+E19*(COV19-COV20)+E20*(COV20-COV21)+E21*(COV21-COV22)+E22*(COV22-COV23)+E23*(COV23-COV24)+E24*COV24
R0BASE = R0START
R0END = E23
R0_0 = R0START ;Startwert
R0_1 = E1 ;Startwert
R0_2 = E2 ;nach Schulschließung
R0_3 = E3 ;nach Kontaktverbot
R0_4 = E4 ;
CP_4 = MTIME(1)
R0_5 = E5 ;
CP_5 = MTIME(2)
R0_6 = E6 ;
CP_6 = MTIME(3)
R0_7 = E7 ;
CP_7 = MTIME(4)
R0_8 = E8 ;
CP_8 = MTIME(5)
R0_9 = E9 ;
CP_9 = MTIME(7)
R0_10 = E10 ;
CP_10 = MTIME(8)
R0_11 = E11 ;
CP_11 = MTIME(9)
R0_12 = E12 ;
CP_12 = MTIME(10)
R0_13 = E13 ;
CP_13 = MTIME(11)
R0_14 = E14 ;
CP_14 = MTIME(12)
R0_15 = E15 ;
CP_15 = MTIME(13)
R0_16 = E16 ;
CP_16 = MTIME(14)
R0_17 = E17 ;
CP_17 = MTIME(15)
R0_18 = E18 ;
CP_18 = MTIME(16)
R0_19 = E19 ;
CP_19 = MTIME(17)
R0_20 = E20 ;
CP_20 = MTIME(18)
R0_21 = E21 ;
CP_21 = MTIME(19)
R0_22 = E22 ;
CP_22 = MTIME(20)
R0_23 = E23 ;
CP_23 = MTIME(21)
R0_24 = E24 ;
CP_24 = MTIME(22)

F2 = fixF2
ALAG2 = 9

isolation = 100
infectious = 7

gamma = 1/(infectious)

MTIME(6) = 1
MTDIFF = 1
REC_TIME = MPAST(6)

A_0(1) = NUMBER

k=0.072
initMutante=0.002
Starttag = 347

$DES
DaysSinceMutante = T-Starttag
IF(T.LT.Starttag) DaysSinceMutante = 0
Mutante = 1/(1+((1-initMutante)/initMutante)*EXP(-k*DaysSinceMutante))
IF(T.LT.Starttag) Mutante = 0
R0FREE = R0FREE_OLD*(1-Mutante)+R0FREE_OLD*1.35*Mutante
beta = R0FREE*gamma

DADT(1) = -beta*A(1)/NUMBER*A(2) 
DADT(2) = beta*A(1)/NUMBER*A(2) -  gamma*A(2) 
DADT(3) = gamma*A(2)

$ERROR
CASES = A(3)+Q_0
IPRED = CASES
W = IPRED
DEL = 0
IF(IPRED.EQ.0) DEL = 0.001
IRES	= DV -  IPRED
IWRES   = IRES/(W+DEL)
Y  = IPRED *EXP(EPS(1)) + EPS(2)

$THETA
(0, 9.32) FIX;1 +MTIME17
(0.5, 0.933) FIX;2 EFFECT17
(0, 8.21) FIX;3 +MTIME18
(0.5, 1.23) FIX;4 EFFECT18
(0, 19) FIX;5 +MTIME19
(0.5, 0.673) FIX;6 EFFECT19
(0, 10.3) FIX;7 +MTIME20
(0.5, 1.23) FIX;8 EFFECT20
(0, 10.3) FIX;9 +MTIME21
(0.5, 0.781) FIX;10 EFFECT21
(0, 10.2) FIX;11 +MTIME22
(0.5, 0.84) FIX;12 EFFECT22
(0, 11.1) FIX;11 +MTIME23
(0.5, 0.728) FIX;12 EFFECT23
(0, 10) FIX;11 +MTIME24
(0, 0.187) FIX;12 EFFECT24

$OMEGA 
 0 FIX;1 EFFECT17
 0 FIX;2 EFFECT18
 0 FIX;3 EFFECT19
 0.0059 FIX;4 EFFECT20
 0.0127 FIX;5 EFFECT21
 0.0307 FIX;6 EFFECT22
 0.0038 FIX;7 EFFECT23
 0 FIX;8 EFFECT24

$SIGMA 
 0.1 ;prop c
10 ;add c

$EST METHOD=1 INTER MAXEVAL=9999 NOABORT NSIG=3 SIGL=6 PRINT=1 POSTHOC

$COV UNCONDITIONAL

$TABLE ID TIME DV OBS CMT MDV IPRED IdLandkreis R0_0 R0_1 R0_2 R0_3 R0_4 R0_5 R0_6 R0_7 R0_8 R0_9 R0_10 R0_11 R0_12 R0_13 R0_14 R0_15 R0_16 R0_17 R0_18 R0_19 R0_20 R0_21 R0_22 R0_23 R0_24 CP_4 CP_5 CP_6 CP_7 CP_8 CP_9 CP_10 CP_11 CP_12 CP_13 CP_14 CP_15 CP_16 CP_17 CP_18 CP_19 CP_20 CP_21 CP_22 CP_23 CP_24 F2 R0FREE R0END ONEHEADER NOPRINT FILE=sdtab0045
