#define a constellation
#'label',total number, number of planes, sat per plane, inc, alt

# STARLINK Gen1
#csl1a1 = ['SL1 B0',1584,72,22, 53., 550.]
#csl1a3 = ['SL1 B3',375 ,10,38, 81., 560.]
#csl1a4 = ['SL1 B4',450, 12,38, 70., 570.]
#csl1a2 = ['SL1 B2',400 ,16,25, 74., 545.]
#csl1a5 = ['SL1 B1',1600,32,50, 54., 540.]
#csl1b1 = ['SL1 A1',2493,42,60, 42., 336.]
#csl1b2 = ['SL1 A2',2478,42,60, 48., 341.]
#csl1b3 = ['SL1 A3',2547,42,60, 53., 346.]

# updated from SAT-MOD-20200417-00037 - this overwrites prev section (12ksat)
csl1a1 = ['SL1.A1', 1584, 72, 22, 53	, 550.0]
csl1a2 = ['SL1.A2', 1584, 72, 22, 53.2	, 540.0]
csl1a3 = ['SL1.A3', 720 , 36, 20, 70	, 570.0]
csl1a4 = ['SL1.A4', 348 ,  6, 58, 97.6	, 560.0]
csl1a5 = ['SL1.A5', 172 ,  4, 43, 97.6	, 560.0]
csl1b1 = ['SL1.B1', 2493, 42, 60, 42	, 335.9]
csl1b2 = ['SL1.B2', 2478, 42, 60, 48	, 340.8]
csl1b3 = ['SL1.B3', 2547, 42, 60, 53	, 345.6]

# STARLINK Gen2
# random orbits grouped in planes
csl2a1 = ['SL2.A1',7178,84,86, 30.,  328.]
csl2a2 = ['SL2.A2',7178,84,86, 40.,  334.]
csl2a3 = ['SL2.A3',7178,84,86, 53.,  345.]
csl2a4 = ['SL2.A4',2000,40,50, 96.9, 360.]
csl2a5 = ['SL2.A5',1998,20,100,75., 373.]
csl2b0 = ['SL2.B' ,4000,40,100,53., 499.]
csl2c1 = ['SL2.C1',144 ,12, 12,148., 604.]
csl2c2 = ['SL2.C2',324 ,18, 18,116., 614.]


#TEST CONSTELLATION FRANZ
#satname,0,planes,sats,inc,alt
fsl1 = ['FSL1',7140,84,85,30.,328.]
fsl2 = ['FSL2',7140,84,85,40.,334.]
fsl3 = ['FSL3',7140,84,85,53.,345.]
fsl4 = ['FSL4',2000,20,100,75.,373.]
fsl5 = ['FSL5',4000,40,100,53.,499.]
fsl6 = ['FSL6',4000,40,100,53.,499.]
fsl7 = ['FSL7',144,12,12,148.,604.]
fsl8 = ['FSL8',324,18,18,116.,614.]
fsl9 = ['FSL9',2000,40,50,97.,360.]
franz = [ fsl1, fsl2, fsl3, fsl4, fsl5, fsl6, fsl7, fsl8, fsl9]


# OneWeb
# OW2 original
cowa  = ['OW2 A',  1764, 36, 49,87.9,1200.]
cowb  = ['OW2 B', 23040,128,180,40.,1200.]
cowc  = ['OW2 C', 23040,128,180,55.,1200.]
# OW2-reduced 2021
cow2a  = ['OW2 A', 1764,36,49,87.9,1200.]
cow2b  = ['OW2 B', 2304,32,72,40.,1200.]
cow2c  = ['OW2 C', 2304,32,72,55.,1200.]



# Guo Wang (China) GW-A59 -2	    12,992 				
cgw11 = ['GW-A59.1',  480,  8, 60, 85.,  590.] 
cgw12 = ['GW-A59.2', 2000, 40, 50, 50.,  600.] 
cgw13 = ['GW-A59.3', 3600, 60, 60, 55.,  508.] 
cgw21 = ['GW.2.1',   1728, 27, 64, 30., 1145.] 
cgw22 = ['GW.2.2',   1728, 27, 64, 40., 1145.] 
cgw23 = ['GW.2.3',   1728, 27, 64, 50., 1145.] 
cgw24 = ['GW.2.4',   1728, 27, 64, 60., 1145.]

# Amazon Kuiper	 Total: 	3236
#'label',total number, number of planes, sat per plane, inc, alt

cak1 = ['AK.',  784, 28, 28, 33., 590.] 
cak2 = ['AK.', 1296, 36, 36, 42., 610.] 
cak3 = ['AK.', 1156, 34, 34, 52., 630.] 

# SPECIAL CONSTELLATIONS:
todayHigh1  =  ['U-LEO1', 1110 ,60,20, 70.,1055.]
todayHigh2  =  ['U-LEO2', 1110 ,60,20, 60.,1055.]
todayLow    =  ['L-LEO',   505, 50,10, 60., 505.]
todaySL     =  ['SL'  ,   1635, 72,22, 53., 550.]
todayOW     =  ['OW2 A',   128, 36, 49,87.9,1200.]


SL1  = [ csl1a1, csl1a2, csl1a3, csl1a4, csl1a5, csl1b1, csl1b2, csl1b3 ]
SL2  = [ csl2a1, csl2a2, csl2a3, csl2a4, csl2a5, csl2b0, csl2c1, csl2c2 ]
OW2  = [ cowa, cowb, cowc]
OW2r = [ cow2a, cow2b, cow2c]
GW   = [ cgw11, cgw12, cgw13, cgw21, cgw22, cgw23, cgw24]
AK   = [cak1, cak2, cak3]


constellations = [ csl2a1, csl2a2, csl2a3, csl2a4, csl2a5, csl2b0, csl2c1, csl2c2,
                   cowa, cowb, cowc] # for testing legacy SLOW2 vs SatCon report and such

yesturday = [todayHigh1, todayHigh2, todayLow]
constoday = [todaySL, todayOW]
SLtoday = [todaySL]
