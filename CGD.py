########################################################################################### 
# Copyright (c),A Vipin Menon and BIG LAB in Hanyang University(HYU) 			  #
# Author A Vipin Menon 									  #
# Date 18th September, 2019								  #	
# Email a.vipin.menon@gmail.com				 				  #
###########################################################################################

import argparse
import sys,os 
import math,RNA		
import  numpy as np		 
from itertools import chain
from collections import OrderedDict
from get_sequence import reverseComp
def CINDELcalculatescore(seq):
	params = [('A',13,0.359277769),('A',28,0.128732671),('A',2,0.089207956),('A',31,0.050508798),('A',7,0.063364247),('A',8,0.193517666),('AA',0,0.074052742),('AA',9,-0.218926976),('AA',10,-0.667432796),('AA',11,-0.521040934),('AA',12,-0.264806378),('AA',19,-0.295058409),('AA',21,-0.17718561),('AA',25,0.186318341),('AA',28,0.176168855),('AA',2,0.11203447),('AC',0,0.119809836),('AC',9,-0.427869295),('AC',13,0.176192535),('AC',18,0.101264344),('AC',1,-0.057025748),('AC',19,0.123838693),('AC',20,0.12864349),('AC',27,-0.115956893),('AC',30,-0.211710227),('AC',31,0.055096611),('AC',32,-0.066466648),('AC',7,0.284818922),('AC',8,0.168980629),('AG',12,0.107615207),('AG',13,0.116406611),('AG',14,0.085343537),('AG',23,-0.052145867),('AG',25,-0.147681142),('AG',32,-0.251946146),
('AG',7,0.269547232),('AT',12,-0.055738967),('AT',14,-0.280193723),('AT',16,-0.072299285),('AT',17,-0.188148058),('AT',21,0.13577806),('AT',24,0.071551024),('AT',28,0.121715007),('AT',30,-0.128900537),('C',15,0.187281105),('C',16,0.050394743),('C',17,0.046877049),
('C',28,-0.046360015),('C',2,-0.139751449),('C',29,-0.108997003),('C',3,-0.58046435),('C',7,0.282845387),('CA',0,-0.036713083),('CA',9,-0.129619149),('CA',10,-0.228614104),('CA',13,-0.095178226),('CA',1,0.163224633),('CA',19,0.056594182),('CA',20,0.197202638),
('CA',21,0.027273837),('CA',22,-0.139136158),('CA',7,0.154463045),('CA',8,-0.112498711),('CC',13,0.129634035),('CC',14,0.053211288),('CC',15,0.056682693),('CC',16,0.099331396),('CC',23,-0.099912427),('CC',28,-0.109688314),('CC',2,-0.207256935),('CC',29,-0.047640371),
('CC',32,-0.05591437),('CC',8,0.65410187),('CG',10,-0.193574201),('CG',12,-0.085785037),('CG',13,-0.231165614),('CG',14,-0.11446628),('CG',15,-0.119233092),('CG',17,-0.188450362),('CG',18,-0.090577903),('CG',1,0.755899442),('CG',25,-0.103394296),('CG',26,-0.235079641),
('CG',27,-0.196914463),('CG',28,-0.102209147),('CG',2,0.885302791),('CG',32,0.072689139),('CG',8,-0.108629339),('CT',9,0.348496508),('CT',11,-0.224637913),('CT',17,0.054202198),('CT',18,0.081685704),('CT',26,0.033291728),('CT',27,0.119643634),('CT',3,-0.183633657),
('G',12,0.027338979),('G',16,-0.053448765),('G',18,-0.073831779),('G',24,-0.137566753),('G',27,-0.350861759),('G',2,0.036045258),('G',3,0.096654517),('G',8,0.826090966),('GA',9,0.359935685),('GA',10,0.078569053),('GA',12,0.169220542),('GA',13,0.058267212),('GA',16,-0.096677117),('GA',25,0.183733491),('GA',30,0.049311266),('GA',31,-0.066188348),('GA',32,0.07418699),('GC',10,-0.372731969),('GC',11,-0.072920116),('GC',12,-0.18347655),('GC',16,0.07684379),('GC',27,-0.207964823),('GC',28,-0.068003699),('GC',30,-0.169119607),
('GC',31,-0.084231075),('GG',10,0.061275537),('GG',11,0.363683574),('GG',12,0.107234473),('GG',13,0.2762344),('GG',15,-0.130283364),('GG',16,-0.210606843),('GG',22,-0.278365734),('GG',23,-0.392061429),('GG',24,-0.375657208),('GG',25,-0.276708759),('GG',26,-0.053345569),('GG',32,-0.039499132),('GG',7,0.056485458),('GT',9,0.191596031),('GT',11,0.166350263),('GT',12,0.087029602),('GT',20,0.086971691),('GT',22,0.055831476),('GT',23,0.04992352),('GT',26,-0.036547343),('GT',27,-0.227481683),('GT',28,-0.086145062),('GT',29,-0.316010288),
('GT',7,-0.162305425),('GT',8,0.191117867),('T',14,-0.174791539),('T',15,-0.064452789),('T',17,-0.124207175),('T',20,-0.035009888),('T',3,0.233384052),('T',7,-3.464862095),('T',8,-0.184797232),('TA',10,0.311652992),('TA',14,-0.085557183),('TA',16,0.097333539),
('TA',18,-0.206798454),('TA',24,0.235579018),('TA',28,0.366518978),('TA',2,-0.457700388),('TA',6,0.064310276),('TA',7,0.460015141),('TC',12,-0.058018063),('TC',16,-0.058701337),('TC',18,0.212081096),('TC',1,-0.108103079),('TC',19,0.134732673),('TC',23,0.045835459),
('TC',26,0.245536487),('TC',2,-0.0720811),('TC',29,0.229540677),('TC',30,0.035413033),('TC',6,0.091569704),('TC',7,-0.651599392),('TG',9,0.312514593),('TG',10,0.321561498),('TG',1,0.312174546),('TG',21,0.041875516),('TG',22,0.088872094),('TG',24,0.039480307),('TG',31,-0.099295609),('TT',0,-0.089958105),('TT',9,-0.746490329),('TT',10,-0.500057876),('TT',11,-0.577720541),('TT',12,-0.806473095),('TT',13,-0.752573367),('TT',14,-0.717943612),('TT',15,-0.648672115),('TT',16,-0.353673401),('TT',17,-0.497720724),
('TT',18,-0.474783261),('TT',1,-0.089333716),('TT',19,-0.557002415),('TT',20,-0.553738685),('TT',21,-0.581588808),('TT',22,-0.442802126),('TT',24,-0.347367156),('TT',25,-0.320309377),('TT',28,0.135942601),('TT',3,0.027108425),('TT',6,-0.304449945),('TT',7,-0.443442269),('TT',8,-1.188734393)]
	intercept = -1.139975209
        Free_Energy = 0.187763689
	Entropy =  1.125236597
	AC = 0.023788317
	AG = 0.013631088
	CC = -0.162320093
	CG = 0.082469431
	CT = 0.039421446
	GC = -0.022381161
	GT = 0.017215267
	TA = 0.084467023
        score  = intercept
	gcLow = -0.00000000000499
	gcHigh = 0.666115596
        #### feature extraction 
	##Free_energy
        Energy = seq[8:31]
        Energycal = RNA.fold(Energy)[-1]
        Energycal = round(Energycal,0)
        score +=(Energycal*Free_Energy)
	entropy = dict()
        Entropycal = []
        entropy_seq = seq[4:27]
        lentseq= len(entropy_seq)
	ntscount = {'A':0, 'G':0, 'T':0, 'C':0}
	nuc = ['A','T','G','C']
        for ant in nuc:
                ntscount[ant]=(entropy_seq.count(ant))/float((lentseq))
        for ant in nuc:
                if ntscount[ant]!=0:
                       entropy[ant] = -(ntscount[ant]*np.log2(ntscount[ant]))
                else:
                        entropy[ant] = 0
	entropySum = sum(entropy.values())
	entropySumR = round(entropySum,1)
	score += (entropySumR*Entropy)
	# independent nucleotide composition
	score  = score + AG*(seq.count('AG'))
	score = score + AC*(seq.count('AC'))
	score = score + CG*(seq.count('CG'))
	score = score + CC*(seq.count('CC'))
	score = score + CT*(seq.count('CT'))
	score = score + GC*(seq.count('GC'))
	score = score + GT*(seq.count('GT'))
	score = score + TA*(seq.count('TA'))

	guideSeq = seq[8:31]
    	gcCount = guideSeq.count("G") + guideSeq.count("C")
    	if gcCount <= 9:
        	gcWeight = gcLow
    	if gcCount > 9:
        	gcWeight = gcHigh
    	score += abs(9-gcCount)*gcWeight
        for bp, position, wt in params:
                subSeq = seq[position:position+len(bp)]
                if subSeq==bp:
                        score += wt
        return (1.0/(1.0 + math.exp(-score)))

def LINDELcalculatescore(seq):
	parameters = [('A',13,0.1352634445),('A',18,0.2266601183),('A',27,-0.102297985),('T',24,-0.1930543379),('G',4,0.407190242),('G',7,0.2783643584),('G',10,0.0697306019),('G',17,-0.3592921504),('G',18,-0.1287897125),('G',20,0.3691969814),('G',22,0.4827249725),('G',23,1.3789539029),('G',24,0.1259575325),('G',27,-0.2625186279),('G',6,-0.1587906028),('C',21,0.4349944451),('C',27,0.38874478),('AA',18,0.2991883708),('AA',21,-0.4959624082),('TA',17,0.4111542835),('TA',19,-0.2819993431),('TA',20,-1.2192705551),('TA',22,-0.3685329773),('GA',18,0.2622360001),('CA',21,0.3077715591),('CA',22,0.4931654331),('CA',27,0.3573515163),('AT',21,-0.3234310359),('TT',2,0.3494023736),('TT',4,-0.3881133926),('TT',8,-0.2643292347),('TT',13,-0.7733015532),('TT',14,-0.5396509519),('TT',15,-0.2632898469),('TT',19,-0.48701914),('TT',20,-1.6781779772),('TT',21,-1.113814519),('TT',22,-0.9512292771),('TT',23,-0.2276175337),('GT',6,0.3553403317),('GT',17,-0.2247142002),('GT',18,-0.2918636039),('GT',22,1.0425529535),('CT',16,0.2859762363),('CT',17,0.1824673344),('CT',22,-0.7220324868),('CT',23,-0.1886224787),('AG',19,0.301297765),('AG',20,-0.5703135186),('TG',0,-0.1225163563),('TG',9,0.2277042429),('TG',12,-0.2042563995),('TG',16,-0.1867730992),('TG',21,0.3007303605),('GG',12,0.2508979853),('GG',19,-0.3387841934),('CG',16,-0.2579818758),('CG',17,-0.1991746456),('CG',20,0.2107592806),('CG',27,0.2430978627),('AC',14,0.1378218027),('AC',20,0.9705751273),('AC',21,0.7595488933),('TC',5,-0.1842673099),('TC',11,-0.1419670918),('TC',22,-0.6059995886),('GC',9,0.2384841836),('GC',10,0.318694388),('GC',17,-0.2323498858),('GC',18,-0.158530008),('GC',19,0.2060297787),('CC',4,-0.2574206023),('CC',9,-0.1012745036),('CC',23,0.296799014),('CC',28,-0.3137204638)]
	intercept = 0.675532755
        Free_Energy = 0.168073832
	Entropy =  0.4037850136
	gchigh = 0.0346278599
	gclow = -2.13E-006
	T = -0.0478914162
	AC = 0.0706769127
	CG = -0.0786235386
	GT = 0.0751977394
	TT = -0.2414317576
	score  = intercept
        #### feature extraction 
	##Free_energy
        Energy = seq[4:24]
        Energycal = RNA.fold(Energy)[-1]
        Energycal = round(Energycal,0)
        score +=(Energycal*Free_Energy)
	entropy = dict()
        Entropycal = []
        entropy_seq = seq[4:24]
        lentseq= len(entropy_seq)
	ntscount = {'A':0, 'G':0, 'T':0, 'C':0}
	nuc = ['A','T','G','C']
        for ant in nuc:
                ntscount[ant]=(entropy_seq.count(ant))/float((lentseq))
        for ant in nuc:
                if ntscount[ant]!=0:
                        entropy[ant] = -(ntscount[ant]*np.log2(ntscount[ant]))
                else:
                        entropy[ant] = 0
	entropySum = sum(entropy.values())
	entropySumR = round(entropySum,1)
	score += (entropySumR*Entropy)
	# independent nucleotide composition
	score  = score + AC*(seq.count('AC'))
	score = score + T*(seq.count('T'))
	score = score + CG*(seq.count('CG'))
	score = score + GT*(seq.count('GT'))
	score = score + TT*(seq.count('TT'))

	guideSeq = seq[4:24]
        gcCount = guideSeq.count("G") + guideSeq.count("C")
        if gcCount <= 10:
                gcWeight = gclow
        if gcCount > 10:
                gcWeight = gchigh
        score += abs(10-gcCount)*gcWeight

        for bp, position, wt in parameters:
                subSeq = seq[position:position+len(bp)]
                if subSeq==bp:
                        score += wt
        return (1.0/(1.0 + math.exp(-score)))

		 
def LINDELacalculatescore(seq):
	parameters = [("A",20,0.137170389),
("T",1,-0.074585823),
("T",7,0.111059881),
("T",29,0.098637193),
("C",7,-0.049469777),
("C",9,0.120651133),
("C",13,-0.076869685),
("C",15,0.095490951),
("AA",13,0.441176916),
("AA",23,0.151540679),
("TA",1,-0.087314896),
("TA",4,0.417718773),
("TA",7,0.105648267),
("TA",8,0.570685229),
("TA",10,-0.665609037),
("TA",14,-0.185567376),
("TA",22,0.869906711),
("TA",27,0.108646417),
("GA",4,-0.130676946),
("GA",8,-0.471092337),
("GA",10,-0.190239216),
("GA",14,0.158674623),
("GA",19,0.251114014),
("GA",21,0.286665059),
("GA",25,0.414530408),
("CA",18,-0.422451619),
("CA",28,0.476939986),
("AT",8,-0.10205527),
("AT",11,-0.135439919),
("AT",12,-0.340822076),
("AT",15,-0.220870944),
("AT",21,0.76460505),
("AT",23,-1.125123608),
("TT",5,-0.258567901),
("TT",8,0.684879543),
("TT",11,0.263153909),
("TT",13,-0.226042759),
("GT",0,-0.145425996),
("GT",2,-0.056128866),
("GT",6,0.293513389),
("GT",7,0.116560365),
("GT",17,-0.557781356),
("CT",4,-0.143478539),
("CT",21,-0.127121675),
("CT",22,0.146064932),
("CT",23,-0.054563533),
("CT",25,0.253908979),
("CT",27,-0.160033658),
("AG",1,0.181299386),
("AG",3,0.088189755),
("AG",16,-0.169852859),
("TG",8,-0.629010616),
("TG",9,-0.532089435),
("TG",13,0.116569465),
("TG",21,-0.088926423),
("GG",0,0.097573308),
("GG",4,0.089778821),
("GG",13,0.043998425),
("GG",17,0.051130828),
("GG",20,0.183442263),
("CG",12,0.584212426),
("AC",5,-0.487374162),
("AC",6,-0.191273931),
("AC",7,0.849177822),
("AC",12,0.277918839),
("AC",19,-0.091508321),
("TC",4,-0.802078354),
("TC",28,-0.144629743),
("GC",7,-0.052976793),
("GC",15,0.360078914),
("GC",22,0.079236972),
("CC",0,0.075396858),
("CC",5,-0.050765507),
("CC",19,0.115084718),
("CC",22,-0.089198303)]
	intercept = -1.671946104
	Entropy =  0.4037850136
	gchigh = 0.240138978
	gclow = -0.016037794
	Free_Energy = -0.000342883
	GG = 0.090650765
	AT = -0.040546895
	MT = 0.005148235
	score  = intercept
        #### feature extraction 
	##Free_energy
	##Free_energy
        Energy = seq[4:24]
        Energycal = RNA.fold(Energy)[-1]
        Energycal = round(Energycal,0)
        score +=(Energycal*Free_Energy)
	# independent nucleotide composition
	score  = score + AT*(seq.count('AT'))
	score = score + GG*(seq.count('GG'))
	guideSeq = seq[4:24]
        gcCount = guideSeq.count("G") + guideSeq.count("C")
        if gcCount <= 10:
                gcWeight = gclow
        if gcCount > 10:
                gcWeight = gchigh
        score += abs(10-gcCount)*gcWeight
	Melt_temp = (float(64.9) + float(41.0) * ((float(guideSeq.count('G')) + float(guideSeq.count('C')) - float(16.4))/(20.0)))
	score+= Melt_temp*MT
        for bp, position, wt in parameters:
                subSeq = seq[position:position+len(bp)]
                if subSeq==bp:
                        score += wt
        return (1.0/(1.0 + math.exp(-score)))

def LINDELicalculatescore(seq):

        parameters = [('A',18,0.0815450954),('T',3,0.1251411781),('T',14,-0.05601527),('T',16,-0.1374303579),('T',26,0.6129942086),('G',3,-0.1584387612),('C',22,-0.0687106001),('C',23,-0.1450156038),('C',27,0.1564440206),('C',28,-0.1007689224),('AA',9,-0.1547718709),('AA',11,-0.3012529992),('AA',22,0.1655011092),('TA',8,-0.1481922885),('GA',13,0.1244079638),('GA',18,0.1880705815),('CA',20,-0.3439729695),('CA',21,0.1416105185),('CA',22,0.1522284385),('CA',27,0.1847169166),('CA',28,-0.1975980794),('AT',2,0.3944690605),('AT',8,-0.1972874534),('AT',11,-0.2336100532),('AT',22,-0.2545299014),('AT',27,-0.3396280329),('TT',5,-0.3974796405),('TT',6,-0.0900360596),('TT',7,-0.2787390725),('TT',9,-0.1649851365),('TT',10,-0.2572279591),('TT',11,-0.302645449),('TT',13,-0.3616339187),('TT',14,-0.265153483),('TT',20,-0.293255811),('TT',21,-0.4450813624),('TT',23,-0.2013282932),('GT',0,0.1335707585),('GT',1,-0.1045879319),('GT',7,0.1046889318),('CT',15,-0.1706747529),('CT',23,-0.2440952346),('AG',13,0.0674383094),('AG',21,-0.0548264226),('TG',16,-0.0621481838),('GG',18,-0.3465682018),('GG',19,-0.1770629408),('GG',26,-0.207419321),('GG',27,-0.2534565649),('CG',1,0.1439811675),('CG',7,0.1884952631),('CG',10,0.1179399121),('CG',16,-0.2754010247),('CG',17,-0.1893117011),('AC',15,0.2027064141),('AC',18,0.3311347746),('AC',20,0.2303312739),('AC',28,0.1772425831),('TC',9,-0.1755042319),('GC',1,0.1250831227),('GC',5,0.1045485972),('GC',11,0.0851007713),('GC',19,-0.1945267787),('GC',20,-0.1458484707),('GC',21,-0.5363020499),('GC',22,-0.5256951798),('GC',23,0.1888301689),('GC',25,-0.7042919152),('CC',18,0.1522234528),('CC',22,-0.1439032381),('CC',28,-0.1504334951)]

        intercept = -1.3484915738
        Free_energy = 0.0584654915
	Entropy = 0.4056274813
        GChigh = 0.7542669585
        GClow = -0.0065689225
        TT = -0.1045974512
        AT = -0.0957803804
        AG = 0.1051405001
        GG = 0.0459548209
	GT = 0.0463509282
	AA = -0.0437729377
        TA = 0.1324070584
        score = intercept 


        #### feature extraction 
        ##Free_energy
        Energy = seq[4:24]
        Energycal = RNA.fold(Energy)[-1]
        Energycal = round(Energycal,0)
        score +=(Energycal*Free_energy)
        # independent nucleotide composition
        score = score + AG*(seq.count('AG'))
        score = score + AT*(seq.count('AT'))
        score = score + GG*(seq.count('GG'))
        score = score + TT*(seq.count('TT'))
        score = score + TA*(seq.count('TA'))
	score = score + AA*(seq.count('AA'))
	score = score + GT*(seq.count('GT'))
	#entropy cal
	entropy = dict()
        Entropycal = []
        entropy_seq = seq[4:24]
        lentseq= len(entropy_seq)
	ntscount = {'A':0, 'G':0, 'T':0, 'C':0}
	nuc = ['A','T','G','C']
        for ant in nuc:
                ntscount[ant]=(entropy_seq.count(ant))/float((lentseq))
        for ant in nuc:
                if ntscount[ant]!=0:
                        entropy[ant] = -(ntscount[ant]*np.log2(ntscount[ant]))
                else:
                        entropy[ant] = 0
	entropySum = sum(entropy.values())
	entropySumR = round(entropySum,1)
	score += (entropySumR*Entropy)
	
	guideSeq = seq[4:24]
        gcCount = guideSeq.count("G") + guideSeq.count("C")
        if gcCount <= 10:
                gcWeight = GClow
        if gcCount > 10:
                gcWeight = GChigh
        score += abs(10-gcCount)*gcWeight

	for bp, position, wt in parameters:
                subSeq = seq[position:position+len(bp)]
                if subSeq==bp:
                        score += wt
        return (1.0/(1.0 + math.exp(-score)))

def LINDELNGcalculatescore(seq):
	parameters = [('AA',7,-0.166151119),('AA',15,-0.116635449),('AA',16,-0.141978972),('AA',20,-0.46209483),('AA',21,-1.159033893),('AA',23,-0.14092585),('AA',28,-0.140710973),('TA',7,0.242102512),('TA',10,-0.051811653),('TA',16,0.462935638),('TA',17,0.587570687),
('TA',21,-0.321626833),('TA',22,0.221199793),('TA',26,-0.308936349),('TA',27,0.228218948),('GA',0,-0.213575878),('GA',6,0.155415387),('GA',11,-0.066555092),('GA',12,0.155827764),('GA',13,0.191500133),('GA',14,0.105536077),('GA',15,0.104278628),
('GA',17,-0.134284261),('GA',18,0.144061168),('GA',25,0.111163684),('CA',2,0.038959214),('CA',7,0.126214472),('CA',8,-0.150300426),('CA',17,0.297525515),('CA',19,0.25467499),('CA',20,-0.293330178),('CA',22,0.238690067),('CA',23,-0.21982359),('CA',27,2.054117089),('CA',28,0.335055713),('AT',1,-0.097386892),('AT',2,0.123999347),('AT',8,0.081022603),('AT',10,-0.108078392),('AT',12,-0.079018399),('AT',18,-0.174383894),('AT',19,0.320177513),('TT',6,-0.221409382),('TT',7,-0.157901062),('TT',8,-0.061772525),('TT',10,-0.150198902),('TT',11,-0.448551556),('TT',12,-0.622624981),('TT',13,-0.139416659),('TT',14,-0.077523842),('TT',15,-0.948417398),('TT',16,-0.080861838),('TT',17,-1.046258319),('TT',18,-0.375026297),('TT',20,-1.119104273),('TT',21,-0.166128832),('TT',22,-0.930378086),('GT',21,0.510522765),('GT',22,0.332401486),('GT',23,0.219019899),('GT',25,-0.154193855),('CT',10,0.202473712),('CT',14,-0.158357686),('CT',15,-0.066814565),('CT',16,0.160473088),('CT',22,-0.552932155),('AG',13,0.160018261),('AG',19,-0.108198604),('AG',20,0.144071434),('AG',22,0.25226124),('AG',23,1.063887668),('AG',24,-0.009612897),('AG',28,0.122018273),('TG',5,0.241363252),('TG',9,0.079465841),('TG',10,0.149806008),('TG',12,0.049946857),('TG',13,0.155486915),('TG',14,0.031926277),('TG',18,-0.277324377),('TG',20,0.691550937),('TG',26,-0.08254848),('GG',1,-0.157621236),('GG',6,0.278519061),('GG',9,0.259115271),('GG',23,2.055614334),('CG',2,-0.436407352),('CG',7,-0.046536109),('CG',16,-0.15453076),('CG',17,-0.205884833),('CG',18,0.191552793),('CG',20,-0.245249633),
('CG',21,-0.542330163),('CG',27,-0.185411057),('AC',0,-0.085648062),('AC',6,0.148699235),('AC',18,0.105302637),('AC',20,0.175012001),('AC',23,0.264370297),('AC',26,-0.184483722),('TC',0,0.112147888),('TC',11,0.217316946),('TC',13,-0.310569559),
('TC',14,-0.121534773),('TC',21,-0.247459342),('TC',22,-1.241334903),('TC',23,-0.628940116),('TC',27,0.536660199),('TC',28,-0.093641005),('GC',1,0.276788432),('GC',4,-0.00869246),('GC',9,0.128560922),('GC',10,0.080325072),('GC',12,-0.219567722),
('GC',18,-0.073212082),('GC',19,-0.210102633),('GC',20,-0.457465406),('GC',21,-0.17993102),('GC',22,-1.458409354),('CC',6,-0.110326586),('CC',9,-0.259044689),('CC',16,0.38220164),('CC',17,0.13474772),('CC',18,0.352748051),('CC',19,0.771334455),('CC',21,0.234536174),('CC',26,0.757158333),('CC',28,-0.136943015),('A',10,-0.04673173),('A',20,0.05239417),('A',24,-0.22018074),('A',26,3.09091971),('T',19,-0.36533529),('T',21,-0.35869458),('T',23,-0.35126663),('T',26,-0.8356666),('G',6,0.05688286),('G',17,-0.10189948),('G',20,-0.05641955),('G',23,0.22952347),('G',24,0.90850369),('G',29,0.07510093),('C',5,-0.11228712),('C',17,0.10768524),('C',18,0.05610607),('C',19,0.05474233),('C',20,-0.02714258),('C',21,0.13026953),('C',23,-0.27420231)]

	intercept = -4.472148961
	score = intercept
	Entropy = 1.459539385
	Free_energy = 0.17806217
	GChigh = 0.012356139
	GClow = -1.05E-13
	AG = 0.09453565
	TG = 0.034280035
	TT = -0.179103124
	
	#### feature extraction 
        ##Free_energy
        Energy = seq[4:24]
        Energycal = RNA.fold(Energy)[-1]
        Energycal = round(Energycal,0)
        score +=(Energycal*Free_energy)
        # independent nucleotide composition
        score = score + AG*(seq.count('AG'))
        score = score + TT*(seq.count('TT'))
        score = score + TG*(seq.count('TA'))
	#entropy cal
	entropy = dict()
        Entropycal = []
        entropy_seq = seq[4:24]
        lentseq= len(entropy_seq)
	ntscount = {'A':0, 'G':0, 'T':0, 'C':0}
	nuc = ['A','T','G','C']
        for ant in nuc:
                ntscount[ant]=(entropy_seq.count(ant))/float((lentseq))
        for ant in nuc:
                if ntscount[ant]!=0:
                        entropy[ant] = -(ntscount[ant]*np.log2(ntscount[ant]))
                else:
                        entropy[ant] = 0
	entropySum = sum(entropy.values())
	entropySumR = round(entropySum,1)
	score += (entropySumR*Entropy)
	
	guideSeq = seq[4:24]
        gcCount = guideSeq.count("G") + guideSeq.count("C")
        if gcCount <= 10:
                gcWeight = GClow
        if gcCount > 10:
                gcWeight = GChigh
        score += abs(10-gcCount)*gcWeight

	for bp, position, wt in parameters:
                subSeq = seq[position:position+len(bp)]
                if subSeq==bp:
                        score += wt
        return (1.0/(1.0 + math.exp(-score)))

def Comprehensive(fastfile):
        cas9 = {}
	cas9i = {}
	cas9a = {}
	cas9NGA = {}
	cas9NGC = {}
	cas9NGT = {}
	cas12a = {}
	ncas9 = {}
	ncas9i = {}
	ncas9a = {}
	ncas9NGA = {}
	ncas9NGC = {}
	ncas9NGT = {}
	ncas12a = {}
	sara = {}
	ast = []
	comb = {}
	ncomb = {}
	can = {}
	ncan = {}
	eff = {}
	nff = {}
	sff = {}
	rff = {}
	file1 = open(fastfile,'r')
	file1  = file1.readlines()
	cas9f = open('CGD.txt','wb')
	cas9f.write('ID' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Strand' + '\t' + 'Sequence' + '\t' + 'CGDi' + '\t' + 'CGDa' + '\t' + 'CGD9' + '\t' + 'CGDNGR' + '\t' + 'CGD12a' +'\n')
	for line in file1:
		if line.startswith(">"):
			newline = line.strip('\r\n')
			nare = newline[1:]
			nare = nare.split(' ')[0]
			ps = nare
			sara[ps] = []
		else:
			newsequence = line.strip('\r\n').split('\n')
			for az in newsequence:
				ast.append(az)
			my_str = ''.join(map(str,ast))		
			sara[ps] = [my_str]
	for key,value in sara.items():
		my_value = value[0]
		if (int(100) <= len(my_value)<= 10000):
			for n in xrange(len(my_value)): 
				st = my_value.find('GG',n)
				if st == n:
					if int(25) < int(st) < len(my_value)-int(5):
						nstart = int(st) - int(25)
						nend = int(st) + int(5)
						sequences = my_value[nstart:nend]
						ids = 'gRNAc_' + str(nstart)
						idi = 'gRNAi_' + str(nstart)
						ida = 'gRNAa_' + str(nstart)
						idcoms = 'gRNAcoms_' + str(nstart)
						idfs = 'gRNAidfs_' + str(nstart)
						IDs = key
						score_LINDELa = round(LINDELacalculatescore(sequences),2)
						score_LINDELi = round(LINDELicalculatescore(sequences),2)
						score_LINDEL = round(LINDELcalculatescore(sequences),2)
                                                cas9i[idi] = [nstart,nend,'+',sequences,float(score_LINDELi),float(score_LINDELa),float(score_LINDEL),float(0.0),float(0.0),str(IDs)]
                                                
                                                cas9a[ida] = [nstart,nend,'+',sequences,float(score_LINDELi),float(score_LINDELa),float(score_LINDEL),float(0.0),float(0.0),str(IDs)]
                                               
                                                cas9[ids] = [nstart,nend,'+',sequences,float(score_LINDELi),float(score_LINDELa),float(score_LINDEL),float(0.0),float(0.0),str(IDs)]
			


			for t in xrange(len(my_value)): 
				st = my_value.find('CC',t)
				if st == t:
					if int(3) < int(st) < len(my_value)-int(27):
                                		nstart = int(st) - int(3) 
                                        	nend = int(st) + int(27)
                                        	sequences = my_value[nstart:nend]
						nsequences = reverseComp(sequences)
                                        	idss = 'gRNAcn_' + str(nstart)
                                        	idii = 'gRNAin_' + str(nstart)
                                        	idaa = 'gRNAan_' + str(nstart)
						id12a = 'gRNA12a_' + str(nstart)
						idcom = 'gRNAcom_' + str(nstart)
						idef = 'gRNnidfs_' + str(nstart)
						IDs = key
                                        	score_LINDELa = round(LINDELacalculatescore(nsequences),2)
                                        	score_LINDELi = round(LINDELicalculatescore(nsequences),2)
                                        	score_LINDEL = round(LINDELcalculatescore(nsequences),2)
						
						ncas9i[idii] = [nstart,nend,'-',nsequences,float(score_LINDELi),float(score_LINDELa),float(score_LINDEL),float(0.0),float(0.0),str(IDs)]
						
						ncas9a[idaa] = [nstart,nend,'-',nsequences,float(score_LINDELi),float(score_LINDELa),float(score_LINDEL),float(0.0),float(0.0),str(IDs)]
						
						ncas9[idss] = [nstart,nend,'-',nsequences,float(score_LINDELi),float(score_LINDELa),float(score_LINDEL),float(0.0),float(0.0),str(IDs)]


			for n in xrange(len(my_value)): 
				st = my_value.find('GA',n)
				if st == n:
					if int(25) < int(st) < len(my_value)-int(5):
						nstart = int(st) - int(25)
						nend = int(st) + int(5)
						sequences = my_value[nstart:nend]
						idga = 'gRNAga_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(sequences),2)
                                                cas9NGA[idga] = [nstart,nend,'+',sequences,float(0.0),float(0.0),float(0.0),float(score_LINDELNG),float(0.0),str(IDs)]

			for n in xrange(len(my_value)): 
				st = my_value.find('TC',n)
				if st == n:
					if int(3) < int(st) < len(my_value)-int(27):
						nstart = int(st) - int(3)
						nend = int(st) + int(27)
						sequences = my_value[nstart:nend]
						revseq = reverseComp(sequences)
						idct = 'gRNAct_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(revseq),2)
                                                ncas9NGA[idct] = [nstart,nend,'-',sequences,float(0.0),float(0.0),float(0.0),float(score_LINDELNG),float(0.0),str(IDs)]
		
			

			
			for n in xrange(len(my_value)): 
				st = my_value.find('GC',n)
				if st == n:
					if int(25) < int(st) < len(my_value)-int(5):
						nstart = int(st) - int(25)
						nend = int(st) + int(5)
						sequences = my_value[nstart:nend]
						idgc = 'gRNAgc_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(sequences),2)
                                                cas9NGC[idgc] = [nstart,nend,'+',sequences,float(0.0),float(0.0),float(0.0),float(score_LINDELNG),float(0.0),str(IDs)]

			for n in xrange(len(my_value)): 
				st = my_value.find('GC',n)
				if st == n:
					if int(3) < int(st) < len(my_value)-int(27):
						nstart = int(st) - int(3)
						nend = int(st) + int(27)
						sequences = my_value[nstart:nend]
						revseq = reverseComp(sequences)
						idcg = 'gRNAcg_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(revseq),2)
                                                ncas9NGC[idcg] = [nstart,nend,'-',sequences,float(0.0),float(0.0),float(0.0),float(score_LINDELNG),float(0.0),str(IDs)]

			for n in xrange(len(my_value)): 
				st = my_value.find('GT',n)
				if st == n:
					if int(25) < int(st) < len(my_value)-int(5):
						nstart = int(st) - int(25)
						nend = int(st) + int(5)
						sequences = my_value[nstart:nend]
						idgt = 'gRNAgt_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(sequences),2)
                                                cas9NGT[idgt] = [nstart,nend,'+',sequences,float(0.0),float(0.0),float(0.0),float(score_LINDELNG),float(0.0),str(IDs)]

			for n in xrange(len(my_value)): 
				st = my_value.find('AC',n)
				if st == n:
					if int(3) < int(st) < len(my_value)-int(27):
						nstart = int(st) - int(3)
						nend = int(st) + int(27)
						sequences = my_value[nstart:nend]
						revseq = reverseComp(sequences)
						idnca = 'gRNAnca_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(revseq),2)
                                                ncas9NGT[idnca] = [nstart,nend,'-',sequences,float(0.0),float(0.0),float(0.0),float(score_LINDELNG),float(0.0),str(IDs)]
	
												
					
			for n in xrange(len(my_value)):
				st = my_value.find('TTT',n)
				if st == n:
					if int(4) < int(st) < len(my_value)-int(30):
                                        	nstart = int(st) - int(4)
                                                nend = int(st) + int(30)
                                                sequences = my_value[nstart:nend]
						idscas = 'gRNAcas_' + str(nstart)
						score_CINDEL = round(CINDELcalculatescore(sequences),2)
						cas12a[idscas] = [nstart,nend,'+',sequences,float(0.0),float(0.0),float(0.0),float(0.0),float(score_CINDEL),str(IDs)]

			for p in xrange(len(my_value)):
                                st = my_value.find('AAA',p)
                                if st == p:
                                        if int(27) < int(st) < len(my_value)-int(7):
                                                nstart = int(st) - int(27)
                                                nend = int(st) + int(7)
                                                sequences = my_value[nstart:nend]
						nsequences = reverseComp(sequences)
                                                idsca = 'gRNAcas_' + str(nstart)
                                                score_CINDEL = round(CINDELcalculatescore(sequences),2)
                                                ncas12a[idsca] = [nstart,nend,'+',sequences,float(0.0),float(0.0),float(0.0),float(0.0),float(score_CINDEL),str(IDs)]
					
		
		else:
			print "length not appropriate"

	aspha = dict(chain(cas9i.items(),ncas9i.items(),cas9a.items(),ncas9a.items(),cas9.items(),ncas9.items(),cas9NGA.items(),ncas9NGA.items(),cas9NGC.items(),ncas9NGC.items(),cas9NGT.items(),ncas9NGT.items(),cas12a.items(),ncas12a.items()))
	dyt = OrderedDict(sorted(aspha.items(), key=lambda x:x[1][6],reverse=True))
	for key,value in dyt.items():
		cas9f.write(str(value[9]) + '\t' + str(value[0]) + '\t' + str(value[1]) + '\t' + str(value[2]) + '\t' + str(value[3]) + '\t'+ str(value[4]) +  '\t' + str(value[5]) + '\t' + str(value[6]) + '\t' + str(value[7]) + '\t' + str(value[8]) +'\n')
	cas9f.close()
	
		

def CGDi(fastfile):
	cas9 = {}
        cas9i = {}
        cas9a = {}
        cas12a = {}
        ncas9 = {}
        ncas9i = {}
        ncas9a = {}
        ncas12a = {}
        sara = {}
        ast = []
        comb = {}
        ncomb = {}
        can = {}
        ncan = {}
	ncas9i = {}
	cas9i = {}
	file1 = open(fastfile,'r')
        file1  = file1.readlines()
	cas9if = open('CGDi.txt','wb')
	cas9if.write('ID' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Strand' + '\t' + 'Sequence' + '\t' + 'CGDi'  + '\n')
	for line in file1:
                if line.startswith(">"):
                        newline = line.strip('\r\n')
                        nare = newline[1:]
                        nare = nare.split(' ')[0]
                        ps = nare[:-2]
                        sara[ps] = []
                else:
                        newsequence = line.strip('\r\n').split('\n')
                        for az in newsequence:
                                ast.append(az)
                        my_str = ''.join(map(str,ast))
                        sara[ps] = [my_str]
        for key,value in sara.items():
                my_value = value[0]
                if (int(100) <= len(my_value) <= int(10000)):
                        for n in xrange(len(my_value)):
                                st = my_value.find('GG',n)
                                if st == n:
                                        if int(25) < int(st) < len(my_value)-int(5):
                                                nstart = int(st) - int(25)
                                                nend = int(st) + int(5)
                                                sequences = my_value[nstart:nend]
                                                idn = 'gRNAi_' + str(nstart)
						ids = key
						score_LINDELi = round(LINDELicalculatescore(sequences),2)
						cas9i[idn] = [nstart,nend,'+',sequences,float(score_LINDELi),str(ids)]

			for p in xrange(len(my_value)):
                                st = my_value.find('CC',p)
                                if st == p:
                                        if int(3) < int(st) < len(my_value)-int(27):
                                                nstart = int(st) - int(3)
                                                nend = int(st) + int(27)
                                                sequences = my_value[nstart:nend]
						nsequences = reverseComp(sequences)
                                                idi = 'gRNAn_' + str(nstart)
						ids = key
                                                score_LINDELi = round(LINDELicalculatescore(nsequences),2)
                                                
                                                ncas9i[idi] = [nstart,nend,'-',nsequences,float(score_LINDELi),str(ids)]

	newcasi = dict(chain(cas9i.items(),ncas9i.items()))
	syt = OrderedDict(sorted(newcasi.items(), key=lambda x:x[1][4],reverse=True))
	for key,value in syt.items():
                cas9if.write(str(value[5]) + '\t' + str(value[0]) + '\t' + str(value[1]) + '\t' + str(value[2]) + '\t' + str(value[3]) + '\t' + str(value[4]) + '\n')
        cas9if.close()
	return cas9if

def CGDa(fastfile):
	cas9 = {}
        cas9i = {}
        cas9a = {}
        cas12a = {}
        ncas9 = {}
        ncas9i = {}
        ncas9a = {}
        ncas12a = {}
        sara = {}
        ast = []
        comb = {}
        ncomb = {}
        can = {}
        ncan = {}
        ncas9i = {}
        cas9i = {}
	file1 = open(fastfile,'r')
        file1  = file1.readlines()
        cas9if = open('CGDa.txt','wb')
        cas9if.write('ID' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Strand' + '\t' + 'Sequence' + '\t' + 'CGDa'  + '\n')
        for line in file1:
                if line.startswith(">"):
                        newline = line.strip('\r\n')
                        nare = newline[1:]
                        nare = nare.split(' ')[0]
                        ps = nare[:-2]
                        sara[ps] = []
                else:
                        newsequence = line.strip('\r\n').split('\n')
                        for az in newsequence:
                                ast.append(az)
                        my_str = ''.join(map(str,ast))
                        sara[ps] = [my_str]
        for key,value in sara.items():
                my_value = value[0]
                if (int(100) <= len(my_value) <= int(10000)):
                        for n in xrange(len(my_value)):
                                st = my_value.find('GG',n)
                                if st == n:
                                        if int(25) < int(st) < len(my_value)-int(5):
                                                nstart = int(st) - int(25)
                                                nend = int(st) + int(5)
                                                sequences = my_value[nstart:nend]
                                                idn = 'gRNAi_' + str(nstart)
						ids = key
                                                score_LINDELa = round(LINDELacalculatescore(sequences),2)
                                                cas9i[idn] = [nstart,nend,'+',sequences,float(score_LINDELa),str(ids)]

                        for p in xrange(len(my_value)):
                                st = my_value.find('CC',p)
                                if st == p:
                                        if int(3) < int(st) < len(my_value)-int(27):
                                                nstart = int(st) - int(3)
                                                nend = int(st) + int(27)
                                                sequences = my_value[nstart:nend]
                                                nsequences = reverseComp(sequences)
                                                idi = 'gRNAn_' + str(nstart)
                                                score_LINDELa = round(LINDELacalculatescore(nsequences),2)
                                                ncas9i[idi] = [nstart,nend,'-',nsequences,float(score_LINDELa),str(ids)]

        newcasi = dict(chain(cas9i.items(),ncas9i.items()))
        syt = OrderedDict(sorted(newcasi.items(), key=lambda x:x[1][4],reverse=True))
        for key,value in syt.items():
                cas9if.write(str(value[5]) + '\t' + str(value[0]) + '\t' + str(value[1]) + '\t' + str(value[2]) + '\t' + str(value[3]) + '\t' + str(value[4]) + '\n')
        cas9if.close()
	return cas9if

def CGD9(fastfile):
	cas9 = {}
        cas9i = {}
        cas9a = {}
        cas12a = {}
        ncas9 = {}
        ncas9i = {}
        ncas9a = {}
        ncas12a = {}
        sara = {}
        ast = []
        comb = {}
        ncomb = {}
        can = {}
        ncan = {}
        ncas9i = {}
        cas9i = {}
	file1 = open(fastfile,'r')
        file1  = file1.readlines()
        cas9if = open('CGD9.txt','wb')
        cas9if.write('ID' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Strand' + '\t' + 'Sequence' + '\t' + 'CGD9'  + '\n')
        for line in file1:
                if line.startswith(">"):
                        newline = line.strip('\r\n')
                        nare = newline[1:]
                        nare = nare.split(' ')[0]
                        ps = nare[:-2]
                        sara[ps] = []
                else:
                        newsequence = line.strip('\r\n').split('\n')
                        for az in newsequence:
                                ast.append(az)
                        my_str = ''.join(map(str,ast))
                        sara[ps] = [my_str]
        for key,value in sara.items():
                my_value = value[0]
                if (int(100) <= len(my_value) <= int(10000)):
                        for n in xrange(len(my_value)):
                                st = my_value.find('GG',n)
                                if st == n:
                                        if int(25) < int(st) < len(my_value)-int(5):
                                                nstart = int(st) - int(25)
                                                nend = int(st) + int(5)
                                                sequences = my_value[nstart:nend]
                                                idn = 'gRNAi_' + str(nstart)
                                                score_LINDEL = round(LINDELcalculatescore(sequences),2)
                                                cas9i[idn] = [nstart,nend,'+',sequences,float(score_LINDEL),str(key)]
                        for p in xrange(len(my_value)):
                                st = my_value.find('CC',p)
                                if st == p:
                                        if int(3) < int(st) < len(my_value)-int(27):
                                                nstart = int(st) - int(3)
                                                nend = int(st) + int(27)
                                                sequences = my_value[nstart:nend]
                                                nsequences = reverseComp(sequences)
                                                idi = 'gRNAn_' + str(nstart)
                                                score_LINDEL = round(LINDELcalculatescore(nsequences),2)
                                                ncas9i[idi] = [nstart,nend,'-',nsequences,float(score_LINDEL),str(key)]

        newcasi = dict(chain(cas9i.items(),ncas9i.items()))
        syt = OrderedDict(sorted(newcasi.items(), key=lambda x:x[1][4],reverse=True))
        for key,value in syt.items():
                cas9if.write(str(value[5]) + '\t' + str(value[0]) + '\t' + str(value[1]) + '\t' + str(value[2]) + '\t' + str(value[3]) + '\t' + str(value[4]) + '\n')
        cas9if.close()
	return cas9if

def CGD12a(fastfile):
	cas9 = {}
        cas9i = {}
        cas9a = {}
        cas12a = {}
        ncas9 = {}
        ncas9i = {}
        ncas9a = {}
        ncas12a = {}
        sara = {}
        ast = []
        comb = {}
        ncomb = {}
        can = {}
        ncan = {}
        ncas9i = {}
        cas9i = {}	
	ncas9i = {}
        cas9i = {}
	file1 = open(fastfile,'r')
        file1  = file1.readlines()
        cas9if = open('CGD12a.txt','wb')
        cas9if.write('ID' + 'Start' + '\t' + 'End' + '\t' + 'Strand' + '\t' + 'Sequence' + '\t' + 'CGD12a'  + '\n')
        for line in file1:
                if line.startswith(">"):
                        newline = line.strip('\r\n')
                        nare = newline[1:]
                        nare = nare.split(' ')[0]
                        ps = nare[:-2]
                        sara[ps] = []
                else:
                        newsequence = line.strip('\r\n').split('\n')
                        for az in newsequence:
                                ast.append(az)
                        my_str = ''.join(map(str,ast))
                        sara[ps] = [my_str]
        for key,value in sara.items():
                my_value = value[0]
                if (int(100) <= len(my_value) <= int(10000)):
                        for n in xrange(len(my_value)):
                                st = my_value.find('TTT',n)
                                if st == n:
                                        if int(4) < int(st) < len(my_value)-int(30):
                                                nstart = int(st) - int(4)
                                                nend = int(st) + int(30)
                                                sequences = my_value[nstart:nend]
                                                idn = 'gRNAi_' + str(nstart)
                                                score_CINDEL = round(CINDELcalculatescore(sequences),2)
                                                cas9i[idn] = [nstart,nend,'+',sequences,float(score_CINDEL),str(key)]

                        for p in xrange(len(my_value)):
                                st = my_value.find('AAA',p)
                                if st == p:
                                        if int(27) < int(st) < len(my_value)-int(7):
                                                nstart = int(st) - int(27)
                                                nend = int(st) + int(7)
                                                sequences = my_value[nstart:nend]
                                                nsequences = reverseComp(sequences)
                                                idi = 'gRNAn_' + str(nstart)
                                                score_CINDEL = round(CINDELcalculatescore(nsequences),2)
                                           	ncas9i[idi] = [nstart,nend,'-',nsequences,float(score_CINDEL),str(key)]

        newcasi = dict(chain(cas9i.items(),ncas9i.items()))
        syt = OrderedDict(sorted(newcasi.items(), key=lambda x:x[1][4],reverse=True))
        for key,value in syt.items():
                cas9if.write(str(value[5]) + '\t' + str(value[0]) + '\t' + str(value[1]) + '\t' + str(value[2]) + '\t' + str(value[3]) + '\t' + str(value[4]) + '\n')
        cas9if.close()
	return cas9if

def CGD9NG(fastfile):
	cas9 = {}
	cas9i = {}
	cas9a = {}
	cas9NGA = {}
	cas9NGC = {}
	cas9NGT = {}
	cas12a = {}
	ncas9 = {}
	ncas9i = {}
	ncas9a = {}
	ncas9NGA = {}
	ncas9NGC = {}
	ncas9NGT = {}
	ncas12a = {}
	sara = {}
	ast = []
	comb = {}
	ncomb = {}
	can = {}
	ncan = {}
	eff = {}
	nff = {}
	sff = {}
	rff = {}
	file1 = open(fastfile,'r')
        file1  = file1.readlines()
        cas9ng = open('CGD9NG.txt','wb')
        cas9ng.write('ID' + 'Start' + '\t' + 'End' + '\t' + 'Strand' + '\t' + 'Sequence' + '\t' + 'CGD9NG'  + '\n')
	for line in file1:
                if line.startswith(">"):
                        newline = line.strip('\r\n')
                        nare = newline[1:]
                        nare = nare.split(' ')[0]
                        ps = nare[:-2]
                        sara[ps] = []
                else:
                        newsequence = line.strip('\r\n').split('\n')
                        for az in newsequence:
                                ast.append(az)
                        my_str = ''.join(map(str,ast))
                        sara[ps] = [my_str]
        for key,value in sara.items():
                my_value = value[0]
                if (int(100) <= len(my_value) <= int(10000)):
			for n in xrange(len(my_value)): 
				st = my_value.find('GA',n)
				if st == n:
					if int(25) < int(st) < len(my_value)-int(5):
						nstart = int(st) - int(25)
						nend = int(st) + int(5)
						sequences = my_value[nstart:nend]
						idga = 'gRNAga_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(sequences),2)
                                                cas9NGA[idga] = [nstart,nend,'+',sequences,float(score_LINDELNG),str(key)]

			for n in xrange(len(my_value)): 
				st = my_value.find('TC',n)
				if st == n:
					if int(3) < int(st) < len(my_value)-int(27):
						nstart = int(st) - int(3)
						nend = int(st) + int(27)
						sequences = my_value[nstart:nend]
						revseq = reverseComp(sequences)
						idct = 'gRNAct_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(revseq),2)
                                                ncas9NGA[idct] = [nstart,nend,'-',revseq,float(score_LINDELNG),str(key)]
		
			

			
			for n in xrange(len(my_value)): 
				st = my_value.find('GC',n)
				if st == n:
					if int(25) < int(st) < len(my_value)-int(5):
						nstart = int(st) - int(25)
						nend = int(st) + int(5)
						sequences = my_value[nstart:nend]
						idgc = 'gRNAgc_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(sequences),2)
                                                cas9NGC[idgc] = [nstart,nend,'+',sequences,float(score_LINDELNG),str(key)]

			for n in xrange(len(my_value)): 
				st = my_value.find('GC',n)
				if st == n:
					if int(3) < int(st) < len(my_value)-int(27):
						nstart = int(st) - int(3)
						nend = int(st) + int(27)
						sequences = my_value[nstart:nend]
						revseq = reverseComp(sequences)
						idcg = 'gRNAcg_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(revseq),2)
                                                ncas9NGC[idcg] = [nstart,nend,'-',revseq,float(score_LINDELNG),str(key)]

			for n in xrange(len(my_value)): 
				st = my_value.find('GT',n)
				if st == n:
					if int(25) < int(st) < len(my_value)-int(5):
						nstart = int(st) - int(25)
						nend = int(st) + int(5)
						sequences = my_value[nstart:nend]
						idgt = 'gRNAgt_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(sequences),2)
                                                cas9NGT[idgt] = [nstart,nend,'+',sequences,float(score_LINDELNG),str(key)]

			for n in xrange(len(my_value)): 
				st = my_value.find('AC',n)
				if st == n:
					if int(3) < int(st) < len(my_value)-int(27):
						nstart = int(st) - int(3)
						nend = int(st) + int(27)
						sequences = my_value[nstart:nend]
						revseq = reverseComp(sequences)
						idnca = 'gRNAnca_' + str(nstart)
						score_LINDELNG = round(LINDELNGcalculatescore(revseq),2)
                                                ncas9NGT[idnca] = [nstart,nend,'-',revseq,float(score_LINDELNG),str(key)]
	newcasi = dict(chain(cas9NGA.items(),ncas9NGA.items(),cas9NGC.items(),ncas9NGC.items(),cas9NGT.items(),ncas9NGT.items()))
        syt = OrderedDict(sorted(newcasi.items(), key=lambda x:x[1][4],reverse=True))
        for key,value in syt.items():
                cas9ng.write(str(value[5]) + '\t' + str(value[0]) + '\t' + str(value[1]) + '\t' + str(value[2]) + '\t' + str(value[3]) + '\t' + str(value[4]) + '\n')
        cas9ng.close()
	return cas9ng	

def helps(ane):

	parser = argparse.ArgumentParser(description = 'working with CGD')
	parser.add_argument('-a',metavar='\b',type=str,help = 'The option generates Comprehensive score using Comprehesive function')
	parser.add_argument('-b',metavar='\b',help = 'The option generates CRISPRi score using CGDi function')
	parser.add_argument('-c',metavar='\b',help = 'The option generates CRISPRa score using CGDa function')
	parser.add_argument('-d',metavar='\b',help = 'The option generates CRISPR-Cas9 score using CGD9 function')
	parser.add_argument('-e',metavar='\b',help = 'The option generates CRISPR-Cas12a score using CGD12a function')
	parser.add_argument('-f',metavar='\b',help = 'The option generates CRISPR-Cas9 non canonical score using CGD9NG function')
	args = parser.parse_args()
	return args


	
func_arg = {"-a": Comprehensive, "-b": CGDi, "-c": CGDa, "-d": CGD9, "-e": CGD12a, "-f": CGD9NG, "-h": helps}

if __name__ == "__main__":
	func_arg[sys.argv[1]](sys.argv[2])

