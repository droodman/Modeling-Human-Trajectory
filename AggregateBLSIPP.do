// Adaptation for Stata/Mata of the Matlab program AggregateBLSIPP.m, in the Bloom et al. (2020) data+code archive
// openicpsr.org/openicpsr/project/111743/version/V2/view?path=/openicpsr/111743/fcr:versions/V2/Aggregate/AggregateBLSIPP.m&type=file
// This version does not generate graphs. It incorporates directly in the code the contents of the Matlab data file WageScientistData.mat.
// It also inserts a CPI-U data series in the WageScientistData data matrix, so that that may be used instead of the high-skilled wage index
// for deflating research investment.
// To undo that change and mimic the original, replace the 2 in line 491 with a 3
// Estimate of beta, of primary interest here, is printed in line 540


// AggregateBLSIPP.m  6/17/2016
//
//   USES Total IPP instead of just R\&D
//
//  Basic Idea TFP calculation using BLS's decadal TFP growth numbers
//  and the NIPA IPP numbers
//
//  NIPA data on R&D and Intellectual Property Products, via FRED (RND data list)



cap mata mata drop delchange()
cap mata mata drop trimr()

mata

real matrix delchange(real matrix X)
  return (X[|2,.\.,.|] - X[|.,.\rows(X)-1,.|])
real matrix trimr(real matrix X, real scalar top, real scalar bot)
  return (X[|1+top,.\rows(X)-bot,.|])

Lambda = 1

Gordondata=
// From Bob Gordon Figure 16-5. Currently a rough guess
//  See GordonCh16tables_151226.xls, "New TFP", cells T154 and following.
// 	New
//Decade	Annual TFP Growth
//1900	0.3747
//1910	0.2516
//1920	0.7521
//1930	1.2460 // The 1930 observation is for the 1920s
1940,	1.8155 \
1950,	3.3888 \
1960,	1.6026 \
1970,	1.4009 \
1980,	0.3445 \
1990,	0.7828 \
2000,	0.7652 \
2014,	0.6800


"Gordon Data"
Gordondata


GordonGrowth=Gordondata[,2]
GordonYears=Gordondata[,1]
DecadeNames=tokens("1930s 1940s 1950s 1960s 1970s 1980s 1990s 2000s")


// Y001RC1A027NBEA
// FRED Graph Observations		FRED Graph Observations			
// Federal Reserve Economic Data		Federal Reserve Economic Data			
// Link: https://fred.stlouisfed.org		Link: https://fred.stlouisfed.org			
// Help: https://fred.stlouisfed.org/help-faq		Help: https://fred.stlouisfed.org/help-faq			
// Economic Research Division		Economic Research Division			
// Federal Reserve Bank of St. Louis		Federal Reserve Bank of St. Louis			
// Gross Private Domestic Investment: Fixed Investment: Nonresidential: Intellectual Property Products, Billions of Dollars, Annual, Not Seasonally Adjusted					
// Y001RC1A027NBEA		Y055RC1A027NBEA	Government Gross Investment: Intellectual Property Products, Billions of Dollars, Annual, Not Seasonally Adjusted		
					
// Frequency: Annual		Frequency: Annual			
// observation_date	Y001RC1A027NBEA	observation_date	Y055RC1A027NBEA		TOTAL GPDI in IPP

datarnd=  //  PvtIPP                 GovtIPP        TotalIPP
1929-01-01,	0.6,  1929-01-01,	0.1	,	0.7 \
1930-01-01,	0.6,  1930-01-01,	0.1	,	0.7 \
1931-01-01,	0.5,  1931-01-01,	0.1	,	0.6 \
1932-01-01,	0.4,  1932-01-01,	0.1	,	0.5 \
1933-01-01,	0.4,  1933-01-01,	0.1	,	0.5 \
1934-01-01,	0.5,  1934-01-01,	0.1	,	0.6 \
1935-01-01,	0.6,  1935-01-01,	0.1	,	0.7 \
1936-01-01,	0.6,  1936-01-01,	0.1	,	0.7 \
1937-01-01,	0.7,  1937-01-01,	0.1	,	0.8 \
1938-01-01,	0.8,  1938-01-01,	0.1	,	0.9 \
1939-01-01,	0.8,  1939-01-01,	0.1	,	0.9 \
1940-01-01,	0.8,  1940-01-01,	0.1	,	0.9 \
1941-01-01,	1.1,  1941-01-01,	0.3	,	1.4 \
1942-01-01,	1.2,  1942-01-01,	0.5	,	1.7 \
1943-01-01,	1.1,  1943-01-01,	0.9	,	2.0 \
1944-01-01,	1.2,  1944-01-01,	1.6	,	2.8 \
1945-01-01,	1.4,  1945-01-01,	1.5	,	2.9 \
1946-01-01,	1.8,  1946-01-01,	1.4	,	3.2 \
1947-01-01,	2.0,  1947-01-01,	1.4	,	3.4 \
1948-01-01,	2.1,  1948-01-01,	1.6	,	3.7 \
1949-01-01,	2.0,  1949-01-01,	1.7	,	3.7 \
1950-01-01,	2.3,  1950-01-01,	1.9	,	4.2 \
1951-01-01,	2.4,  1951-01-01,	2.1	,	4.5 \
1952-01-01,	3.0,  1952-01-01,	2.4	,	5.4 \
1953-01-01,	3.7,  1953-01-01,	2.7	,	6.4 \
1954-01-01,	3.9,  1954-01-01,	3.0	,	6.9 \
1955-01-01,	4.3,  1955-01-01,	3.6	,	7.9
datarnd = datarnd \
1956-01-01,	5.2,  1956-01-01,	4.8	,	10.0 \
1957-01-01,	5.6,  1957-01-01,	5.9	,	11.5 \
1958-01-01,	6.0,  1958-01-01,	6.4	,	12.4 \
1959-01-01,	6.6,  1959-01-01,	7.0	,	13.6 \
1960-01-01,	7.1,  1960-01-01,	7.8	,	14.9 \
1961-01-01,	8.0,  1961-01-01,	8.8	,	16.8 \
1962-01-01,	8.4,  1962-01-01,	9.9	,	18.3 \
1963-01-01,	9.2,  1963-01-01,	11.7,		20.9 \
1964-01-01,	9.8,  1964-01-01,	12.9,		22.7 \
1965-01-01,	11.1, 1965-01-01,	13.8,		24.9 \
1966-01-01,	12.8, 1966-01-01,	15.4,		28.2 \
1967-01-01,	14.0, 1967-01-01,	16.2,		30.2 \
1968-01-01,	15.6, 1968-01-01,	17.0,		32.6 \
1969-01-01,	17.2, 1969-01-01,	17.7,		34.9 \
1970-01-01,	17.9, 1970-01-01,	17.6,		35.5 \
1971-01-01,	18.7, 1971-01-01,	18.1,		36.8 \
1972-01-01,	20.6, 1972-01-01,	19.3,		39.9 \
1973-01-01,	22.7, 1973-01-01,	20.3,		43.0 \
1974-01-01,	25.5, 1974-01-01,	21.5,		47.0 \
1975-01-01,	27.8, 1975-01-01,	23.3,		51.1 \
1976-01-01,	32.2, 1976-01-01,	25.5,		57.7 \
1977-01-01,	35.8, 1977-01-01,	27.9,		63.7 \
1978-01-01,	40.4, 1978-01-01,	31.0,		71.4 \
1979-01-01,	48.1, 1979-01-01,	35.0,		83.1 \
1980-01-01,	54.4, 1980-01-01,	39.6,		94.0 \
1981-01-01,	64.8, 1981-01-01,	45.0,		109.8 \
1982-01-01,	72.7, 1982-01-01,	49.7,		122.4
datarnd = datarnd \
1983-01-01,	81.3, 1983-01-01,	55.2,		136.5 \
1984-01-01,	95.1, 1984-01-01,	62.2,		157.3 \
1985-01-01,	105., 1985-01-01,	71.0,		176.3 \
1986-01-01,	113., 1986-01-01,	75.2,		188.7 \
1987-01-01,	120., 1987-01-01,	81.7,		201.8 \
1988-01-01,	132., 1988-01-01,	85.1,		217.8 \
1989-01-01,	150., 1989-01-01,	87.8,		237.9 \
1990-01-01,	164., 1990-01-01,	91.1,		255.5 \
1991-01-01,	179., 1991-01-01,	91.4,		270.5 \
1992-01-01,	187., 1992-01-01,	91.6,		279.3 \
1993-01-01,	196., 1993-01-01,	91.4,		288.3 \
1994-01-01,	205., 1994-01-01,	92.0,		297.7 \
1995-01-01,	226., 1995-01-01,	94.0,		320.8 \
1996-01-01,	253., 1996-01-01,	95.5,		348.8 \
1997-01-01,	288., 1997-01-01,	98.3,		386.3 \
1998-01-01,	317., 1998-01-01,	102.3,		420.0 \
1999-01-01,	364., 1999-01-01,	106.5,		470.5 \
2000-01-01,	409., 2000-01-01,	113.2,		522.7 \
2001-01-01,	412., 2001-01-01,	119.7,		532.3 \
2002-01-01,	406., 2002-01-01,	126.0,		532.4 \
2003-01-01,	420., 2003-01-01,	136.5,		557.4 \
2004-01-01,	442., 2004-01-01,	145.5,		587.6 \
2005-01-01,	475., 2005-01-01,	153.9,		629.0 \
2006-01-01,	504., 2006-01-01,	160.6,		665.2 \
2007-01-01,	537., 2007-01-01,	169.0,		706.9 \
2008-01-01,	563., 2008-01-01,	177.2,		740.6 \
2009-01-01,	550., 2009-01-01,	179.8,		730.7
datarnd = datarnd \
2010-01-01,	564., 2010-01-01,	187.4,		751.7 \
2011-01-01,	592., 2011-01-01,	191.6,		783.8 \
2012-01-01,	621., 2012-01-01,	189.2,		810.9 \
2013-01-01,	649., 2013-01-01,	188.1,		838.0 \
2014-01-01,	690., 2014-01-01,	187.2,		877.2 \
2015-01-01,	728., 2015-01-01,	190.8,		919.4


// Title:               Gross Domestic Product: Research and Development
// Series ID:           Y694RC1A027NBEA
// Source:              US. Bureau of Economic Analysis
// Release:             Gross Domestic Product
// Seasonal Adjustment: Not Seasonally Adjusted
// Frequency:           Annual
// Units:               Billions of Dollars
// Date Range:          1929-01-01 to 2015-01-01
// Last Updated:        2016-03-25 1:31 PM CDT
// Notes:               BEA Account Code: Y694RC1
// https://research.stlouisfed.org/fred2/series/Y694RC1A027NBEA

databasic=
1929,    0.3 \
1930,    0.3 \
1931,    0.3 \
1932,    0.3 \
1933,    0.3 \
1934,    0.3 \
1935,    0.3 \
1936,    0.4 \
1937,    0.4 \
1938,    0.4 \
1939,    0.4 \
1940,    0.5 \
1941,    0.9 \
1942,    1.1 \
1943,    1.4 \
1944,    2.1 \
1945,    2.1 \
1946,    2.2 \
1947,    2.4 \
1948,    2.7 \
1949,    2.7 \
1950,    3.1 \
1951,    3.5 \
1952,    4.3 \
1953,    5.1 \
1954,    5.6 \
1955,    6.3 \
1956,    8.4 \
1957,    9.6 \
1958,   10.4 \
1959,   11.4 \
1960,   12.7 \
1961,   13.9 \
1962,   15.4 \
1963,   17.6 \
1964,   19.3 \
1965,   20.8 \
1966,   23.2 \
1967,   24.9 \
1968,   26.6 \
1969,   28.2 \
1970,   28.3 \
1971,   29.1 \
1972,   31.3 \
1973,   33.8 \
1974,   36.4 \
1975,     39 \
1976,   43.1 \
1977,   47.5 \
1978,   53.5 \
1979,   60.9 \
1980,   69.8 \
1981,   79.9 \
1982,     89 \
1983,   98.4 \
1984,  111.4 \
1985,  125.1 \
1986,  132.0 \
1987,  140.2 \
1988,  148.2 \
1989,  155.9 \
1990,  164.2 \
1991,  170.0 \
1992,  171.9 \
1993,  172.0 \
1994,  175.0 \
1995,  187.5 \
1996,  199.8 \
1997,  212.0 \
1998,  222.4 \
1999,  237.2 \
2000,  255.6 \
2001,  263.2 \
2002,  265.7 \
2003,  277.6 \
2004,  291.2 \
2005,  313.1 \
2006,  334.6 \
2007,  360.0 \
2008,  380.5 \
2009,  374.8 \
2010,  392.1 \
2011,  404.1 \
2012,  413.5 \
2013,  427.9 \
2014,  444.8 \
2015,  468.3


annualrnd=datarnd[,5]
annualyrs=datarnd[,1]:+2
annualbasic=databasic[,2]
yr0=1928


// BLS TFP Indexes
// From mfp_tables_historical-2017-02-17.xls
// Downloaded from BLS on 2/17/17
//   https://www.bls.gov/mfp/tables.htm
// Note: We will add back in the contributions of R&D and IPP

//Year	Multifactor Productivity	Contribution of Research and Development  Intensity	Contribution of All Other  Intellectual Property Products Intensity
// Year   TFP    RND     IPP
blsdata=			
1948,	44.186 ,	.,	. \
1949,	44.323 ,	.,	. \
1950,	47.495 ,	.,	. \
1951,	48.629 ,	.,	. \
1952,	49.466 ,	.,	. \
1953,	50.755 ,	.,	. \
1954,	50.728 ,	.,	. \
1955,	52.889 ,	.,	. \
1956,	52.523 ,	.,	. \
1957,	53.310 ,	.,	. \
1958,	53.528 ,	.,	. \
1959,	55.923 ,	.,	. \
1960,	56.267 ,	.,	. \
1961,	57.470 ,	.,	. \
1962,	59.560 ,	.,	. \
1963,	61.302 ,	.,	. \
1964,	63.698 ,	.,	. \
1965,	65.763 ,	.,	. \
1966,	67.773 ,	.,	. \
1967,	67.853 ,	.,	. \
1968,	69.562 ,	.,	. \
1969,	69.288 ,	.,	. \
1970,	69.134 ,	.,	. \
1971,	71.302 ,	.,	. \
1972,	73.337 ,	.,	. \
1973,	75.357 ,	.,	. \
1974,	72.747 ,	.,	. \
1975,	73.430 ,	.,	. \
1976,	76.096 ,	.,	. \
1977,	77.356 ,	.,	. \
1978,	78.401 ,	.,	. \
1979,	78.094 ,	.,	. \
1980,	76.292 ,	.,	. \
1981,	76.455 ,	.,	. \
1982,	73.887 ,	.,	. \
1983,	76.026 ,	.,	. \
1984,	78.305 ,	.,	. \
1985,	79.309 ,	.,	. \
1986,	80.563 ,	.,	. \
1987,	80.844 ,	98.040	, 95.473 \
1988,	81.622 ,	98.109	, 95.587 \
1989,	81.978 ,	98.175	, 95.763 \
1990,	82.686 ,	98.326	, 96.153 \
1991,	82.304 ,	98.519	, 96.599 \
1992,	85.038 ,	98.644	, 96.988 \
1993,	84.721 ,	98.674	, 97.295 \
1994,	85.246 ,	98.654	, 97.567 \
1995,	84.797 ,	98.671	, 97.909 \
1996,	86.421 ,	98.740	, 98.347 \
1997,	87.418 ,	98.778	, 98.818 \
1998,	88.517 ,	98.853	, 99.291 \
1999,	90.378 ,	98.922	, 99.617 \
2000,	91.957 ,	99.015	, 99.820 \
2001,	92.435 ,	99.165	, 99.967 \
2002,	94.397 ,	99.303	, 99.942 \
2003,	96.638 ,	99.386	, 99.675 \
2004,	99.089 ,	99.407	, 99.300 \
2005,	100.593,	99.420	, 99.064 \
2006,	100.917,	99.436	, 98.969 \
2007,	101.381,	99.516	, 99.052 \
2008,	100.174,	99.678	, 99.454 \
2009,	100.000,	100.000	, 100.00
blsdata = blsdata \ 
2010,	102.869,	100.094	, 99.952 \
2011,	102.986,	100.110	, 99.846 \
2012,	103.594,	100.116	, 99.798 \
2013,	103.841,	100.152	, 99.839 \
2014,	104.558,	100.176	, 99.927 \
2015,	104.798,	.,	.

blsyrs=blsdata[,1]
blslevel=blsdata[,2]
bls_rnd=blsdata[,3]
bls_ipp=blsdata[,4]

// Add back the RND and IPP contributions, which we do not want taken out
// (only starts taking out in 1987)
rawgrowth=delchange(log(blslevel))
rnd_contribution=delchange(log(bls_rnd))
ipp_contribution=delchange(log(bls_ipp))
_editmissing(rnd_contribution, 0)
_editmissing(ipp_contribution, 0)
tfpgrowth_rndipp=rawgrowth+rnd_contribution+ipp_contribution
tfpindex=1 \ exp(quadrunningsum(tfpgrowth_rndipp))
tfpindex=tfpindex/tfpindex[1]*blslevel[1]

""
"BLS Data for adding back RND and IPP"
blah=trimr((blslevel, tfpindex),0,1)
"Year BLS RND IPP Fixed BLSIndex TFPIndex"
trimr(blsyrs,0,1), rawgrowth, rnd_contribution, ipp_contribution, tfpgrowth_rndipp, blah
""

blsAvgYears= 1950 \ 
             1960 \ 
             1970 \ 
             1980 \ 
             1990 \ 
             2000 \ 
             2014
blsgrowth=delchange(log(blslevel[blsAvgYears:-1947])) :/ delchange(blsAvgYears)*100
tfpgrowth_fixed=delchange(log(tfpindex[blsAvgYears:-1947])) :/ delchange(blsAvgYears)*100
""
"Average Annual TFP Growth: Private Business Sector (BLS, by ending year)"
"Year BLSRaw BLS+RND"
trimr(blsAvgYears,1,0), blsgrowth, tfpgrowth_fixed


""
"Merging Gordon data for 1930s and 1940s with BLS data for 1950s and beyond"
years=GordonYears[1::2] \ trimr(blsAvgYears,1,0)
tfpgrowth=GordonGrowth[1::2] \ tfpgrowth_fixed
"Merged data:"
years, tfpgrowth

WageScientistData =  // with Roodman insertion of CPI-U in col. 2 https://www.usinflationcalculator.com/inflation/consumer-price-index-and-annual-percent-changes-from-1913-to-2008/
1929, 17.1,	2599     \
1930, 16.7,	2266     \
1931, 15.2,	1887     \
1932, 13.7,	1442     \
1933, 13,	1378     \
1934, 13.4,	1599     \
1935, 13.7,	1766     \
1936, 13.9,	2005     \
1937, 14.4,	2184     \
1938, 14.1,	2036     \
1939, 13.9,	2163     \
1940, 14,	2330.184 \
1941, 14.7,	2510.29  \
1942, 16.3,	2704.316 \
1943, 17.3,	2913.34  \
1944, 17.6,	3138.519 \
1945, 18,	3381.103 \
1946, 19.5,	3642.437 \
1947, 22.3,	3923.97  \
1948, 24.1,	4227.264 \
1949, 23.8,	4554     \
1950, 24.1,	4924.798 \
1951, 26,	5325.787 \
1952, 26.5,	5759.426 \
1953, 26.7,	6228.373 \
1954, 26.9,	6735.503 \
1955, 26.8,	7283.925 \
1956, 27.2,	7877     \
1957, 28.1,	8251.116 \
1958, 28.9,	8643     \
1959, 29.1,	9017.842 \
1960, 29.6,	9408.94  \
1961, 29.9,	9817     \
1962, 30.2,	9814     \
1963, 30.6,	9811     \
1964, 31,	10284    \
1965, 31.5,	10987.442\
1966, 32.4,	11739    \
1967, 33.4,	11924    \
1968, 34.8,	12938    \
1969, 36.7,	14079    \
1970, 38.8,	14434    \
1971, 40.5,	15133    \
1972, 41.8,	16201    \
1973, 44.4,	17064    \
1974, 49.3,	18265    \
1975, 53.8,	19111    \
1976, 56.9,	20516    \
1977, 60.6,	22125    \
1978, 65.2,	23724    \
1979, 72.6,	25544    \
1980, 82.4,	27216    \
1981, 90.9,	29278    \
1982, 96.5,	31055    \
1983, 99.6,	32472    \
1984, 103.9,	34736    \
1985, 107.6,	37570    \
1986, 109.6,	39773    \
1987, 113.6,	40840    \
1988, 118.3,	42861    \
1989, 124,	46932    \
1990, 130.7,	46961    \
1991, 136.2,	47350    \
1992, 140.3,	49116    \
1993, 144.5,	54682    \
1994, 148.2,	56298    \
1995, 152.4,	57018    \
1996, 156.9,	58527    \
1997, 160.5,	62718    \
1998, 163,	65444    \
1999, 166.6,	70174    \
2000, 172.2,	73156    \
2001, 177.1,	74027    \
2002, 179.9,	73229    \
2003, 184,	73244    \
2004, 188.9,	76376    \
2005, 195.3,	80130    \
2006, 201.6,	82827    \
2007, 207.3,	83127    \
2008, 215.303,	84694    \
2009, 214.537,	82301    \
2010, 218.056,	81653    \
2011, 224.939,	84112
WageScientistData = WageScientistData \
2012, 229.594,	85602    \
2013, 232.957,	87681    \
2014, 236.736,	90276    \
2015, 237.017,	92947.801


WageYears = WageScientistData[,1]
WageScientist = WageScientistData[,2]  // change 2 to 3 to use scientist wages instead of CPI
WageSci=WageScientist[annualyrs:-WageYears[1]:+1]

yrs= 1930, 1940 \
     1940, 1950 \
     1950, 1960 \
     1960, 1970 \
     1970, 1980 \
     1980, 1990 \
     1990, 2000 \
     2000, 2015

T=rows(yrs)
IdeaTFP=J(T,1,0)


//RND=exp(average(log(annualrnd),yrs-yr0)) // Geometric average
//Scientists=RND:/WageSci
RND=annualrnd
ScientistsAnnual=RND:/WageSci:*10^9
RNDBasic=annualbasic
ScientistsBasic=RNDBasic:/WageSci:*10^9

Scientists=exp(panelsum(log(ScientistsAnnual), yrs:-yr0) :/ (yrs[,2] :- yrs[,1] :+ 1)) // Geometric Average
Scientists=Scientists/Scientists[1]
Scientists=Scientists:^Lambda // Adjusting effective research by lambda
IdeaTFP=tfpgrowth:/Scientists
IdeaTFP=IdeaTFP/IdeaTFP[1]

""
"Underlying annual data"
"Year IPP($b) WageSci ScientistsIpp SciBasic"
annualyrs, RND, WageSci, ScientistsAnnual, ScientistsBasic

""
"RESULTS"
"TFPGrowth Scientists IdeaTFP"
tfpgrowth, Scientists, IdeaTFP*100
""
printf("Mean tfpgrowth = %f\n",mean(tfpgrowth))

""
TT=.5*(2000+2015) - .5*(1940+1930)
gSci =100*log(Scientists[rows(Scientists)]/Scientists[1])/TT
giTFP=100*log(IdeaTFP[rows(IdeaTFP)]/IdeaTFP[1])/TT
printf("The implied average growth rate of scientists is %f percent\n",gSci)
//fprintf("                            scientists^lambda is %f\n",gSci*Lambda)
printf("The implied average growth rate of idea TFP   is %f\n",giTFP)
printf("      Half life of idea TFP = %f\n",-log(2)/giTFP*100)
printf("      Implied beta =%f\n",-giTFP/mean(tfpgrowth))

""
printf("Idea TFP Factor decline, 1930s - 2000s = %f\n",IdeaTFP[1]/IdeaTFP[rows(IdeaTFP)])
printf("Idea TFP Factor decline, 1960s - 2000s = %f\n",IdeaTFP[4]/IdeaTFP[rows(IdeaTFP)])
printf("Idea TFP Factor decline, 1980s - 2000s = %f\n",IdeaTFP[6]/IdeaTFP[rows(IdeaTFP)])
""
printf("Scientists factor increase, 1930s - 2000s = %f\n",Scientists[rows(Scientists)]/Scientists[1])
printf("Scientists factor increase, 1960s - 2000s = %f\n",Scientists[rows(Scientists)]/Scientists[4])
printf("Scientists factor increase, 1980s - 2000s = %f\n",Scientists[rows(Scientists)]/Scientists[6])
end
