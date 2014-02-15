/*****************************************************
* This code is distributed under a Non-commercial use 
* license.  For details see LICENSE.  Use of this
* code must be properly attributed to its author
* Gregory E. Sims provided that its use or derivative 
* use is non-commercial in nature.  Proper attribution        
* can be made by citing:
*
* Sims GE, et al (2009) Alignment-free genome 
* comparison with feature frequency profiles (FFP) and 
* optimal resolutions. Proc. Natl. Acad. Sci. USA.
* 106, 2677-82.
*
* Gregory E. Sims (C) 2010-2012
*
*****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "hashroll.h"
#include "utils.h"
#include "../config.h"

#define A 16807
#define IAMOD2P16 22039
#define MOD2P16 0x0000ffff
extern char *weightVector;


/**
 *
 * Returns an unmasked rolling hash value for an RY coded feature
 * using a pre-existing hash value.
 * 
 * This is the hashing function used to calculate
 * a rolling hash.  Given a pre-calculated hash
 * of length l, we can quickly compute what the
 * new hash value would be by shifting the characters
 * to the left and adding a new character to the right
 *
 *    Example:
 *    ATATAT   has a precaclucated hash value
 *    TATATG   What is the hash of this?
 *
 *
 * We use the Rabin-Karp rolling hash
 * Karb RM Rabin MO (1987) Efficient randomized pattern-matching
 * algorithms.  IBM Journal of R&D. 31:2,249-260.
 *
 * H=s1a^k-1 + s2a^k-2 + ...   ska^0
 *
 * In this case a=3
 *
 * To find hash of TATATG multiply H
 * by a and subtract out s1a^k and add
 * sk
 *
 * @param s An RY-coded character string.
 * @param hash Pre-existing hash value
 * @param length of string
 */

// make this a macro

void resetHash(HASH * h)
{
    h->st_hash = h->rt_hash = 0;
    h->numChar = 0;
}


void init(HASH * h, int isMasked, int mode, bool isClass, bool reverse, int k)
{
    h->k = k;
    h->s = (char *)chkmalloc(sizeof(char),k+1);
    h->r = (char *)chkmalloc(sizeof(char),k+1);
    h->s[h->k] = '\0';
    h->r[h->k] = '\0';
    h->st_hash = h->rt_hash = 0;
    h->numChar = 0;
    h->isClass = isClass;
    h->mode=mode;
    h->reverse = reverse;
    memset(h->table, 0, sizeof(NODE *) * BUCKETS);
    h->keyN = 0;

    if (mode == nucleotide) {
	if (isClass)
	    h->hashi = base_ry_hash_values;
	else
	    h->hashi = base_hash_values;
    } else if (mode == amino) {
	if (isClass)
	    h->hashi = aac_values;
	else
	    h->hashi = aa_hash_values;
    } else if (mode == text) {
	h->hashi = txt_hash_values;
    }


}



/* This is a static table of hash increments 
 * Eeach term to be subtracted or
 * added to the hash has been pre-calculated 
 * in mod form. 
 *
 * ci*a^k  
 * 
 * ci= the new position index value = {1..4}
 *  a= 16807
 *  k= the position of the new character in
 *     the k-mer.
 *
 * All possibilities have been precalculated
 * and then mod 2^16, using the bc language
 * arbitrary precision calculator.
 * 
 * Element apnmod[1][3] for example contains
 * the results of 1*16807^3%2^16 which gives
 * 42039.  This should save plenty of 
 * calculation time.
 *
 */

static const unsigned apnmod[5][41] = { {0},
{1, 16807, 15089, 42039, 5857, 3527, 33745, 3671, 29121, 13799, 53425, 5239,
 37025, 14855, 41361, 13975, 62337, 39463, 30321, 62647, 6753, 54855, 53073,
 52951, 34113, 28263, 11313, 17655, 46113, 57991, 3345, 55063, 9985, 45735,
 61937, 1335, 24033, 24263, 23249, 20311, 55489},
{2, 33614, 30178, 18542, 11714, 7054, 1954, 7342, 58242, 27598, 41314, 10478,
 8514, 29710, 17186, 27950, 59138, 13390, 60642, 59758, 13506, 44174, 40610,
 40366, 2690, 56526, 22626, 35310, 26690, 50446, 6690, 44590, 19970, 25934,
 58338, 2670, 48066, 48526, 46498, 40622, 45442},
{3, 50421, 45267, 60581, 17571, 10581, 35699, 11013, 21827, 41397, 29203, 15717,
 45539, 44565, 58547, 41925, 55939, 52853, 25427, 56869, 20259, 33493, 28147,
 27781, 36803, 19253, 33939, 52965, 7267, 42901, 10035, 34117, 29955, 6133,
 54739,
 4005, 6563, 7253, 4211, 60933, 35395},
{4, 1692, 60356, 37084, 23428, 14108, 3908, 14684, 50948, 55196, 17092, 20956,
 17028, 59420, 34372, 55900, 52740, 26780, 55748, 53980, 27012, 22812, 15684,
 15196, 5380, 47516, 45252, 5084, 53380, 35356, 13380, 23644, 39940, 51868,
 51140, 5340, 30596, 31516, 27460, 15708, 25348}
};


static const unsigned apnmodaa[21][41] = { {0},
{1, 16807, 15089, 42039, 5857, 3527, 33745, 3671, 29121, 13799, 53425, 5239,
 37025, 14855, 41361, 13975, 62337, 39463, 30321, 62647, 6753, 54855, 53073,
 52951, 34113, 28263, 11313, 17655, 46113, 57991, 3345, 55063, 9985, 45735,
 61937, 1335, 24033, 24263, 23249, 20311, 55489},
{2, 33614, 30178, 18542, 11714, 7054, 1954, 7342, 58242, 27598, 41314, 10478,
 8514, 29710, 17186, 27950, 59138, 13390, 60642, 59758, 13506, 44174, 40610,
 40366, 2690, 56526, 22626, 35310, 26690, 50446, 6690, 44590, 19970, 25934,
 58338, 2670, 48066, 48526, 46498, 40622, 45442},
{3, 50421, 45267, 60581, 17571, 10581, 35699, 11013, 21827, 41397, 29203, 15717,
 45539, 44565, 58547, 41925, 55939, 52853, 25427, 56869, 20259, 33493, 28147,
 27781, 36803, 19253, 33939, 52965, 7267, 42901, 10035, 34117, 29955, 6133,
 54739,
 4005, 6563, 7253, 4211, 60933, 35395},
{4, 1692, 60356, 37084, 23428, 14108, 3908, 14684, 50948, 55196, 17092, 20956,
 17028, 59420, 34372, 55900, 52740, 26780, 55748, 53980, 27012, 22812, 15684,
 15196, 5380, 47516, 45252, 5084, 53380, 35356, 13380, 23644, 39940, 51868,
 51140, 5340, 30596, 31516, 27460, 15708, 25348},
{5, 18499, 9909, 13587, 29285, 17635, 37653, 18355, 14533,
 3459, 4981, 26195, 54053, 8739, 10197, 4339, 49541, 707,
 20533, 51091, 33765, 12131, 3221, 2611, 39493, 10243,
 56565, 22739, 33957, 27811, 16725, 13171, 49925, 32067,
 47541, 6675, 54629, 55779, 50709, 36019},
{6, 35306, 24998, 55626, 35142, 21162, 5862, 22026, 43654,
 17258, 58406, 31434, 25542, 23594, 51558, 18314, 46342,
 40170, 50854, 48202, 40518, 1450, 56294, 55562, 8070,
 38506, 2342, 40394, 14534, 20266, 20070, 2698, 59910,
 12266, 43942, 8010, 13126, 14506, 8422, 56330},
{7, 52113, 40087, 32129, 40999, 24689, 39607, 25697, 7239,
 31057, 46295, 36673, 62567, 38449, 27383, 32289, 43143,
 14097, 15639, 45313, 47271, 56305, 43831, 42977, 42183,
 1233, 13655, 58049, 60647, 12721, 23415, 57761, 4359,
 58001, 40343, 9345, 37159, 38769, 31671, 11105},
{8, 3384, 55176, 8632, 46856, 28216, 7816, 29368, 36360,
 44856, 34184, 41912, 34056, 53304, 3208, 46264, 39944,
 53560, 45960, 42424, 54024, 45624, 31368, 30392, 10760,
 29496, 24968, 10168, 41224, 5176, 26760, 47288, 14344,
 38200, 36744, 10680, 61192, 63032, 54920, 31416},
{9, 20191, 4729, 50671, 52713, 31743, 41561, 33039, 65481,
 58655, 22073, 47151, 5545, 2623, 44569, 60239, 36745,
 27487, 10745, 39535, 60777, 34943, 18905, 17807, 44873,
 57759, 36281, 27823, 21801, 63167, 30105, 36815, 24329,
 18399, 33145, 12015, 19689, 21759, 12633, 51727},
{10, 36998, 19818, 27174, 58570, 35270, 9770, 36710,
 29066, 6918, 9962, 52390, 42570, 17478, 20394, 8678,
 33546, 1414, 41066, 36646, 1994, 24262, 6442, 5222,
 13450, 20486, 47594, 45478, 2378, 55622, 33450, 26342,
 34314, 64134, 29546, 13350, 43722, 46022, 35882, 6502},
{11, 53805, 34907, 3677, 64427, 38797, 43515, 40381, 58187,
 20717, 63387, 57629, 14059, 32333, 61755, 22653, 30347,
 40877, 5851, 33757, 8747, 13581, 59515, 58173, 47563,
 48749, 58907, 63133, 48491, 48077, 36795, 15869, 44299,
 44333, 25947, 14685, 2219, 4749, 59131, 26813},
{12, 5076, 49996, 45716, 4748, 42324, 11724, 44052, 21772,
 34516, 51276, 62868, 51084, 47188, 37580, 36628, 27148,
 14804, 36172, 30868, 15500, 2900, 47052, 45588, 16140,
 11476, 4684, 15252, 29068, 40532, 40140, 5396, 54284,
 24532, 22348, 16020, 26252, 29012, 16844, 47124},
{13, 21883, 65085, 22219, 10605, 45851, 45469, 47723,
 50893, 48315, 39165, 2571, 22573, 62043, 13405, 50603,
 23949, 54267, 957, 27979, 22253, 57755, 34589, 33003,
 50253, 39739, 15997, 32907, 9645, 32987, 43485, 60459,
 64269, 4731, 18749, 17355, 50285, 53275, 40093, 1899},
{14, 38690, 14638, 64258, 16462, 49378, 13678, 51394,
 14478, 62114, 27054, 7810, 59598, 11362, 54766, 64578,
 20750, 28194, 31278, 25090, 29006, 47074, 22126, 20418,
 18830, 2466, 27310, 50562, 55758, 25442, 46830, 49986,
 8718, 50466, 15150, 18690, 8782, 12002, 63342, 22210},
{15, 55497, 29727, 40761, 22319, 52905, 47423, 55065,
 43599, 10377, 14943, 13049, 31087, 26217, 30591, 13017,
 17551, 2121, 61599, 22201, 35759, 36393, 9663, 7833,
 52943, 30729, 38623, 2681, 36335, 17897, 50175, 39513,
 18703, 30665, 11551, 20025, 32815, 36265, 21055, 42521},
{16, 6768, 44816, 17264, 28176, 56432, 15632, 58736, 7184,
 24176, 2832, 18288, 2576, 41072, 6416, 26992, 14352, 41584,
 26384, 19312, 42512, 25712, 62736, 60784, 21520, 58992,
 49936, 20336, 16912, 10352, 53520, 29040, 28688, 10864,
 7952, 21360, 56848, 60528, 44304, 62832},
{17, 23575, 59905, 59303, 34033, 59959, 49377, 62407, 36305,
 37975, 56257, 23527, 39601, 55927, 47777, 40967, 11153,
 15511, 56705, 16423, 49265, 15031, 50273, 48199, 55633,
 21719, 61249, 37991, 63025, 2807, 56865, 18567, 38673, 56599,
 4353, 22695, 15345, 19255, 2017, 17607},
{18, 40382, 9458, 35806, 39890, 63486, 17586, 542, 65426, 51774,
 44146, 28766, 11090, 5246, 23602, 54942, 7954, 54974, 21490,
 13534, 56018, 4350, 37810, 35614, 24210, 49982, 7026, 55646,
 43602, 60798, 60210, 8094, 48658, 36798, 754, 24030, 39378,
 43518, 25266, 37918},
{19, 57189, 24547, 12309, 45747, 1477,
 51331, 4213, 29011, 37, 32035, 34005, 48115, 20101, 64963, 3381,
 4755, 28901, 51811, 10645, 62771, 59205, 25347, 23029, 58323,
 12709, 18339, 7765, 24179, 53253, 63555, 63157, 58643, 16997,
 62691, 25365, 63411, 2245, 48515, 58229},
{20, 8460, 39636, 54348, 51604, 5004, 19540, 7884, 58132, 13836,
 19924, 39244, 19604, 34956, 40788, 17356, 1556, 2828, 16596,
 7756, 3988, 48524, 12884, 10444, 26900, 40972, 29652, 25420,
 4756, 45708, 1364, 52684, 3092, 62732, 59092, 26700, 21908,
 26508, 6228, 13004}
};


// can replace aa and nuc table with this master table.

static const unsigned apnmodtxt[27][41] = { {0},
{1,16807,15089,42039,5857,3527,33745,3671,
 29121,13799,53425,5239,37025,14855,41361,
 13975,62337,39463,30321,62647,6753,54855,
 53073,52951,34113,28263,11313,17655,46113,
 57991,3345,55063,9985,45735,61937,1335,
 24033,24263,23249,20311,55489},
{2,33614,30178,18542,11714,7054,1954,7342,
 58242,27598,41314,10478,8514,29710,17186,
 27950,59138,13390,60642,59758,13506,44174,
 40610,40366,2690,56526,22626,35310,26690,
 50446,6690,44590,19970,25934,58338,2670,
 48066,48526,46498,40622,45442},
{3,50421,45267,60581,17571,10581,35699,11013,
 21827,41397,29203,15717,45539,44565,58547,
 41925,55939,52853,25427,56869,20259,33493,
 28147,27781,36803,19253,33939,52965,7267,
 42901,10035,34117,29955,6133,54739,4005,
 6563,7253,4211,60933,35395},
{4,1692,60356,37084,23428,14108,3908,14684,
 50948,55196,17092,20956,17028,59420,34372,
55900,52740,26780,55748,53980,27012,22812,
15684,15196,5380,47516,45252,5084,53380,
35356,13380,23644,39940,51868,51140,5340,
30596,31516,27460,15708,25348},
{5,18499,9909,13587,29285,17635,37653,18355,
14533,3459,4981,26195,54053,8739,10197,4339,
49541,707,20533,51091,33765,12131,3221,2611,
39493,10243,56565,22739,33957,27811,16725,
13171,49925,32067,47541,6675,54629,55779,
50709,36019,15301},
{6,35306,24998,55626,35142,21162,5862,22026,
43654,17258,58406,31434,25542,23594,51558,
18314,46342,40170,50854,48202,40518,1450,
56294,55562,8070,38506,2342,40394,14534,
20266,20070,2698,59910,12266,43942,8010,
13126,14506,8422,56330,5254},
{7,52113,40087,32129,40999,24689,39607,25697,
7239,31057,46295,36673,62567,38449,27383,32289,
43143,14097,15639,45313,47271,56305,43831,42977,
42183,1233,13655,58049,60647,12721,23415,57761,
4359,58001,40343,9345,37159,38769,31671,11105,
60743},
{8,3384,55176,8632,46856,28216,7816,29368,36360,
44856,34184,41912,34056,53304,3208,46264,39944,
53560,45960,42424,54024,45624,31368,30392,10760,
29496,24968,10168,41224,5176,26760,47288,14344,
38200,36744,10680,61192,63032,54920,31416,50696},
{9,20191,4729,50671,52713,31743,41561,33039,65481,
58655,22073,47151,5545,2623,44569,60239,36745,
27487,10745,39535,60777,34943,18905,17807,44873,
57759,36281,27823,21801,63167,30105,36815,24329,
18399,33145,12015,19689,21759,12633,51727,40649},
{10,36998,19818,27174,58570,35270,9770,36710,
29066,6918,9962,52390,42570,17478,20394,8678,
33546,1414,41066,36646,1994,24262,6442,5222,
13450,20486,47594,45478,2378,55622,33450,26342,
34314,64134,29546,13350,43722,46022,35882,6502,
30602},
{11,53805,34907,3677,64427,38797,43515,40381,
58187,20717,63387,57629,14059,32333,61755,22653,
30347,40877,5851,33757,8747,13581,59515,58173,
47563,48749,58907,63133,48491,48077,36795,15869,
44299,44333,25947,14685,2219,4749,59131,26813,
20555},
{12,5076,49996,45716,4748,42324,11724,44052,21772,
34516,51276,62868,51084,47188,37580,36628,27148,
14804,36172,30868,15500,2900,47052,45588,16140,
11476,4684,15252,29068,40532,40140,5396,54284,
24532,22348,16020,26252,29012,16844,47124,10508},
{13,21883,65085,22219,10605,45851,45469,47723,
50893,48315,39165,2571,22573,62043,13405,50603,
23949,54267,957,27979,22253,57755,34589,33003,
50253,39739,15997,32907,9645,32987,43485,60459,
64269,4731,18749,17355,50285,53275,40093,1899,461},
{14,38690,14638,64258,16462,49378,13678,51394,
14478,62114,27054,7810,59598,11362,54766,64578,
20750,28194,31278,25090,29006,47074,22126,20418,
18830,2466,27310,50562,55758,25442,46830,49986,
8718,50466,15150,18690,8782,12002,63342,22210,55950},
{15,55497,29727,40761,22319,52905,47423,55065,
43599,10377,14943,13049,31087,26217,30591,13017,
17551,2121,61599,22201,35759,36393,9663,7833,
52943,30729,38623,2681,36335,17897,50175,39513,
18703,30665,11551,20025,32815,36265,21055,42521,
45903},
{16,6768,44816,17264,28176,56432,15632,58736,
7184,24176,2832,18288,2576,41072,6416,26992,
14352,41584,26384,19312,42512,25712,62736,60784,
21520,58992,49936,20336,16912,10352,53520,29040,
28688,10864,7952,21360,56848,60528,44304,62832,
35856},
{17,23575,59905,59303,34033,59959,49377,62407,
36305,37975,56257,23527,39601,55927,47777,40967,
11153,15511,56705,16423,49265,15031,50273,48199,
55633,21719,61249,37991,63025,2807,56865,18567,
38673,56599,4353,22695,15345,19255,2017,17607,
25809},
{18,40382,9458,35806,39890,63486,17586,542,
65426,51774,44146,28766,11090,5246,23602,54942,
7954,54974,21490,13534,56018,4350,37810,35614,
24210,49982,7026,55646,43602,60798,60210,8094,
48658,36798,754,24030,39378,43518,25266,37918,
15762},
{19,57189,24547,12309,45747,1477,51331,4213,
29011,37,32035,34005,48115,20101,64963,3381,
4755,28901,51811,10645,62771,59205,25347,
23029,58323,12709,18339,7765,24179,53253,
63555,63157,58643,16997,62691,25365,63411,
2245,48515,58229,5715},
{20,8460,39636,54348,51604,5004,19540,7884,
58132,13836,19924,39244,19604,34956,40788,
17356,1556,2828,16596,7756,3988,48524,12884,
10444,26900,40972,29652,25420,4756,45708,
1364,52684,3092,62732,59092,26700,21908,
26508,6228,13004,61204},
{21,25267,54725,30851,57461,8531,53285,11555,
21717,27635,7813,44483,56629,49811,16613,31331,
63893,42291,46917,4867,10741,37843,421,63395,
61013,3699,40965,43075,50869,38163,4709,42211,
13077,42931,55493,28035,45941,50771,29477,
33315,51157},
{22,42074,4278,7354,63318,12058,21494,15226,
50838,41434,61238,49722,28118,64666,57974,
45306,60694,16218,11702,1978,17494,27162,
53494,50810,29590,31962,52278,60730,31446,
30618,8054,31738,23062,23130,51894,29370,
4438,9498,52726,53626,41110},
{23,58881,19367,49393,3639,15585,55239,18897,
14423,55233,49127,54961,65143,13985,33799,
59281,57495,55681,42023,64625,24247,16481,
41031,38225,63703,60225,63591,12849,12023,
23073,11399,21265,33047,3329,48295,30705,
28471,33761,10439,8401,31063},
{24,10152,34456,25896,9496,19112,23448,22568,
43544,3496,37016,60200,36632,28840,9624,7720,
54296,29608,6808,61736,31000,5800,28568,25640,
32280,22952,9368,30504,58136,15528,14744,10792,
43032,49064,44696,32040,52504,58024,33688,28712,
21016},
{25,26959,49545,2399,15353,22639,57193,26239,7129,
17295,24905,65439,8121,43695,50985,21695,51097,3535,
37129,58847,37753,60655,16105,13055,857,51215,20681,
48159,38713,7983,18089,319,53017,29263,41097,33375,
11001,16751,56937,49023,10969},
{26,43766,64634,44438,21210,26166,25402,29910,36250,
31094,12794,5142,45146,58550,26810,35670,47898,42998,
1914,55958,44506,49974,3642,470,34970,13942,31994,278,
19290,438,21434,55382,63002,9462,37498,34710,35034,
41014,14650,3798,922}
};

//use config.h and configure to take care of these define statements
//On Cygwin GNU this is part of the library
char *strupr(char *);

#ifndef HAVE_STRUPR
//On Linux 2.6.18-194.11.3.el5 this is not part of the library
char *strupr(char *s)
{
    char *t = (char *)chkmalloc(sizeof(char),strlen(s));
    char *r;
    strcpy(t, s);
    r = t;
    while (*t != '\0') {
	*t = toupper((int) *t);
	t++;
    }
    return r;
}
#endif



// buckets is 20013 
// inverse of 5 mod buckets is 12008 mod buckets

// It seems possible to implement the building of
// the hash value in the reverse direction so that
// we don't need to know the length of the string.

int hashAdd(HASH * h, char *s, unsigned val)
{
    int i = 0;
    int k = h->k;
    NODE *ptr;
    NODE *r;
    unsigned int idx = 0;

    h->st_hash = h->rt_hash = 0;
    // isolate in its own hash function.

    if (h->isClass)
	ry(s, h->k);

    while (k > 0) {
	k--;
	idx += apnmod[h->hashi[(unsigned char) s[i]]][k];
	idx &= MOD2P16;
	i++;
    }


    if (h->reverse) {
	strcpy(h->r, s);
	rev(h->r, h->k);
	complement(h->r, h->k);
	k = h->k;
	i = 0;
	while (k > 0) {
	    k--;
	    h->rt_hash += apnmod[h->hashi[(unsigned char) h->r[i]]][k];
	    h->rt_hash &= MOD2P16;
	    i++;
	}
	if (h->rt_hash < idx) {
	    idx = h->rt_hash;
	    s = h->r;
	}
    }

    ptr = h->table[idx];
    while (ptr != NULL) {
	if (strcmp(ptr->key, s) == 0) {
	    ptr->value += val;
	    return 1;
	}
	ptr = ptr->next;
    }

    r = (NODE *)chkmalloc(sizeof(NODE),1);
    r->key = (char *)chkmalloc(sizeof(char), h->k + 1);
    strcpy(r->key, strupr(s));
    r->value = val;
    r->next = h->table[idx];
    h->table[idx] = r;
    h->keyN++;
    return -1;
}

//@todo test that string reversal w/ mask is working properly

int hashAddw(HASH * h, char *s, unsigned val)
{
    int i = 0;
    int k = h->k;
    NODE *ptr;
    NODE *r;
    unsigned int idx = 0;
    int j;

    // reseting the hash might be unnecessary.
    // aslo using the hash types internal values 
    // might not be needed, that way pushes and lookups can be mixed.
    if (h->isClass)
	ry(s, h->k);

    while (k > 0) {
	k--;
	idx += apnmod[h->hashi[(unsigned char) s[i]]][k];
	idx &= MOD2P16;
	i++;
    }
    if (h->reverse) {
	h->rt_hash = 0;
	strcpy(h->r, s);
	rev(h->r, h->k);
	complement(h->r, h->k);
	k = h->k;
	i = 0;
	while (k > 0) {
	    k--;
	    h->rt_hash += apnmod[h->hashi[(unsigned char) h->r[i]]][k];
	    h->rt_hash &= MOD2P16;
	    i++;
	}
	if (h->rt_hash < idx) {
	    idx = h->rt_hash;
	    s = h->r;
	}
    }


    ptr = h->table[idx];

    for (j = 0; j < h->k; j++)
	if (weightVector[j] == '0') {
	    idx -= apnmod[h->hashi[(unsigned char) s[j]]][h->k - j - 1];
	    idx &= MOD2P16;
	}


    while (ptr != NULL) {
	if (new_strcmpw(ptr->key, s) == 0) {
	    ptr->value += val;
	    return 1;
	}
	ptr = ptr->next;
    }

    r = (NODE *) chkmalloc(sizeof(NODE),1);
    r->key = (char *) chkmalloc(sizeof(char), h->k + 1);
    strcpy(r->key, strupr(s));
    r->value = val;
    r->next = h->table[idx];
    h->table[idx] = r;
    h->keyN++;
    return -1;

}


int hashAddaa(HASH * h, char *s, unsigned val)
{
    int i = 0;
    int k = h->k;
    NODE *ptr;
    NODE *r;
    unsigned int idx = 0;

    h->st_hash = 0;
    // isolate in its own hash function.

    if (h->isClass)
	_class(s, h->k);

    while (k > 0) {
	k--;
	idx += apnmodaa[h->hashi[(unsigned char) s[i]]][k];
	idx &= MOD2P16;
	i++;
    }


    ptr = h->table[idx];
    while (ptr != NULL) {
	if (strcmp(ptr->key, s) == 0) {
	    ptr->value += val;
	    return 1;
	}
	ptr = ptr->next;
    }

    r = (NODE *) chkmalloc(sizeof(NODE),1);
    r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
    strcpy(r->key, strupr(s));
    r->value = val;
    r->next = h->table[idx];
    h->table[idx] = r;
    h->keyN++;
    return -1;
}




int hashAddtxt(HASH * h, char *s, unsigned val)
{
    int i = 0;
    int k = h->k;
    NODE *ptr;
    NODE *r;
    unsigned int idx = 0;

    h->st_hash = 0;

    while (k > 0) {
	k--;
	idx += apnmodtxt[h->hashi[(unsigned char) s[i]]][k];
	idx &= MOD2P16;
	i++;
    }


    ptr = h->table[idx];
    while (ptr != NULL) {
	if (strcmp(ptr->key, s) == 0) {
	    ptr->value += val;
	    return 1;
	}
	ptr = ptr->next;
    }

    r = (NODE *) chkmalloc(sizeof(NODE),1);
    r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
    strcpy(r->key, strupr(s));
    r->value = val;
    r->next = h->table[idx];
    h->table[idx] = r;
    h->keyN++;
    return -1;
}





//Note these should be changed to Assign val.
int hashAddwaa(HASH * h, char *s, unsigned val)
{
    int i = 0;
    int k = h->k;
    NODE *ptr;
    NODE *r;
    unsigned int idx = 0;
    int j;

    if (h->isClass)
	_class(s, h->k);

    while (k > 0) {
	k--;
	idx += apnmodaa[h->hashi[(unsigned char) s[i]]][k];
	idx &= MOD2P16;
	i++;
    }


    for (j = 0; j < h->k; j++)
	if (weightVector[j] == '0') {
	    idx -= apnmodaa[h->hashi[(unsigned char) s[j]]][h->k - j - 1];
	    idx &= MOD2P16;
	}

    ptr = h->table[idx];
    while (ptr != NULL) {
	if (new_strcmpw(ptr->key, s) == 0) {
	    ptr->value += val;
	    return 1;
	}
	ptr = ptr->next;
    }

    r = (NODE *) chkmalloc(sizeof(NODE),1);
    r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
    strcpy(r->key, strupr(s));
    r->value = val;
    r->next = h->table[idx];
    h->table[idx] = r;
    h->keyN++;
    return 0;
}


unsigned int hashValNuc(HASH * h, register char *s)
{
    int i = 0;
    int k = h->k;
    NODE *ptr;
    h->st_hash = h->rt_hash = 0;
    unsigned int idx;
    // isolate in its own hash function.

    while (k > 0) {
	k--;
	h->st_hash += apnmod[h->hashi[(unsigned char) s[i]]][k];
	h->st_hash &= MOD2P16;
	i++;
    }
    idx = h->st_hash;

    //only applies to Nucleotides.
    if (h->reverse) {
	strcpy(h->r, s);
	rev(h->r, h->k);
	complement(h->r, h->k);
	k = h->k;
	i = 0;
	while (k > 0) {
	    k--;
	    h->rt_hash += apnmod[h->hashi[(unsigned char) h->r[i]]][k];
	    h->rt_hash &= MOD2P16;
	    i++;
	}
	if (h->rt_hash < idx) {
	    idx = h->rt_hash;
	    s = h->r;
	}
    }
    ptr = h->table[idx];

    while (ptr != NULL) {
	if (strcmp(ptr->key, s) == 0)
	    return ptr->value;
	ptr = ptr->next;
    }

    return 0;
}


int hashAssign(HASH * h, register char *s, unsigned int val)
{
    int i = 0;
    int k = h->k;
    NODE *ptr;
    h->st_hash = 0;
    static const unsigned mask = 0x0000ffff;

    // isolate in its own hash function.

    while (k > 0) {
	k--;
	h->st_hash += apnmod[h->hashi[(unsigned char) s[i]]][k];
	h->st_hash &= mask;
	i++;
    }

    ptr = h->table[h->st_hash];

    while (ptr != NULL) {
	if ((*(h->strcmpf)) (ptr->key, s) == 0) {
	    ptr->value = val;
	    return 1;
	}
	ptr = ptr->next;
    }

    return 0;
}





/* Established a new convention:  We no longer make lookups to check whether
 * the reverse complement feature is stored in  
 * the hash, we simply do this:  Scanning and hashing of the sequence is done
 * in the forward direction. 
 * We make the decision to physically store a key as in the reverse or forward direction
 * within the hash table: using this criterion:
 * Whichever word would come first if alphabetically sorted is placed in the hash.  
 * Therefore ATG is printed instead of its reverse complement word CAT. With a reverse complement
 * palindrome of course both are equivalent.  By doing this we can half the number
 * of lookups. */

// Use define templates to reduce the overall code size

void pushatgc(HASH * h, char *c, int n,bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    NODE *r;
    static int inHeader = 0;
    unsigned char index;
    char *s;
    unsigned *idx;
    extern bool mflag; // global value indicating process multiple headers
    int i;
    
    for (i = 0; i < n; i++) {
	// Buffer and process several bases at once.
	// Skip over header defline.
	// 
	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) 
			i++;

         // check to see if we're done processing header
	    if (i < n) 
		inHeader = 0;
	    else 
		return;
	
         //@todo add warning for not finding any keys.	    
	 // In this case no warnings will be produced when
	 // no keys are found.
	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;

	}

	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = atgc_to_ry(c[i]);
	//invalid character reset hash

	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;
	}


	c[i] = toupper((int) c[i]);

	if (h->numChar == h->k) {
	    h->st_hash -= apnmod[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	    if (h->reverse)
		h->rt_hash -=
		    apnmod[h->hashi[(unsigned char) h->r[h->k - 1]]][0];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];



	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;

	if (h->reverse) {	//save time by using function w/o test

	    memmove(&h->r[1], h->r, h->k - 1);
	    h->r[0] = flip(c[i]);

	    h->rt_hash *= IAMOD2P16;
	    h->rt_hash &= MOD2P16;
	    h->rt_hash += apnmod[h->hashi[(unsigned char) c[i]]][h->k - 1];
	    h->rt_hash &= MOD2P16;
	}


	if (h->numChar == h->k) {

	    s = h->s;
	    idx = &h->st_hash;

	    if (h->reverse) {
		if (h->st_hash > h->rt_hash) {
		    s = h->r;
		    idx = &h->rt_hash;
		}
	    }


	    ptr = h->table[*idx];

	    while (ptr != NULL) {
		    if (strcmp(ptr->key, s) == 0) {
			ptr->value++;	// to make generic.
			break;
		    }
		    ptr = ptr->next;
		}

	     //correct this error in all instances.

	     if (ptr == NULL) {  // Possibly redesign
	     r = (NODE *) chkmalloc(sizeof(NODE),1);
	     r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
	     strcpy(r->key, s);
	     r->value = 1;
	     r->next = h->table[*idx];
             h->table[*idx] = r;
	     h->keyN++;
	     }
	     
	}
    }
    //should never see
    return;
}


void chkpushatgc(HASH * h, char *c, int n,bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    static int inHeader = 0;
    unsigned char index;
    char *s;
    unsigned *idx;
    extern bool mflag;
    int i;
    for (i = 0; i < n; i++) {
	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) {
		i++;
	    }

	    if (i < n)
		inHeader = 0;
	    else
		return;

	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;

	}
	


	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = atgc_to_ry(c[i]);
	//invalid character reset hash
	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;
	}


	if (h->numChar == h->k) {
	    if (h->reverse)
		h->rt_hash -=
		    apnmod[h->hashi[(unsigned char) h->r[h->k - 1]]][0];
	    h->st_hash -= apnmod[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];

	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;


	if (h->reverse) {
	    memmove(&h->r[1], h->r, h->k - 1);
	    h->r[0] = flip(c[i]);
	    h->rt_hash *= IAMOD2P16;
	    h->rt_hash &= MOD2P16;
	    h->rt_hash += apnmod[h->hashi[(unsigned char) c[i]]][h->k - 1];
	    h->rt_hash &= MOD2P16;
	}


	if (h->numChar == h->k) {
	    s = h->s;
	    idx = &h->st_hash;


	    if (h->reverse) {
		if (h->st_hash > h->rt_hash) {
		    s = h->r;
		    idx = &h->rt_hash;
		}
	    }

	    ptr = h->table[*idx];
	    while (ptr != NULL) {
		if (strcmp(ptr->key, s) == 0) {
		    ptr->value++;	// to make generic.
		    break;
		}
		ptr = ptr->next;
	    }
	}
    }
    //should never see
    return;
}



void chkpushatgcw(HASH * h, char *c, int n,bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    static int inHeader = 0;
    unsigned char index;
    char *s;
    unsigned int idx;
    extern bool mflag;
    int i, j;
    for (i = 0; i < n; i++) {
	// It might be possible to buffer this so that we process several bases at once.
	// w/o the expense of a function call every base.
	// avoid the unsigned char casts, just declare as unsigned char
	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) {
		i++;
	    }

	    if (i < n)
		inHeader = 0;
	    else
		return;

	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;

	}



	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = atgc_to_ry(c[i]);
	//invalid character reset hash
	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;
	}


	if (h->numChar == h->k) {
	    h->st_hash -= apnmod[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	    if (h->reverse)
		h->rt_hash -=
		    apnmod[h->hashi[(unsigned char) h->r[h->k - 1]]][0];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];

	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;

	if (h->reverse) {
	    memmove(&h->r[1], h->r, h->k - 1);
	    h->r[0] = flip(c[i]);
	    h->rt_hash *= IAMOD2P16;
	    h->rt_hash &= MOD2P16;
	    h->rt_hash += apnmod[h->hashi[(unsigned char) c[i]]][h->k - 1];
	    h->rt_hash &= MOD2P16;
	}


	if (h->numChar == h->k) {
	    s = h->s;
	    idx = h->st_hash;

	    if (h->reverse) {
		if (h->st_hash > h->rt_hash) {
		    s = h->r;
		    idx = h->rt_hash;
		}
	    }

	    for (j = 0; j < h->k; j++)
		if (weightVector[j] == '0') {
		    idx -= apnmod[h->hashi[(unsigned char) s[j]]][h->k - j - 1];
		    idx &= MOD2P16;
		}


	    ptr = h->table[idx];
	    while (ptr != NULL) {
		if (new_strcmpw(ptr->key, s) == 0) {
		    ptr->value++;	// to make generic.
		    break;
		}
		ptr = ptr->next;
	    }
	}
    }
    //should never see
    return;
}





void pushatgcw(HASH * h, char *c, int n,bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    NODE *r;
    static int inHeader = 0;
    unsigned char index;
    char *s;
    unsigned idx;
    extern bool mflag;
    int i, j;


    for (i = 0; i < n; i++) {


	// It might be possible to buffer this so that we process several bases at once.
	// w/o the expense of a function call every base.
	// avoid the unsigned char casts, just declare as unsigned char
	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) {
		i++;
	    }

	    if (i < n)
		inHeader = 0;
	    else
		return;

	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;
	}

	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = atgc_to_ry(c[i]);

	//invalid character reset hash
	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;
	}



	if (h->numChar == h->k) {
	    h->st_hash -= apnmod[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	    if (h->reverse)
		h->rt_hash -=
		    apnmod[h->hashi[(unsigned char) h->r[h->k - 1]]][0];
	} else
	    h->numChar++;


	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];

	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;



	if (h->reverse) {
	    memmove(&h->r[1], h->r, h->k - 1);
	    h->r[0] = flip(c[i]);
	    h->rt_hash *= IAMOD2P16;
	    h->rt_hash &= MOD2P16;
	    h->rt_hash += apnmod[h->hashi[(unsigned char) c[i]]][h->k - 1];
	    h->rt_hash &= MOD2P16;
	}
	// now perform mask transformation
	if (h->numChar == h->k) {
	    s = h->s;
	    idx = h->st_hash;




	    if (h->reverse) {
		if (h->st_hash > h->rt_hash) {
		    s = h->r;
		    idx = h->rt_hash;
		}
	    }
	    //apply weight mask
	    for (j = 0; j < h->k; j++)
		if (weightVector[j] == '0') {
		    idx -= apnmod[h->hashi[(unsigned char) s[j]]][h->k - j - 1];
		    idx &= MOD2P16;
		}


	    ptr = h->table[idx];

		while (ptr != NULL) {
		    if (new_strcmpw(ptr->key, s) == 0) {
			ptr->value++;	// to make generic.
			break;
		    }
		    ptr = ptr->next;
		}

	    if (ptr == NULL) {
		r = (NODE *) chkmalloc(sizeof(NODE),1);
		r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
		strcpy(r->key, s);
		r->value = 1;
		r->next = h->table[idx];
		h->table[idx] = r;
		h->keyN++;
	    }
	}
    }
    return;
}


/* Established a new convention:  We no longer make lookups to check whether
 * the reverse complement feature is stored in  
 * the hash, we simply do this:  Scanning and hashing of the sequence is done
 * in the forward direction. 
 * We make the decision to physically store a key as in the reverse or forward direction
 * within the hash table: using this criterion:
 * Whichever word would come first if alphabetically sorted is placed in the hash.  
 * Therefore ATG is printed instead of its reverse complement word CAT. With a reverse complement
 * palindrome of course both are equivalent.  By doing this we can half the number
 * of lookups. */

// Use define templates to reduce the overall code size



void pushaa(HASH * h, char *c, int n,bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    NODE *r;
    static int inHeader = 0;
    unsigned char index;
    char *s;
    unsigned *idx;
    extern bool mflag;
    int i;
    for (i = 0; i < n; i++) {
	// It might be possible to buffer this so that we process several bases at once.
	// w/o the expense of a function call every base.
	// avoid the unsigned char casts, just declare as unsigned char
	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) {
		i++;
	    }

	    if (i < n)
		inHeader = 0;
	    else
		return;

	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;

	}



	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = aa_to_class(c[i]);
	//invalid character reset hash

	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = 0;
	    continue;
	}

	c[i] = toupper((int) c[i]);

	if (h->numChar == h->k) {
	    h->st_hash -= apnmodaa[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];



	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;


	if (h->numChar == h->k) {
	    s = h->s;		//no longer needed
	    idx = &h->st_hash;	//ditto

	    ptr = h->table[*idx];
		while (ptr != NULL) {
		    if (strcmp(ptr->key, s) == 0) {
			ptr->value++;	// to make generic.
			break;
		    }
		    ptr = ptr->next;
		}

	    if (ptr == NULL) {
		r = (NODE *) chkmalloc(sizeof(NODE),1);
		r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
		strcpy(r->key, s);
		r->value = 1;
		r->next = h->table[*idx];
		h->table[*idx] = r;
		h->keyN++;
	    }
	}
    }
    return;
}



void chkpushaa(HASH * h, char *c, int n,bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    static int inHeader = 0;
    unsigned char index;
    char *s;
    unsigned *idx;
    extern bool mflag;
    int i;
    for (i = 0; i < n; i++) {
	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) {
		i++;
	    }

	    if (i < n)
		inHeader = 0;
	    else
		return;

	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;

	}



	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = aa_to_class(c[i]);
	//invalid character reset hash
	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = 0;
	    continue;
	}


	if (h->numChar == h->k) {
	    h->st_hash -= apnmodaa[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];

	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;

	if (h->numChar == h->k) {
	    s = h->s;		// remove
	    idx = &h->st_hash;	//ditto

	    ptr = h->table[*idx];
	    while (ptr != NULL) {
		if (strcmp(ptr->key, s) == 0) {
		    ptr->value++;	// to make generic.
		    break;
		}
		ptr = ptr->next;
	    }
	}
    }
    //should never see
    return;
}



void chkpushaaw(HASH * h, char *c, int n, bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    static int inHeader = 0;
    unsigned char index;
    unsigned int idx;
    extern bool mflag;

    int i, j;
    for (i = 0; i < n; i++) {
	// It might be possible to buffer this so that we process several bases at once.
	// w/o the expense of a function call every base.
	// avoid the unsigned char casts, just declare as unsigned char
	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) {
		i++;
	    }

	    if (i < n)
		inHeader = 0;
	    else
		return;

	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;

	}



	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = aa_to_class(c[i]);

	//invalid character reset hash
	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = 0;
	    continue;
	}


	if (h->numChar == h->k) {
	    h->st_hash -= apnmodaa[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];

	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;

	if (h->numChar == h->k) {
	    idx = h->st_hash;
	    for (j = 0; j < h->k; j++)
		if (weightVector[j] == '0') {
		    idx -=
			apnmodaa[h->hashi[(unsigned char) h->s[j]]][h->k - j -
								    1];
		    idx &= MOD2P16;
		}

	    ptr = h->table[idx];
	    while (ptr != NULL) {
		if (new_strcmpw(ptr->key, h->s) == 0) {
		    ptr->value++;	// to make generic.
		    break;
		}
		ptr = ptr->next;
	    }
	}
    }
    return;
}


void chkpushtxt(HASH * h, char *c, int n)
{
    //calculate hash
    NODE *ptr;
    unsigned char index;
    int i;
    for (i = 0; i < n; i++) {



	if (isspace((int) c[i]))
	    continue;

	//invalid character reset hash
	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = 0;
	    continue;
	}


	if (h->numChar == h->k) {
	    h->st_hash -= apnmodtxt[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];

	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;

	if (h->numChar == h->k) {
	    ptr = h->table[h->st_hash];
	    while (ptr != NULL) {
		if (strcmp(ptr->key, h->s) == 0) {
		    ptr->value++;	// to make generic.
		    break;
		}
		ptr = ptr->next;
	    }
	}
    }
    //should never see
    return;
}





void pushaaw(HASH * h, char *c, int n,bool firstRecord)
{
    //calculate hash
    NODE *ptr;
    NODE *r;
    static int inHeader = 0;
    unsigned char index;
    char *s;
    unsigned idx;
    extern bool mflag;
    int i, j;


    for (i = 0; i < n; i++) {

	if (c[i] == '>' || inHeader) {
	    inHeader = 1;
	    // gulp up header
	    while (c[i] != '\n' && i < n) {
		i++;
	    }

	    if (i < n)
		inHeader = 0;
	    else
		return;

	    if (mflag && !firstRecord)  
		printFeatures(h);

	    if (firstRecord) 
		firstRecord=false;

	    h->numChar = h->st_hash = h->rt_hash = 0;
	    continue;

	}

	if (isspace((int) c[i]))
	    continue;

	if (h->isClass)
	    c[i] = aa_to_class(c[i]);

	//invalid character reset hash
	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash;
	    continue;
	}



	if (h->numChar == h->k) {
	    h->st_hash -= apnmodaa[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	} else
	    h->numChar++;


	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];

	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;

	// now perform mask transformation
	if (h->numChar == h->k) {
	    s = h->s;
	    idx = h->st_hash;


	    //apply weight mask
	    for (j = 0; j < h->k; j++)
		if (weightVector[j] == '0') {
		    idx -=
			apnmodaa[h->hashi[(unsigned char) s[j]]][h->k - j - 1];
		    idx &= MOD2P16;
		}


	    ptr = h->table[idx];

		while (ptr != NULL) {
		    if (new_strcmpw(ptr->key, s) == 0) {
			ptr->value++;	// to make generic.
			break;
		    }
		    ptr = ptr->next;
		}

	    if (ptr == NULL) {
		r = (NODE *) chkmalloc(sizeof(NODE),1);
		r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
		strcpy(r->key, s);
		r->value = 1;
		r->next = h->table[idx];
		h->table[idx] = r;
		h->keyN++;
	    }
	}
    }
    return;
}


void pushtxt(HASH * h, char *c, int n)
{
    //calculate hash
    NODE *ptr;
    NODE *r;
    unsigned char index;
    int i;
    for (i = 0; i < n; i++) {

	if (isspace((int) c[i]))
	    continue;

	if (!(index = h->hashi[(unsigned char) c[i]])) {
	    h->numChar = h->st_hash = 0;
	    continue;
	}

	c[i] = toupper((int) c[i]);

	if (h->numChar == h->k) {
	    h->st_hash -= apnmodtxt[h->hashi[(unsigned char) h->s[0]]][h->k - 1];
	} else
	    h->numChar++;

	//push character
	memmove(h->s, &h->s[1], h->k - 1);
	h->s[h->k - 1] = c[i];



	h->st_hash *= A;
	h->st_hash &= MOD2P16;
	h->st_hash += index;
	h->st_hash &= MOD2P16;


	if (h->numChar == h->k) {

	    ptr = h->table[h->st_hash];
		while (ptr != NULL) {
		    if (strcmp(ptr->key,h->s) == 0) {
			ptr->value++;	// to make generic.
			break;
		    }
		    ptr = ptr->next;
		}

	    if (ptr == NULL) {
		r = (NODE *) chkmalloc(sizeof(NODE),1);
		r->key = (char *) chkmalloc(sizeof(char),h->k + 1);
		strcpy(r->key, h->s);
		r->value = 1;
		r->next = h->table[h->st_hash];
		h->table[h->st_hash] = r;
		h->keyN++;
	    }
	}
    }
    return;
}




/**
 *
 * Frees the memory allocated to store the hash
 *
 * All memory allocated with malloc is freed from the
 * hash.  All values are deleted. returns 0 on success.
 * The number of hash keys is also returned to zero.
 *
 * @param None
 * @retval 1 on success
 *
 *
 */


int freeHash(HASH * h)
{
    int i;
    NODE *ptr;
    NODE *last;
    for (i = 0; i < BUCKETS; i++) {
	ptr = h->table[i];
	while (ptr != NULL) {
	    last = ptr;
	    ptr = ptr->next;
	    free(last->key);
	    free(last);
	}
	h->table[i] = NULL;
    }
    h->keyN = 0;
    return 1;
}


/**
 *
 * Masked String comparison method used internally by the hash table functions.
 *
 * This function differs slightly from the library version of strcmp
 * It short circuits when it finds the first difference in strings
 * It is also declared inline so that the compiler has the option of
 * substituting this code in directly.  It is declared static so that
 * it can only be used by the hash functions.
 * Uses a feature mask, ignoring masked out features
 *
 * @param s a null terminated string
 * @param t a null terminated string
 * @retval 1 if s and t are equal
 * @retval 0 if s and t are different
 *
 */

//make static again.
int new_strcmpw(register const char *s, register const char *t)
{
    char *w = weightVector;
    while (*s != '\0')
	if ((*w++) == '1') {
	    if ((*s++) != (*t++))
		return 1;
	} else {
	    s++;
	    t++;
	}
    return 0;
}


/**
 *
 * String comparison method used internally by the hash table functions.
 *
 * This function differs slightly from the library version of strcmp
 * It short circuits when it finds the first difference in strings
 * It is also declared inline so that the compiler has the option of
 * substituting this code in directly.  It is declared static so that
 * it can only be used by the hash functions.
 *
 * @param s a null terminated string
 * @param t a null terminated string
 * @retval 1 if s and t are equal 
 * @ret val 0 if s and t are different
 *
 */





//make static again
int new_strcmp(const char *s, const char *t)
{
    while (*s != '\0') {
	if ((*s++) != (*t++))
	    return 1;
    }
    return 0;
}

/**
 *
 * returns a list of keys stored in the hash table
 *
 * This is function is used to return the values of the
 * keys hashed into the table.  No error checking is
 * provided to make sure that enough memory has been
 * allocated to store all the keys into the array of
 * strings.
 *
 * @param s a pointer to an array of character arrays.
 * @return None
 *
 */


void hashKeys(HASH * h, char ***s)
{
    int i, j;
    NODE *ptr;
    *s = (char **) chkmalloc(sizeof(char *),h->keyN);
    for (i = 0; i < h->keyN; i++)
	(*s)[i] = (char *) chkmalloc(sizeof(char),h->k + 1);


    i = 0;
    for (j = 0; j < BUCKETS; j++)
	if (h->table[j] != NULL) {
	    ptr = h->table[j];
	    while (ptr != NULL) {
		strcpy((*s)[i++], ptr->key);
		//printf("%s\n",ptr->key);
		ptr = ptr->next;
	    }
	}
}


/**
 *
 * returns a list of values stored in the hash table
 *
 * This is function is used to return the values of the
 * keys hashed into the table.  No error checking is
 * provided to make sure that enough memory has been
 * allocated to store all the values into an array
 *
 * @param s a pointer to an array of character arrays.
 * @return None
 * @todo implement error checking.
 */


void hashKeysAndValues(HASH * h, char ***s, unsigned **d)
{
    unsigned i, j;
    NODE *ptr;
    char **sp;
    unsigned * dd;

    dd = (unsigned *) chkmalloc(sizeof(unsigned),h->keyN);
    sp = (char **) chkmalloc(sizeof(char *), h->keyN);
    for (i = 0; i < h->keyN; i++) {
	sp[i] = (char *) malloc(sizeof(char)*(h->k + 1));
	}
    i = 0;
    for (j = 0; j < BUCKETS; j++)
	if (h->table[j] != NULL) {
	    ptr = h->table[j];
	    while (ptr != NULL) {
		strcpy(sp[i], ptr->key);
		dd[i++] = ptr->value;
		ptr = ptr->next;
	    }
	}
    *s=sp;
    *d=dd;
}


// Returns sum of all values stored in hash


unsigned int sumValues(HASH * h)
{
    int j;
    NODE *ptr;
    unsigned sum = 0;

    for (j = 0; j < BUCKETS; j++)
	if (h->table[j] != NULL) {
	    ptr = h->table[j];
	    while (ptr != NULL) {
		sum += ptr->value;
		ptr = ptr->next;
	    }
	}
    return sum;
}




/** Extract hash values to an array and i
 *  set hash values to zero.
 */


void hashValuesAndSet(HASH * h, unsigned **d)
{
    int i, j;
    NODE *ptr;
    *d = (unsigned *) chkmalloc(sizeof(unsigned),h->keyN);

    i = 0;
    for (j = 0; j < BUCKETS; j++)
	if (h->table[j] != NULL) {
	    ptr = h->table[j];
	    while (ptr != NULL) {
		(*d)[i++] = ptr->value;
		ptr->value = 0;
		ptr = ptr->next;
	    }
	}
}




/* Prints featuers, resets and frees hash memory */


void printFeatures(HASH * h)
{
    int i;
    char **keys;
    unsigned *values;
    extern char fflag;

    if (h->keyN == 0)
	warn_msg("Warning: No keys of length %d found.\n", h->k);

    if (fflag)
	hashValuesAndSet(h, &values);
    else
	hashKeysAndValues(h, &keys, &values);


    for (i = 0; i < h->keyN-1; i++) {
	if (fflag)
	    printf("%d\t", values[i]);
	else {
	    printf("%s\t%d\t", keys[i], values[i]);
	    free(keys[i]);
	}
    }
	if (fflag)
	    printf("%d\n", values[i]);
	else {
	    printf("%s\t%d\n", keys[i], values[i]);
	    free(keys[i]);
	}

    resetHash(h);
    if (!fflag)
	freeHash(h);
    free(values);
}
