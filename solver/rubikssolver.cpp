

#include <string.h>
#include <stdlib.h>
#include <Slice_Flip_Prun.h>
#include <Slice_Twist_Prun.h>
#include <flipMove.h>
#include <FRtoBR_Move.h>
#include <MergeURtoULandUBtoDF2.h>
#include <Slice_URFtoDLF_Parity_Prun.h>
#include <Slice_URtoDF_Parity_Prun.h>
#include <twistMove.h>
#include <UBtoDF_Move.h>

#include <URtoUL_Move.h>
#include <time.h>


#define N_TWIST     2187
#define N_FLIP      2048
#define N_SLICE1    495
#define N_SLICE2    24
#define N_PARITY    2
#define N_URFtoDLF  20160
#define N_FRtoBR    11880
#define N_URtoUL    1320
#define N_UBtoDF    1320
#define N_URtoDF    20160
#define N_URFtoDLB  40320
#define N_URtoBR    479001600
#define N_MOVE      18

#define CORNER_COUNT 8
#define EDGE_COUNT 12

typedef enum {
    URF, UFL, ULB, UBR, DFR, DLF, DBL, DRB
} corner_t;

typedef enum {
    UR, UF, UL, UB, DR, DF, DL, DB, FR, FL, BL, BR
} edge_t;

struct cubiecube {
    // initialize to Id-Cube
    // corner permutation
    corner_t cp[8];
    // corner orientation
    signed char co[8];
    // edge permutation
    edge_t ep[12];
    // edge orientation
    signed char eo[12];
};
typedef struct cubiecube cubiecube_t;

cubiecube_t* get_cubiecube()
{
    cubiecube_t* result = (cubiecube_t *) calloc(1, sizeof(cubiecube_t));

    static const corner_t   cp[8]   = { URF, UFL, ULB, UBR, DFR, DLF, DBL, DRB };
    static const signed char       co[8]   = { 0, 0, 0, 0, 0, 0, 0, 0 };
    static const edge_t     ep[12]  = { UR, UF, UL, UB, DR, DF, DL, DB, FR, FL, BL, BR };
    static const signed char       eo[12]  = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    memcpy(result->cp, cp, sizeof(cp));
    memcpy(result->co, co, sizeof(co));
    memcpy(result->ep, ep, sizeof(ep));
    memcpy(result->eo, eo, sizeof(eo));

    return result;
}

cubiecube_t * get_moveCube()
{
    static cubiecube_t moveCube[6];
    static int moveCube_initialized = 0;
    static const corner_t     cpU[8]  = { UBR, URF, UFL, ULB, DFR, DLF, DBL, DRB };
    static const signed char  coU[8]  = { 0, 0, 0, 0, 0, 0, 0, 0 };
    static const edge_t       epU[12] = { UB, UR, UF, UL, DR, DF, DL, DB, FR, FL, BL, BR };
    static const signed char  eoU[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    static const corner_t     cpR[8]  = { DFR, UFL, ULB, URF, DRB, DLF, DBL, UBR };
    static const signed char  coR[8]  = { 2, 0, 0, 1, 1, 0, 0, 2 };
    static const edge_t       epR[12] = { FR, UF, UL, UB, BR, DF, DL, DB, DR, FL, BL, UR };
    static const signed char  eoR[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    static const corner_t     cpF[8]  = { UFL, DLF, ULB, UBR, URF, DFR, DBL, DRB };
    static const signed char  coF[8]  = { 1, 2, 0, 0, 2, 1, 0, 0 };
    static const edge_t       epF[12] = { UR, FL, UL, UB, DR, FR, DL, DB, UF, DF, BL, BR };
    static const signed char  eoF[12] = { 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0 };
    static const corner_t     cpD[8]  = { URF, UFL, ULB, UBR, DLF, DBL, DRB, DFR };
    static const signed char  coD[8]  = { 0, 0, 0, 0, 0, 0, 0, 0 };
    static const edge_t       epD[12] = { UR, UF, UL, UB, DF, DL, DB, DR, FR, FL, BL, BR };
    static const signed char  eoD[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    static const corner_t     cpL[8]  = { URF, ULB, DBL, UBR, DFR, UFL, DLF, DRB };
    static const signed char  coL[8]  = { 0, 1, 2, 0, 0, 2, 1, 0 };
    static const edge_t       epL[12] = { UR, UF, BL, UB, DR, DF, FL, DB, FR, UL, DL, BR };
    static const signed char  eoL[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    static const corner_t     cpB[8]  = { URF, UFL, UBR, DRB, DFR, DLF, ULB, DBL };
    static const signed char  coB[8]  = { 0, 0, 1, 2, 0, 0, 2, 1 };
    static const edge_t       epB[12] = { UR, UF, UL, BR, DR, DF, DL, BL, FR, FL, UB, DB };
    static const signed char  eoB[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1 };

    if (!moveCube_initialized) {
        memcpy(moveCube[0].cp, cpU, sizeof(cpU));
        memcpy(moveCube[0].co, coU, sizeof(coU));
        memcpy(moveCube[0].ep, epU, sizeof(epU));
        memcpy(moveCube[0].eo, eoU, sizeof(eoU));
        memcpy(moveCube[1].cp, cpR, sizeof(cpR));
        memcpy(moveCube[1].co, coR, sizeof(coR));
        memcpy(moveCube[1].ep, epR, sizeof(epR));
        memcpy(moveCube[1].eo, eoR, sizeof(eoR));
        memcpy(moveCube[2].cp, cpF, sizeof(cpF));
        memcpy(moveCube[2].co, coF, sizeof(coF));
        memcpy(moveCube[2].ep, epF, sizeof(epF));
        memcpy(moveCube[2].eo, eoF, sizeof(eoF));
        memcpy(moveCube[3].cp, cpD, sizeof(cpD));
        memcpy(moveCube[3].co, coD, sizeof(coD));
        memcpy(moveCube[3].ep, epD, sizeof(epD));
        memcpy(moveCube[3].eo, eoD, sizeof(eoD));
        memcpy(moveCube[4].cp, cpL, sizeof(cpL));
        memcpy(moveCube[4].co, coL, sizeof(coL));
        memcpy(moveCube[4].ep, epL, sizeof(epL));
        memcpy(moveCube[4].eo, eoL, sizeof(eoL));
        memcpy(moveCube[5].cp, cpB, sizeof(cpB));
        memcpy(moveCube[5].co, coB, sizeof(coB));
        memcpy(moveCube[5].ep, epB, sizeof(epB));
        memcpy(moveCube[5].eo, eoB, sizeof(eoB));
    }

    return moveCube;
}

void rotateRight_corner(corner_t* arr, int l, int r)
// Right rotation of all array elements between l and r
{
    int i;
    corner_t temp = arr[r];
    for (i = r; i > l; i--)
        arr[i] = arr[i - 1];
    arr[l] = temp;
}

int Cnk(int n, int k) {
    int i, j, s;
    if (n < k)
        return 0;
    if (k > n / 2)
        k = n - k;
    for (s = 1, i = n, j = 1; i != n - k; i--, j++) {
        s *= i;
        s /= j;
    }
    return s;
}

void setURFtoDLF(cubiecube_t* cubiecube, short idx)
{
    int x;
    corner_t corner6[6] = { URF, UFL, ULB, UBR, DFR, DLF };
    corner_t otherCorner[2] = { DBL, DRB };
    int b = idx % 720; // Permutation
    int a = idx / 720; // Combination
    int c, j, k;
    for(c = 0; c < CORNER_COUNT; c++)
        cubiecube->cp[c] = DRB;// Use DRB to invalidate all corners

    for (j = 1; j < 6; j++)// generate permutation from index b
    {
        k = b % (j + 1);
        b /= j + 1;
        while (k-- > 0)
            rotateRight_corner(corner6, 0, j);
    }
    x = 5;// generate combination and set corners
    for (j = DRB; j >= 0; j--)
        if (a - Cnk(j, x + 1) >= 0) {
            cubiecube->cp[j] = corner6[x];
            a -= Cnk(j, x-- + 1);
        }
    x = 0;
    for (j = URF; j <= DRB; j++)
        if (cubiecube->cp[j] == DRB)
            cubiecube->cp[j] = otherCorner[x++];
}

void cornerMultiply(cubiecube_t* cubiecube, cubiecube_t* b)
{
    int corn;
    signed char oriA, oriB, ori;
    corner_t cPerm[8] = {URF};
    signed char cOri[8] = {0};
    for (corn = 0; corn < CORNER_COUNT; corn++) {
        cPerm[corn] = cubiecube->cp[b->cp[corn]];

        oriA = cubiecube->co[b->cp[corn]];
        oriB = b->co[corn];
        ori = 0;

        if (oriA < 3 && oriB < 3) // if both cubes are regular cubes...
        {
            ori = oriA + oriB; // just do an addition modulo 3 here
            if (ori >= 3)
                ori -= 3; // the composition is a regular cube

            // +++++++++++++++++++++not used in this implementation +++++++++++++++++++++++++++++++++++
        } else if (oriA < 3 && oriB >= 3) // if cube b is in a mirrored
            // state...
        {
            ori = oriA + oriB;
            if (ori >= 6)
                ori -= 3; // the composition is a mirrored cube
        } else if (oriA >= 3 && oriB < 3) // if cube a is an a mirrored
            // state...
        {
            ori = oriA - oriB;
            if (ori < 3)
                ori += 3; // the composition is a mirrored cube
        } else if (oriA >= 3 && oriB >= 3) // if both cubes are in mirrored
            // states...
        {
            ori = oriA - oriB;
            if (ori < 0)
                ori += 3; // the composition is a regular cube
            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        }
        cOri[corn] = ori;
    }
    for(corn = 0; corn < CORNER_COUNT; corn++) {
        cubiecube->cp[corn] = cPerm[corn];
        cubiecube->co[corn] = cOri[corn];
    }
}

void rotateLeft_corner(corner_t* arr, int l, int r)
// Left rotation of all array elements between l and r
{
    int i;
    corner_t temp = arr[l];
    for (i = l; i < r; i++)
        arr[i] = arr[i + 1];
    arr[r] = temp;
}

short getURFtoDLF(cubiecube_t* cubiecube)
{
    int a = 0, x = 0, j, b = 0;
    corner_t corner6[6] = {URF};
    // compute the index a < (8 choose 6) and the corner permutation.
    for (j = URF; j <= DRB; j++)
        if (cubiecube->cp[j] <= DLF) {
            a += Cnk(j, x + 1);
            corner6[x++] = cubiecube->cp[j];
        }

    for (j = 5; j > 0; j--)// compute the index b < 6! for the
        // permutation in corner6
    {
        int k = 0;
        while (corner6[j] != j) {
            rotateLeft_corner(corner6, 0, j);
            k++;
        }
        b = (j + 1) * b + k;
    }
    return (short) (720 * a + b);
}

void rotateRight_edge(edge_t* arr, int l, int r)
// Right rotation of all array elements between l and r
{
    int i;
    edge_t temp = arr[r];
    for (i = r; i > l; i--)
        arr[i] = arr[i - 1];
    arr[l] = temp;
}

void setURtoDF(cubiecube_t* cubiecube, int idx)
{
    int x, e, j, k;
    edge_t edge6[6] = { UR, UF, UL, UB, DR, DF };
    edge_t otherEdge[6] = { DL, DB, FR, FL, BL, BR };
    int b = idx % 720; // Permutation
    int a = idx / 720; // Combination

    for(e = 0; e < EDGE_COUNT; e++)
        cubiecube->ep[e] = BR;// Use BR to invalidate all edges

    for (j = 1; j < 6; j++)// generate permutation from index b
    {
        k = b % (j + 1);
        b /= j + 1;
        while (k-- > 0)
            rotateRight_edge(edge6, 0, j);
    }
    x = 5;// generate combination and set edges
    for (j = BR; j >= 0; j--)
        if (a - Cnk(j, x + 1) >= 0) {
            cubiecube->ep[j] = edge6[x];
            a -= Cnk(j, x-- + 1);
        }
    x = 0; // set the remaining edges DL..BR
    for (j = UR; j <= BR; j++)
        if (cubiecube->ep[j] == BR)
            cubiecube->ep[j] = otherEdge[x++];
}

void edgeMultiply(cubiecube_t* cubiecube, cubiecube_t* b)
{
    int edge;
    edge_t ePerm[12] = {UR};
    signed char eOri[12] = {0};

    for(edge = 0; edge < EDGE_COUNT; edge++) {
        ePerm[edge] = cubiecube->ep[b->ep[edge]];
        eOri[edge] = (b->eo[edge] + cubiecube->eo[b->ep[edge]]) % 2;
    }
    for(edge = 0; edge < EDGE_COUNT; edge++) {
        cubiecube->ep[edge] = ePerm[edge];
        cubiecube->eo[edge] = eOri[edge];
    }
}

void rotateLeft_edge(edge_t* arr, int l, int r)
// Left rotation of all array elements between l and r
{
    int i;
    edge_t temp = arr[l];
    for (i = l; i < r; i++)
        arr[i] = arr[i + 1];
    arr[r] = temp;
}

int getURtoDF(cubiecube_t* cubiecube)
{
    int a = 0, x = 0;
    int b = 0, j;
    edge_t edge6[6] = {UR};
    // compute the index a < (12 choose 6) and the edge permutation.
    for (j = UR; j <= BR; j++)
        if (cubiecube->ep[j] <= DF) {
            a += Cnk(j, x + 1);
            edge6[x++] = cubiecube->ep[j];
        }

    for (j = 5; j > 0; j--)// compute the index b < 6! for the
        // permutation in edge6
    {
        int k = 0;
        while (edge6[j] != j) {
            rotateLeft_edge(edge6, 0, j);
            k++;
        }
        b = (j + 1) * b + k;
    }
    return 720 * a + b;
}

void setPruning(signed char *table, int index, signed char value) {
    if ((index & 1) == 0)
        table[index / 2] &= 0xf0 | value;
    else
        table[index / 2] &= 0x0f | (value << 4);
}

signed char getPruning(signed char *table, int index) {
    signed char res;

    if ((index & 1) == 0)
        res = (table[index / 2] & 0x0f);
    else
        res = ((table[index / 2] >> 4) & 0x0f);

    return res;
}

typedef enum {U, R, F, D, L, B} color_t;

typedef struct {
    int ax[31];             // The axis of the move
    int po[31];             // The power of the move
    int flip[31];           // phase1 coordinates
    int twist[31];
    int slice[31];
    int parity[31];         // phase2 coordinates
    int URFtoDLF[31];
    int FRtoBR[31];
    int URtoUL[31];
    int UBtoDF[31];
    int URtoDF[31];
    int minDistPhase1[31];  // IDA* distance do goal estimations
    int minDistPhase2[31];
} search_t;


struct facecube {
    color_t f[54];
};
typedef struct facecube facecube_t;

typedef struct {

    // All coordinates are 0 for a solved cube except for UBtoDF, which is 114
    short twist;
    short flip;
    short parity;
    short FRtoBR;
    short URFtoDLF;
    short URtoUL;
    short UBtoDF;
    int URtoDF;
} coordcube_t;

facecube_t* get_facecube_fromstring(int* cube_colors)
{
    int i;
    facecube_t* res = (facecube_t *) calloc(1, sizeof(facecube_t));
    for (i = 0; i < 54; i++) {
        switch(cube_colors[i]) {
            case 0:
                res->f[i] = U;
                break;
            case 1:
                res->f[i] = R;
                break;
            case 2:
                res->f[i] = F;
                break;
            case 3:
                res->f[i] = D;
                break;
            case 4:
                res->f[i] = L;
                break;
            case 5:
                res->f[i] = B;
                break;
        }
    }
    return res;
}

typedef enum {
    U1, U2, U3, U4, U5, U6, U7, U8, U9, R1, R2, R3, R4, R5, R6, R7, R8, R9, F1, F2, F3, F4, F5, F6, F7, F8, F9, D1, D2, D3, D4, D5, D6, D7, D8, D9, L1, L2, L3, L4, L5, L6, L7, L8, L9, B1, B2, B3, B4, B5, B6, B7, B8, B9
} facelet_t;

#define FACELET_COUNT 54

facelet_t cornerFacelet[8][3] = { { U9, R1, F3 }, { U7, F1, L3 }, { U1, L1, B3 }, { U3, B1, R3 },
                                  { D3, F9, R7 }, { D1, L9, F7 }, { D7, B9, L7 }, { D9, R9, B7 } };

facelet_t edgeFacelet[12][2] = { { U6, R2 }, { U8, F2 }, { U4, L2 }, { U2, B2 }, { D6, R8 }, { D2, F8 },
                                 { D4, L8 }, { D8, B8 }, { F6, R4 }, { F4, L6 }, { B6, L4 }, { B4, R6 } };

color_t cornerColor[8][3] = { { U, R, F }, { U, F, L }, { U, L, B }, { U, B, R }, { D, F, R }, { D, L, F },
                              { D, B, L }, { D, R, B } };

color_t edgeColor[12][2] = { { U, R }, { U, F }, { U, L }, { U, B }, { D, R }, { D, F }, { D, L }, { D, B },
                             { F, R }, { F, L }, { B, L }, { B, R } };


cubiecube_t* toCubieCube(facecube_t* facecube)
{
    int i, j;
    signed char ori;
    color_t col1, col2;
    cubiecube_t* ccRet = (cubiecube_t*) calloc(1, sizeof(cubiecube_t));
    for (i = 0; i < 8; i++)
        ccRet->cp[i] = URF;// invalidate corners
    for (i = 0; i < 12; i++)
        ccRet->ep[i] = UR;// and edges

    for(i = 0; i < CORNER_COUNT; i++) {
        // get the colors of the cubie at corner i, starting with U/D
        for (ori = 0; ori < 3; ori++)
            if (facecube->f[cornerFacelet[i][ori]] == U || facecube->f[cornerFacelet[i][ori]] == D)
                break;
        col1 = facecube->f[cornerFacelet[i][(ori + 1) % 3]];
        col2 = facecube->f[cornerFacelet[i][(ori + 2) % 3]];

        for (j = 0; j < CORNER_COUNT; j++) {
            if (col1 == cornerColor[j][1] && col2 == cornerColor[j][2]) {
                // in cornerposition i we have cornercubie j
                ccRet->cp[(corner_t)i] = (corner_t)j;
                ccRet->co[i] = ori % 3;
                break;
            }
        }
    }

    for (i = 0; i < EDGE_COUNT; i++) {
        for (j = 0; j < EDGE_COUNT; j++) {
            if (facecube->f[edgeFacelet[i][0]] == edgeColor[j][0]
                && facecube->f[edgeFacelet[i][1]] == edgeColor[j][1]) {
                ccRet->ep[(edge_t)i] = (edge_t)j;
                ccRet->eo[i] = 0;
                break;
            }
            if (facecube->f[edgeFacelet[i][0]] == edgeColor[j][1]
                && facecube->f[edgeFacelet[i][1]] == edgeColor[j][0]) {
                ccRet->ep[(edge_t)i] = (edge_t)j;
                ccRet->eo[i] = 1;
                break;
            }
        }
    }
    return ccRet;
}

short cornerParity(cubiecube_t* cubiecube)
{
    int i, j;
    int s = 0;
    for (i = (int)DRB; i >= (int)URF + 1; i--)
        for (j = i - 1; j >= (int)URF; j--)
            if (cubiecube->cp[j] > cubiecube->cp[i])
                s++;
    return (short) (s % 2);
}

short edgeParity(cubiecube_t* cubiecube)
{
    int i, j;
    int s = 0;
    for (i = (int)BR; i >= (int)UR + 1; i--)
        for (j = i - 1; j >= (int)UR; j--)
            if (cubiecube->ep[j] > cubiecube->ep[i])
                s++;
    return (short) (s % 2);
}

int verify(cubiecube_t* cubiecube)
{
    int sum = 0, e, i, c;
    int edgeCount[12] = {0};
    int cornerCount[8] = {0};

    for(e = 0; e < EDGE_COUNT; e++)
        edgeCount[cubiecube->ep[e]]++;
    for (i = 0; i < 12; i++)
        if (edgeCount[i] != 1)
            return -2;

    for (i = 0; i < 12; i++)
        sum += cubiecube->eo[i];
    if (sum % 2 != 0)
        return -3;

    for(c = 0; c < CORNER_COUNT; c++)
        cornerCount[cubiecube->cp[c]]++;
    for (i = 0; i < 8; i++)
        if (cornerCount[i] != 1)
            return -4;// missing corners

    sum = 0;
    for (i = 0; i < 8; i++)
        sum += cubiecube->co[i];
    if (sum % 3 != 0)
        return -5;// twisted corner

    if ((edgeParity(cubiecube) ^ cornerParity(cubiecube)) != 0)
        return -6;// parity error

    return 0;// cube ok
}

short getTwist(cubiecube_t* cubiecube)
{
    short ret = 0;
    int i;
    for (i = (int)URF; i < (int)DRB; i++)
        ret = (short) (3 * ret + cubiecube->co[i]);
    return ret;
}

short getFlip(cubiecube_t* cubiecube)
{
    int i;
    short ret = 0;
    for (i = (int)UR; i < (int)BR; i++)
        ret = (short) (2 * ret + cubiecube->eo[i]);
    return ret;
}

short getFRtoBR(cubiecube_t* cubiecube)
{
    int a = 0, x = 0, j;
    int b = 0;
    edge_t edge4[4] = {(edge_t)0};
    // compute the index a < (12 choose 4) and the permutation array perm.
    for (j = (int)BR; j >= (int)UR; j--)
        if ((int)FR <= cubiecube->ep[j] && cubiecube->ep[j] <= (int)BR) {
            a += Cnk(11 - j, x + 1);
            edge4[3 - x++] = cubiecube->ep[j];
        }

    for (j = 3; j > 0; j--)// compute the index b < 4! for the
        // permutation in perm
    {
        int k = 0;
        while (edge4[j] != j + 8) {
            rotateLeft_edge(edge4, 0, j);
            k++;
        }
        b = (j + 1) * b + k;
    }
    return (short) (24 * a + b);
}

short getURtoUL(cubiecube_t* cubiecube)
{
    int a = 0, b = 0, x = 0, j;
    edge_t edge3[3] = {(edge_t)0};
    // compute the index a < (12 choose 3) and the edge permutation.
    for (j = (int)UR; j <= (int)BR; j++)
        if (cubiecube->ep[j] <= (int)UL) {
            a += Cnk(j, x + 1);
            edge3[x++] = cubiecube->ep[j];
        }

    for (j = 2; j > 0; j--)// compute the index b < 3! for the
        // permutation in edge3
    {
        int k = 0;
        while (edge3[j] != j) {
            rotateLeft_edge(edge3, 0, j);
            k++;
        }
        b = (j + 1) * b + k;
    }
    return (short) (6 * a + b);
}

short getUBtoDF(cubiecube_t* cubiecube)
{
    int a = 0, x = 0, b = 0, j;
    edge_t edge3[3] = {(edge_t)0};
    // compute the index a < (12 choose 3) and the edge permutation.
    for (j = (int)UR; j <= (int)BR; j++)
        if ((int)UB <= cubiecube->ep[j] && cubiecube->ep[j] <= (int)DF) {
            a += Cnk(j, x + 1);
            edge3[x++] = cubiecube->ep[j];
        }

    for (j = 2; j > 0; j--)// compute the index b < 3! for the
        // permutation in edge3
    {
        int k = 0;
        while (edge3[j] != (int)UB + j) {
            rotateLeft_edge(edge3, 0, j);
            k++;
        }
        b = (j + 1) * b + k;
    }
    return (short) (6 * a + b);
}

coordcube_t* get_coordcube(cubiecube_t* cubiecube)
{
    coordcube_t* result = (coordcube_t *) calloc(1, sizeof(coordcube_t));

    result->twist       = getTwist(cubiecube);
    result->flip        = getFlip(cubiecube);
    result->parity      = cornerParity(cubiecube);
    result->FRtoBR      = getFRtoBR(cubiecube);
    result->URFtoDLF    = getURFtoDLF(cubiecube);
    result->URtoUL      = getURtoUL(cubiecube);
    result->UBtoDF      = getUBtoDF(cubiecube);
    result->URtoDF      = getURtoDF(cubiecube);// only needed in phase2

    return result;
}

#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))



short parityMove[2][18] = {
        { 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1 },
        { 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 }
};

short URFtoDLF_Move[N_URFtoDLF][N_MOVE] = {{0}};
short URtoDF_Move[N_URtoDF][N_MOVE] = {{0}};

int totalDepth(search_t* search, int depthPhase1, int maxDepth)
{
    int mv = 0, d1 = 0, d2 = 0, i;
    int maxDepthPhase2 = MIN(10, maxDepth - depthPhase1);// Allow only max 10 moves in phase2
    int depthPhase2;
    int n;
    int busy;
    for (i = 0; i < depthPhase1; i++) {
        mv = 3 * search->ax[i] + search->po[i] - 1;
        // System.out.format("%d %d %d %d\n", i, mv, ax[i], po[i]);
        search->URFtoDLF[i + 1] = URFtoDLF_Move[search->URFtoDLF[i]][mv];
        search->FRtoBR[i + 1] = FRtoBR_Move[search->FRtoBR[i]][mv];
        search->parity[i + 1] = parityMove[search->parity[i]][mv];
    }

    if ((d1 = getPruning(Slice_URFtoDLF_Parity_Prun,
                         (N_SLICE2 * search->URFtoDLF[depthPhase1] + search->FRtoBR[depthPhase1]) * 2 + search->parity[depthPhase1])) > maxDepthPhase2)
        return -1;

    for (i = 0; i < depthPhase1; i++) {
        mv = 3 * search->ax[i] + search->po[i] - 1;
        search->URtoUL[i + 1] = URtoUL_Move[search->URtoUL[i]][mv];
        search->UBtoDF[i + 1] = UBtoDF_Move[search->UBtoDF[i]][mv];
    }
    search->URtoDF[depthPhase1] = MergeURtoULandUBtoDF[search->URtoUL[depthPhase1]][search->UBtoDF[depthPhase1]];

    if ((d2 = getPruning(Slice_URtoDF_Parity_Prun,
                         (N_SLICE2 * search->URtoDF[depthPhase1] + search->FRtoBR[depthPhase1]) * 2 + search->parity[depthPhase1])) > maxDepthPhase2)
        return -1;

    if ((search->minDistPhase2[depthPhase1] = MAX(d1, d2)) == 0)// already solved
        return depthPhase1;

    // now set up search

    depthPhase2 = 1;
    n = depthPhase1;
    busy = 0;
    search->po[depthPhase1] = 0;
    search->ax[depthPhase1] = 0;
    search->minDistPhase2[n + 1] = 1;// else failure for depthPhase2=1, n=0
    // +++++++++++++++++++ end initialization +++++++++++++++++++++++++++++++++
    do {
        do {
            if ((depthPhase1 + depthPhase2 - n > search->minDistPhase2[n + 1]) && !busy) {

                if (search->ax[n] == 0 || search->ax[n] == 3)// Initialize next move
                {
                    search->ax[++n] = 1;
                    search->po[n] = 2;
                } else {
                    search->ax[++n] = 0;
                    search->po[n] = 1;
                }
            } else if ((search->ax[n] == 0 || search->ax[n] == 3) ? (++search->po[n] > 3) : ((search->po[n] = search->po[n] + 2) > 3)) {
                do {// increment axis
                    if (++search->ax[n] > 5) {
                        if (n == depthPhase1) {
                            if (depthPhase2 >= maxDepthPhase2)
                                return -1;
                            else {
                                depthPhase2++;
                                search->ax[n] = 0;
                                search->po[n] = 1;
                                busy = 0;
                                break;
                            }
                        } else {
                            n--;
                            busy = 1;
                            break;
                        }

                    } else {
                        if (search->ax[n] == 0 || search->ax[n] == 3)
                            search->po[n] = 1;
                        else
                            search->po[n] = 2;
                        busy = 0;
                    }
                } while (n != depthPhase1 && (search->ax[n - 1] == search->ax[n] || search->ax[n - 1] - 3 == search->ax[n]));
            } else
                busy = 0;
        } while (busy);
        // +++++++++++++ compute new coordinates and new minDist ++++++++++
        mv = 3 * search->ax[n] + search->po[n] - 1;

        search->URFtoDLF[n + 1] = URFtoDLF_Move[search->URFtoDLF[n]][mv];
        search->FRtoBR[n + 1] = FRtoBR_Move[search->FRtoBR[n]][mv];
        search->parity[n + 1] = parityMove[search->parity[n]][mv];
        search->URtoDF[n + 1] = URtoDF_Move[search->URtoDF[n]][mv];

        search->minDistPhase2[n + 1] = MAX(getPruning(Slice_URtoDF_Parity_Prun, (N_SLICE2
                                                                                 * search->URtoDF[n + 1] + search->FRtoBR[n + 1])
                                                                                * 2 + search->parity[n + 1]), getPruning(Slice_URFtoDLF_Parity_Prun, (N_SLICE2
                                                                                                                                                      * search->URFtoDLF[n + 1] + search->FRtoBR[n + 1])
                                                                                                                                                     * 2 + search->parity[n + 1]));
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    } while (search->minDistPhase2[n + 1] != 0);
    return depthPhase1 + depthPhase2;
}

int rotations[70] = {0};

char* solutionToString(search_t* search, int length, int depthPhase1)
{
    char* s = (char*) calloc(length * 3 + 5, 1);
    int cur = 0, i;
    int j = 0;
    for (i = 0; i < length; i++) {
        switch (search->ax[i]) {
            case 0:
                s[cur++] = 'U';
                rotations[j] = 6;
                break;
            case 1:
                s[cur++] = 'R';
                rotations[j] = 3;
                break;
            case 2:
                s[cur++] = 'F';
                rotations[j] = 9;
                break;
            case 3:
                s[cur++] = 'D';
                rotations[j] = -4;
                break;
            case 4:
                s[cur++] = 'L';
                rotations[j] = -1;

                break;
            case 5:
                s[cur++] = 'B';
                rotations[j] = -7;

                break;
        }
        switch (search->po[i]) {
            case 1:
                s[cur++] = ' ';
                break;
            case 2:
                j++;
                rotations[j] = rotations[j-1];
                s[cur++] = '2';
                s[cur++] = ' ';
                break;
            case 3:
                rotations[j] = -rotations[j];
                s[cur++] = '\'';
                s[cur++] = ' ';
                break;
        }
        if (i == depthPhase1 - 1) {
            s[cur++] = '.';
            s[cur++] = ' ';
        }
        j++;
    }
    return s;
}

int* solvecube(int* cube_colors, long ptimeout, int pmaxdepth)
{
	long timeOut = ptimeout;
    int maxDepth = pmaxdepth;
	
    search_t* search = (search_t*) calloc(1, sizeof(search_t));
    facecube_t* fc;
    cubiecube_t* cc;
    coordcube_t* c;

    int s, i;
    int mv, n;
    int busy;
    int depthPhase1;
    time_t tStart;

    int int_cube_colors[54];



    for(int i = 0;i<70;++i){
        rotations[i] = 0;
    }

    for(int i=0;i<54;++i){
        int_cube_colors[i] = cube_colors[i];
    }

    fc = get_facecube_fromstring(int_cube_colors);
    cc = toCubieCube(fc);
    if ((s = verify(cc)) != 0) {
        free(search);
        return nullptr;
    }

    c = get_coordcube(cc);

    search->po[0] = 0;
    search->ax[0] = 0;
    search->flip[0] = c->flip;
    search->twist[0] = c->twist;
    search->parity[0] = c->parity;
    search->slice[0] = c->FRtoBR / 24;
    search->URFtoDLF[0] = c->URFtoDLF;
    search->FRtoBR[0] = c->FRtoBR;
    search->URtoUL[0] = c->URtoUL;
    search->UBtoDF[0] = c->UBtoDF;

    search->minDistPhase1[1] = 1;// else failure for depth=1, n=0
    mv = 0;
    n = 0;
    busy = 0;
    depthPhase1 = 1;

    int useSeparator = 0;
    char* res;

    bool stopSearch = false;

    tStart = time(NULL);

    do {
        do {
            if ((depthPhase1 - n > search->minDistPhase1[n + 1]) && !busy) {

                if (search->ax[n] == 0 || search->ax[n] == 3)// Initialize next move
                    search->ax[++n] = 1;
                else
                    search->ax[++n] = 0;
                search->po[n] = 1;
            } else if (++search->po[n] > 3) {
                do {// increment axis
                    if (++search->ax[n] > 5) {

                        if (time(NULL) - tStart > timeOut)
                            return NULL;

                        if (n == 0) {
                            if (depthPhase1 >= maxDepth)
                                return NULL;
                            else {
                                depthPhase1++;
                                search->ax[n] = 0;
                                search->po[n] = 1;
                                busy = 0;
                                break;
                            }
                        } else {
                            n--;
                            busy = 1;
                            break;
                        }

                    } else {
                        search->po[n] = 1;
                        busy = 0;
                    }
                } while (n != 0 && (search->ax[n - 1] == search->ax[n] || search->ax[n - 1] - 3 == search->ax[n]));
            } else
                busy = 0;
        } while (busy);

        // +++++++++++++ compute new coordinates and new minDistPhase1 ++++++++++
        // if minDistPhase1 =0, the H subgroup is reached
        mv = 3 * search->ax[n] + search->po[n] - 1;
        search->flip[n + 1] = flipMove[search->flip[n]][mv];
        search->twist[n + 1] = twistMove[search->twist[n]][mv];
        search->slice[n + 1] = FRtoBR_Move[search->slice[n] * 24][mv] / 24;
        search->minDistPhase1[n + 1] = MAX(
                getPruning(Slice_Flip_Prun, N_SLICE1 * search->flip[n + 1] + search->slice[n + 1]),
                getPruning(Slice_Twist_Prun, N_SLICE1 * search->twist[n + 1] + search->slice[n + 1])
        );
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // System.out.format("%d %d\n", n, depthPhase1);
        if (search->minDistPhase1[n + 1] == 0 && n >= depthPhase1 - 5) {
            search->minDistPhase1[n + 1] = 10;// instead of 10 any value >5 is possible
            if (n == depthPhase1 - 1 && (s = totalDepth(search, depthPhase1, maxDepth)) >= 0) {
                if (s == depthPhase1
                    || (search->ax[depthPhase1 - 1] != search->ax[depthPhase1] && search->ax[depthPhase1 - 1] != search->ax[depthPhase1] + 3)) {

                    free((void*) fc);
                    free((void*) cc);
                    free((void*) c);
                    if (useSeparator) {
                        res = solutionToString(search, s, depthPhase1);
                    } else {
                        res = solutionToString(search, s, -1);
                    }
                    free((void*) search);
                    free((void*)res);
                    stopSearch = true;
                }
            }

        }
    } while (!stopSearch);





    return rotations;
}

void initarrays()
{
    cubiecube_t* a;
    cubiecube_t* moveCube = get_moveCube();

    int i;
    int k, j;
    a = get_cubiecube();

    for (i = 0; i < N_URFtoDLF; i++) {
        setURFtoDLF(a, (short)i);
        for (j = 0; j < 6; j++) {
            for (k = 0; k < 3; k++) {
                cornerMultiply(a, &moveCube[j]);
                URFtoDLF_Move[i][3 * j + k] = getURFtoDLF(a);
            }
            cornerMultiply(a, &moveCube[j]);
        }

    }
    free(a);


    a = get_cubiecube();
    for (i = 0; i < N_URtoDF; i++) {
        setURtoDF(a, i);

        for (j = 0; j < 6; j++) {
            for (k = 0; k < 3; k++) {
                edgeMultiply(a, &moveCube[j]);
                URtoDF_Move[i][3 * j + k] = (short) getURtoDF(a);
                // Table values are only valid for phase 2 moves!
                // For phase 1 moves, casting to short is not possible.
            }
            edgeMultiply(a, &moveCube[j]);
        }
    }


    free(a);
}