#ifndef GENERAL_DEFINITIONS_DEFINED
#define GENERAL_DEFINITIONS_DEFINED

#define LST_4_1     {  1,   2,   3,   4}
#define LST_4_3     {501, 502, 503, 504}

#define LST_5_1     {  1,   2,   3,   4,  5}

using cx_float  = std::complex<float>;
using cx_double = std::complex<double>;

using vi_init_list = std::initializer_list<int>;
using vf_init_list = std::initializer_list<float>;
using vd_init_list = std::initializer_list<double>;


static const vi_init_list   il_4_1 = LST_4_1;
static const vf_init_list   fl_4_1 = LST_4_1;
static const vd_init_list   db_4_1 = LST_4_1;

static const vi_init_list   il_4_3 = LST_4_3;
static const vf_init_list   fl_4_3 = LST_4_3;
static const vd_init_list   db_4_3 = LST_4_3;

static const vi_init_list   il_5_1 = LST_5_1;
static const vf_init_list   fl_5_1 = LST_5_1;
static const vd_init_list   db_5_1 = LST_5_1;



#define LST_14_2       {{11, 12, 13, 14}}
#define LST_15_1       {{1, 2, 3, 4, 5}}
#define LST_15_CPX       {{{1,2},{2,0},{3,-1},{4,0},{5,-8}}}

#define LST_31_1       {{1}, {2}, {3}}
#define LST_41_1       {{1}, {2}, {3}, {4}}
#define LST_51_1       {{1}, {2}, {3}, {4}, {5}}
#define LST_51_CPX     {{{1,2}}, {{2,0}}, {{3,-1}}, {{4,8}}, {{5,-3}}}

#define ZEROS          {{0, 0, 0},          \
                        {0, 0, 0},          \
                        {0, 0, 0}}

#define IDENTITY       {{1, 0, 0},          \
                        {0, 1, 0},          \
                        {0, 0, 1}}

#define LST_33_M       {{1, 2, 4, 0, 0, 0},          \
                        {6, 8, 10, 0, 0, 0},         \
                        {12, 14, 16, 0, 0, 0}}

#define LST_33_2       {{ 10, 20, 30 },     \
                        { 40, 50, 60 },     \
                        { 70, 80, 90 }}

#define LST_33_2_T     {{ 10, 40, 70 },     \
                        { 20, 50, 80 },     \
                        { 30, 60, 90 }}

#define LST_33_1       {{ 1, -2, 3 },       \
                        { 2, -3, -4 },      \
                        { -3, 4, 5 }}

#define LST_33_CPX     {{{1, 0}, {-2, 0}, {3, 1}}, \
                        {{2, 3}, {2, 5}, {5, 8}}, \
                        {{-3, -9}, {4, -1}, {5, 0}}}

#define LST_34_2       {{ 11, 12, 13, 14 }, \
                        { 21, 22, 23, 24 }, \
                        { 31, 32, 33, 34 }}

#define LST_34_2_N     {{ -11, -12, -13, -14 }, \
                        { -21, -22, -23, -24 }, \
                        { -31, -32, -33, -34 }}

#define LST_34_2_T     {{ 11, 21, 31 },         \
                        { 12, 22, 32 },         \
                        { 13, 23, 33 },         \
                        { 14, 24, 34 }}

#define LST_34_2_T_N   {{ -11, -21, -31 },      \
                        { -12, -22, -32 },      \
                        { -13, -23, -33 },      \
                        { -14, -24, -34 }}

#define LST_43_1       {{  1,  2,  3 },         \
                        {  5,  6,  7 },         \
                        {  1,  3,  2 },         \
                        {  0,  2,  1 }}

#define LST_44_2       {{ 11, 12, 13, 14 },     \
                        { 21, 22, 23, 24 },     \
                        { 31, 32, 33, 34 },     \
                        { 41, 42, 43, 44 }}

#define LST_55_2       {{ 11, 12, 13, 14, 15 }, \
                        { 21, 22, 23, 24, 25 }, \
                        { 31, 32, 33, 34, 35 }, \
                        { 41, 42, 43, 44, 45 }, \
                        { 51, 52, 53, 54, 55 }}

using mi_init_list = std::initializer_list<std::initializer_list<int>>;
using mf_init_list = std::initializer_list<std::initializer_list<float>>;
using md_init_list = std::initializer_list<std::initializer_list<double>>;
using mcd_init_list = std::initializer_list<std::initializer_list<cx_double>>;


static const mi_init_list   il_14_2 = LST_14_2;
static const mf_init_list   fl_14_2 = LST_14_2;
static const md_init_list   db_14_2 = LST_14_2;
static const mi_init_list   il_15_1 = LST_15_1;
static const mf_init_list   fl_15_1 = LST_15_1;
static const md_init_list   db_15_1 = LST_15_1;
static const mcd_init_list  cd_15_1 = LST_15_CPX;

static const mi_init_list   il_41_1 = LST_41_1;
static const mf_init_list   fl_41_1 = LST_41_1;
static const md_init_list   db_41_1 = LST_41_1;
static const md_init_list   db_31_1 = LST_31_1;

static const mi_init_list   il_51_1 = LST_51_1;
static const mf_init_list   fl_51_1 = LST_51_1;
static const md_init_list   db_51_1 = LST_51_1;
static const mcd_init_list  cd_51_2 = LST_51_CPX;

static const mi_init_list   il_33_1   = LST_33_1;
static const mf_init_list   fl_33_1   = LST_33_1;
static const md_init_list   db_33_1   = LST_33_1;
static const mcd_init_list  cd_33_1   = LST_33_CPX;
static const md_init_list   db_33_m   = LST_33_M;
static const md_init_list   db_zero   = ZEROS;
static const mi_init_list   il_ID     = IDENTITY;

static const mi_init_list   il_33_2 = LST_33_2;
static const mf_init_list   fl_33_2 = LST_33_2;
static const md_init_list   db_33_2 = LST_33_2;
static const mi_init_list   il_33_2_t = LST_33_2_T;
static const mf_init_list   fl_33_2_t = LST_33_2_T;
static const md_init_list   db_33_2_t = LST_33_2_T;

static const mf_init_list   fl_34_2     = LST_34_2;
static const mf_init_list   fl_34_2_n   = LST_34_2_N;
static const mf_init_list   fl_34_2_t   = LST_34_2_T;
static const mf_init_list   fl_34_2_t_n = LST_34_2_T_N;

static const mi_init_list   il_43_1 = LST_43_1;
static const mf_init_list   fl_43_1 = LST_43_1;
static const md_init_list   db_43_1 = LST_43_1;

static const mi_init_list   il_44_2 = LST_44_2;
static const mf_init_list   fl_44_2 = LST_44_2;
static const md_init_list   db_44_2 = LST_44_2;

static const mi_init_list   il_55_2 = LST_55_2;
static const mf_init_list   fl_55_2 = LST_55_2;
static const md_init_list   db_55_2 = LST_55_2;


#endif