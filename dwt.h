#ifndef DWT_H_
#define DWT_H_
#include <stdint.h>
#include <string.h>

typedef float float32_t;

#define MAX_DECOMP_LEVELS 10
#define POW2(x) (1 << x)

// figure out needed work buffer sizes at compile time...
// we also store the flipped filters in the work buffer
// there should be a nicer way but this will do
#define DECSIZE_LVL(x,f,n) ((x+(f-1)*(POW2(n)-1))/POW2(n))
#define DWTBUFSIZE_1(x,f) (4*f + (2*((x+2*f)-2))+(2*(DECSIZE_LVL(x,f,2))))
#define DWTBUFSIZE_2(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + 2*(DECSIZE_LVL(x,f,2))))
#define DWTBUFSIZE_3(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + 2*(DECSIZE_LVL(x,f,3))))
#define DWTBUFSIZE_4(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + DECSIZE_LVL(x,f,3) + 2*(DECSIZE_LVL(x,f,4))))
#define DWTBUFSIZE_5(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + DECSIZE_LVL(x,f,3) + DECSIZE_LVL(x,f,4) + 2*(DECSIZE_LVL(x,f,5))))
#define DWTBUFSIZE_6(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + DECSIZE_LVL(x,f,3) + DECSIZE_LVL(x,f,4) + DECSIZE_LVL(x,f,5) + 2*(DECSIZE_LVL(x,f,6))))
#define DWTBUFSIZE_7(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + DECSIZE_LVL(x,f,3) + DECSIZE_LVL(x,f,4) + DECSIZE_LVL(x,f,5) + DECSIZE_LVL(x,f,6) + 2*(DECSIZE_LVL(x,f,7))))
#define DWTBUFSIZE_8(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + DECSIZE_LVL(x,f,3) + DECSIZE_LVL(x,f,4) + DECSIZE_LVL(x,f,5) + DECSIZE_LVL(x,f,6) + DECSIZE_LVL(x,f,7) + 2*(DECSIZE_LVL(x,f,8))))
#define DWTBUFSIZE_9(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + DECSIZE_LVL(x,f,3) + DECSIZE_LVL(x,f,4) + DECSIZE_LVL(x,f,5) + DECSIZE_LVL(x,f,6) + DECSIZE_LVL(x,f,7) + DECSIZE_LVL(x,f,8) + 2*(DECSIZE_LVL(x,f,9))))
#define DWTBUFSIZE_10(x,f) (4*f + (2*((x+2*f)-2))+(DECSIZE_LVL(x,f,1) + DECSIZE_LVL(x,f,2) + DECSIZE_LVL(x,f,3) + DECSIZE_LVL(x,f,4) + DECSIZE_LVL(x,f,5) + DECSIZE_LVL(x,f,6) + DECSIZE_LVL(x,f,7) + DECSIZE_LVL(x,f,8) + DECSIZE_LVL(x,f,9) + 2*(DECSIZE_LVL(x,f,10))))


// these will implement vector multiplication and copy according to current platform
void vector_mult_f32(float* A, float* B, float* res, uint32_t size);
void vector_copy_f32(float* src, float* dest, uint32_t size);

// only care about symmetric padding at this point
typedef enum
{
	DWT_PADDING_SYMMETRIC,
} DWT_Padding;

// we use a single block buffer for all the process data, but these pointers come in handy
typedef struct
{
	DWT_Padding padding;
	uint32_t decompLevel;
	uint32_t filterSize;
	uint32_t workBufferSize;
	float32_t* highPassFilter;
	float32_t* lowPassFilter;
	float32_t* flippedHighPassFilter;
	float32_t* flippedLowPassFilter;
	float32_t* decompBuffer;
	float32_t* recompBuffer;
	float32_t* multABuffer;
	float32_t* multDBuffer;
	float32_t* workBuffer;
} DWT_ctx;


typedef struct
{
	uint32_t size;
	float32_t* DCoeffs;
	float32_t* ACoeffs;
} DWT_Coefficients;

typedef struct
{
	float32_t* coeffs;
	uint32_t size;
} DWT_OutCoeffs;

typedef struct
{
	DWT_OutCoeffs coeffs[MAX_DECOMP_LEVELS + 1 ]; // +1 because last level has both cA and cD
	uint32_t decompLevel;
} DWT_Out;


void DWT_Init(DWT_ctx* ctx,	uint32_t filtersLength, float32_t* pHighFilter, float32_t* pLowFilter, DWT_Padding padding, uint32_t decLevel, float32_t* pInternalBuffer, uint32_t internalBufferSize, uint32_t numSamples);

void DWT_Decomposition(DWT_ctx* ctx, float32_t* pSignal, uint32_t numSamples, DWT_Out* pOut);

void DWT_Recomposition(DWT_ctx* ctx,DWT_Out* recOut, float32_t* outSignal,uint32_t numSamples);
void DWT_SingleDec(DWT_ctx* ctx, float32_t* inputSignal, uint32_t blockSize, DWT_Coefficients* outputCoeffs);
void DWT_SingleRec(DWT_ctx* ctx, DWT_Coefficients* inputCoeffs, float32_t* outSignal, uint32_t blockSize);

uint32_t DWT_GetOutputSize(uint32_t inputSize, uint32_t decLevel, uint32_t filterLength);
uint32_t DWT_GetRequiredBufferSize(uint32_t inputSize, uint32_t decLevel, uint32_t filterLength);
#endif 
