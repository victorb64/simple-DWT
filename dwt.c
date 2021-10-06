#include "dwt.h"

void DWT_Init(DWT_ctx* ctx, uint32_t filtersLength, float32_t* pHighFilter, float32_t* pLowFilter, DWT_Padding padding, uint32_t decLevel, float32_t* workBuffer, uint32_t workBufferSize, uint32_t dataSize)
{
	ctx->filterSize = filtersLength;
	ctx->highPassFilter = pHighFilter;
	ctx->lowPassFilter = pLowFilter;
	ctx->padding = padding;
	ctx->decompLevel = decLevel;
	uint32_t datasize = DWT_GetOutputSize(dataSize, decLevel, filtersLength);
	uint32_t extendedbuffersize = ((dataSize + 2 * filtersLength) - 2);
	ctx->workBuffer = workBuffer;
	ctx->workBufferSize = workBufferSize;
	ctx->flippedLowPassFilter = workBuffer + datasize;
	ctx->flippedHighPassFilter = workBuffer + datasize + filtersLength;
	ctx->decompBuffer = workBuffer + datasize + (filtersLength * 2);
	ctx->recompBuffer = workBuffer + datasize + (filtersLength * 2) + extendedbuffersize;
	ctx->multABuffer = workBuffer + datasize + (filtersLength * 2) + extendedbuffersize + extendedbuffersize;
	ctx->multDBuffer = workBuffer + datasize + (filtersLength * 2) + extendedbuffersize + extendedbuffersize + filtersLength;

	for (uint32_t i = 0; i < ctx->filterSize; i++)
	{
		ctx->flippedLowPassFilter[i] = ctx->lowPassFilter[ctx->filterSize - i - 1];
		ctx->flippedHighPassFilter[i] = ctx->highPassFilter[ctx->filterSize - i - 1];
	}
}

void DWT_Recomposition(DWT_ctx* ctx, DWT_Out* recOut, float32_t* outSignal, uint32_t numSamples)
{
	DWT_Coefficients cf;
	// erm...ugh...
	if (recOut->decompLevel == 1)
	{
		cf.ACoeffs = recOut->coeffs[1].coeffs;
		cf.DCoeffs = recOut->coeffs[0].coeffs;
		cf.size = recOut->coeffs[1].size;
		DWT_SingleRec(ctx, &cf, outSignal, numSamples);
	}
	else if (recOut->decompLevel >= 2)
	{
		int currentlevel = recOut->decompLevel;
		cf.ACoeffs = recOut->coeffs[currentlevel].coeffs;
		cf.DCoeffs = recOut->coeffs[currentlevel - 1].coeffs;
		cf.size = recOut->coeffs[currentlevel].size;
		DWT_SingleRec(ctx, &cf, outSignal, recOut->coeffs[currentlevel - 2].size);
		currentlevel -= 2;
		while (currentlevel > 0)
		{
			cf.ACoeffs = outSignal;
			cf.DCoeffs = recOut->coeffs[currentlevel].coeffs;
			cf.size = recOut->coeffs[currentlevel].size;
			DWT_SingleRec(ctx, &cf, outSignal, recOut->coeffs[currentlevel - 1].size);
			currentlevel--;
		}
		if (currentlevel == 0)
		{
			cf.ACoeffs = outSignal;
			cf.DCoeffs = recOut->coeffs[0].coeffs;
			cf.size = recOut->coeffs[0].size;
			DWT_SingleRec(ctx, &cf, outSignal, numSamples);
		}
	}
	int a = 0;
	a++;
}

void DWT_Decomposition(
	DWT_ctx* ctx,
	float32_t* pSignal,
	uint32_t numSamples,
	DWT_Out* output)
{
	uint32_t i;
	uint32_t outOffset = 0;
	uint32_t decOutputSize = 1 + ((numSamples - 1) / 2) + (ctx->filterSize / 2 - 1);

	DWT_Coefficients decompResult;
	decompResult.ACoeffs = ctx->recompBuffer; // we reuse this bit of mem
	decompResult.DCoeffs = &(ctx->workBuffer[outOffset]);

	// first decomp...
	DWT_SingleDec(ctx, pSignal, numSamples, &decompResult);

	output->coeffs[0].coeffs = &(ctx->workBuffer[outOffset]);
	output->coeffs[0].size = decompResult.size;

	outOffset += decompResult.size;

	// decompose decompLevel times...
	for (i = 1; i < ctx->decompLevel; i++)
	{
		decompResult.DCoeffs = &(ctx->workBuffer[outOffset]);
		// next decompositions
		DWT_SingleDec(ctx, decompResult.ACoeffs, decompResult.size, &decompResult);

		output->coeffs[i].coeffs = &(ctx->workBuffer[outOffset]);
		output->coeffs[i].size = decompResult.size;

		// offset inc
		outOffset += decompResult.size;
	}

	// results
	vector_copy_f32(decompResult.ACoeffs, &(ctx->workBuffer[outOffset]), decompResult.size);

	output->coeffs[ctx->decompLevel].coeffs = &(ctx->workBuffer[outOffset]);
	output->coeffs[ctx->decompLevel].size = decompResult.size;

	outOffset += decompResult.size;

	output->decompLevel = ctx->decompLevel;

}
void DWT_SingleRec(DWT_ctx* ctx, DWT_Coefficients* inputCoeffs, float32_t* outSignal, uint32_t blockSize)
{
	uint32_t i, j;
	float32_t* multResCA = ctx->multABuffer;
	float32_t* multResCD = ctx->multDBuffer;

	// upsample signal
	float32_t* paddedACoeff = ctx->decompBuffer; // we reuse this bit
	float32_t* paddedDCoeff = ctx->recompBuffer;

	memset(paddedACoeff, 0, sizeof(float32_t)*(inputCoeffs->size * 2) + ctx->filterSize);
	memset(paddedDCoeff, 0, sizeof(float32_t)*(inputCoeffs->size * 2) + ctx->filterSize);
	for (uint32_t i = 0; i < (inputCoeffs->size * 2); i++)
	{
		if (i % 2 == 0)
		{
			paddedACoeff[i] = 0;
			paddedDCoeff[i] = 0;
		}
		else
		{
			paddedACoeff[i] = inputCoeffs->ACoeffs[i / 2];
			paddedDCoeff[i] = inputCoeffs->DCoeffs[i / 2];
		}
	}

	for (i = 0; i < blockSize; i++)
	{
		vector_mult_f32(&(paddedACoeff[i]), ctx->lowPassFilter, multResCA, ctx->filterSize);
		vector_mult_f32(&(paddedDCoeff[i]), ctx->highPassFilter, multResCD, ctx->filterSize);
		outSignal[i] = 0;
		for (j = 0; j < ctx->filterSize; j++)
		{
			outSignal[i] += multResCA[j] + multResCD[j];
		}
	}
}

void DWT_SingleDec(DWT_ctx* ctx, float32_t* inputSignal, uint32_t blockSize, DWT_Coefficients* outputCoeffs)
{
	uint32_t i, j;
	// pad with filterlength-1 and decimate/2
	uint32_t  outSize = (blockSize + ctx->filterSize - 1) / 2;

	// work buffer needs filter-1 on left and right side
	uint32_t  workBufferSize = blockSize + 2 * ctx->filterSize - 2;

	switch (ctx->padding)
	{
	case DWT_PADDING_SYMMETRIC:
		for (i = 0; i < ctx->filterSize - 2; i++)
		{
			ctx->decompBuffer[i] = inputSignal[ctx->filterSize - 3 - i];
		}
		for (uint32_t i = 0; i < blockSize; i++)
		{
			ctx->decompBuffer[(ctx->filterSize - 2) + i] = inputSignal[i];
		}
		for (i = 0; i < ctx->filterSize; i++)
		{
			ctx->decompBuffer[(ctx->filterSize - 2) + blockSize + i] = inputSignal[blockSize - 1 - i];
		}
		break;
	default:
		break;
	}

	for (i = 0; i < outSize; i++)
	{
		vector_mult_f32(&(ctx->decompBuffer[2 * i]), ctx->flippedLowPassFilter, ctx->multABuffer, ctx->filterSize);
		vector_mult_f32(&(ctx->decompBuffer[2 * i]), ctx->flippedHighPassFilter, ctx->multDBuffer, ctx->filterSize);

		outputCoeffs->ACoeffs[i] = 0;
		outputCoeffs->DCoeffs[i] = 0.0;

		for (j = 0; j < ctx->filterSize; j++)
		{
			outputCoeffs->ACoeffs[i] += ctx->multABuffer[j];
			outputCoeffs->DCoeffs[i] += ctx->multDBuffer[j];
		}
	}

	outputCoeffs->size = outSize;
}

uint32_t DWT_GetOutputSize(uint32_t inputSize, uint32_t decLevel, uint32_t filterLength)
{
	uint32_t res = 0;

	uint32_t decOutputSize = inputSize;
	for (uint32_t i = 0; i < decLevel; i++)
	{
		decOutputSize = (decOutputSize + (filterLength - 1)) / 2;
		res += decOutputSize;
	}
	res += decOutputSize;
	return res;
}
uint32_t DWT_GetRequiredBufferSize(uint32_t inputSize, uint32_t decLevel, uint32_t filterLength)
{
	uint32_t outsize = DWT_GetOutputSize(inputSize, decLevel, filterLength);
	outsize += 4 * filterLength;
	outsize += 2 * ((inputSize + 2 * filterLength) - 2);
	return outsize;
}
void vector_mult_f32(float* A, float* B, float* res, uint32_t size)
{
	while (size > 0U)
	{
		*res++ = (*A++) * (*B++);
		size--;
	}
}

void vector_copy_f32(float* src, float* dest, uint32_t size)
{
	while (size > 0U)
	{
		*dest++ = *src++;
		size--;
	}
}
