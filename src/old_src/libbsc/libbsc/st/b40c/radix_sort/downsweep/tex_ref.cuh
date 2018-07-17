/******************************************************************************
 * 
 * Copyright (c) 2010-2012, Duane Merrill.  All rights reserved.
 * Copyright (c) 2011-2012, NVIDIA CORPORATION.  All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. 
 * 
 ******************************************************************************/

/******************************************************************************
 * Texture references for downsweep kernels
 ******************************************************************************/

#pragma once

#include "../../util/error_utils.cuh"
#include "../../util/basic_utils.cuh"
#include "../../util/tex_vector.cuh"
#include "../../util/ns_umbrella.cuh"

B40C_NS_PREFIX
namespace b40c {
namespace radix_sort {
namespace downsweep {

// Anonymous namespace in host code for linking
#if (__CUB_CUDA_ARCH__ == 0)
namespace {
#endif

/******************************************************************************
 * Key textures
 ******************************************************************************/

/**
 * Templated texture reference for downsweep keys
 */
template <typename KeyVectorType>
struct TexKeys
{
	// Texture reference type
	typedef texture<KeyVectorType, cudaTextureType1D, cudaReadModeElementType> TexRef;

	static TexRef ref0;
	static TexRef ref1;

	/**
	 * Bind textures
	 */
	static cudaError_t BindTexture(
		void *d0,
		void *d1,
		size_t bytes)
	{
		cudaError_t retval = cudaSuccess;
		do {
			cudaChannelFormatDesc tex_desc = cudaCreateChannelDesc<KeyVectorType>();

			size_t offset;
			if (d0)
			{
				// Bind key texture ref0
				retval = cudaBindTexture(&offset, ref0, d0, tex_desc, bytes);
				if (retval = util::B40CPerror(retval, "cudaBindTexture failed", __FILE__, __LINE__)) break;

				// We need texture-segment aligned input
				if (offset)
				{
					retval = util::B40CPerror(cudaErrorTextureNotBound , "cudaBindTexture failed", __FILE__, __LINE__);
					break;
				}
			}
			if (d1)
			{
				// Bind key texture ref1
				retval = cudaBindTexture(0, ref1, d1, tex_desc, bytes);
				if (retval = util::B40CPerror(retval, "cudaBindTexture failed", __FILE__, __LINE__)) break;

				// We need texture-segment aligned input
				if (offset)
				{
					retval = util::B40CPerror(cudaErrorTextureNotBound , "cudaBindTexture failed", __FILE__, __LINE__);
					break;
				}
			}
		} while (0);

		return retval;
	}

	/**
	 * Unbind textures
	 */
	static cudaError_t UnbindTexture()
	{
		cudaError_t retval = cudaSuccess;
		do {
			retval = cudaUnbindTexture(ref0);
			if (util::B40CPerror(retval , "cudaUnbindTexture failed", __FILE__, __LINE__)) break;

			retval = cudaUnbindTexture(ref1);
			if (util::B40CPerror(retval , "cudaUnbindTexture failed", __FILE__, __LINE__)) break;
		} while (0);

		return retval;
	}
};

// Texture reference definitions
template <typename KeyVectorType>
typename TexKeys<KeyVectorType>::TexRef TexKeys<KeyVectorType>::ref0 = 0;

template <typename KeyVectorType>
typename TexKeys<KeyVectorType>::TexRef TexKeys<KeyVectorType>::ref1 = 0;



/******************************************************************************
 * Value textures
 ******************************************************************************/

/**
 * Templated texture reference for downsweep values
 */
template <typename ValueVectorType>
struct TexValues
{
	// Texture reference type
	typedef texture<ValueVectorType, cudaTextureType1D, cudaReadModeElementType> TexRef;

	static TexRef ref0;
	static TexRef ref1;

	/**
	 * Bind textures
	 */
	static cudaError_t BindTexture(
		void *d0,
		void *d1,
		size_t bytes)
	{
		cudaError_t retval = cudaSuccess;
		do {
			cudaChannelFormatDesc tex_desc = cudaCreateChannelDesc<ValueVectorType>();

			size_t offset;
			if (d0)
			{
				// Bind key texture ref0
				retval = cudaBindTexture(&offset, ref0, d0, tex_desc, bytes);
				if (retval = util::B40CPerror(retval, "cudaBindTexture failed", __FILE__, __LINE__)) break;

				// We need texture-segment aligned input
				if (offset)
				{
					retval = util::B40CPerror(cudaErrorTextureNotBound , "cudaBindTexture failed", __FILE__, __LINE__);
					break;
				}
			}
			if (d1)
			{
				// Bind key texture ref1
				retval = cudaBindTexture(&offset, ref1, d1, tex_desc, bytes);
				if (retval = util::B40CPerror(retval, "cudaBindTexture failed", __FILE__, __LINE__)) break;

				// We need texture-segment aligned input
				if (offset)
				{
					retval = util::B40CPerror(cudaErrorTextureNotBound , "cudaBindTexture failed", __FILE__, __LINE__);
					break;
				}
			}
		} while (0);

		return retval;
	}

	/**
	 * Unbind textures
	 */
	static cudaError_t UnbindTexture()
	{
		cudaError_t retval = cudaSuccess;
		do {
			retval = cudaUnbindTexture(ref0);
			if (util::B40CPerror(retval , "cudaUnbindTexture failed", __FILE__, __LINE__)) break;

			retval = cudaUnbindTexture(ref1);
			if (util::B40CPerror(retval , "cudaUnbindTexture failed", __FILE__, __LINE__)) break;
		} while (0);

		return retval;
	}
};

// Texture reference definitions
template <typename ValueVectorType>
typename TexValues<ValueVectorType>::TexRef TexValues<ValueVectorType>::ref0 = 0;

template <typename ValueVectorType>
typename TexValues<ValueVectorType>::TexRef TexValues<ValueVectorType>::ref1 = 0;



/******************************************************************************
 * Texture types for downsweep kernel
 ******************************************************************************/

template <typename KeyType, typename ValueType, int THREAD_ELEMENTS>
struct Textures
{
	// Elements per texture load
	enum {
		KEY_ELEMENTS_PER_TEX		= util::TexVector<KeyType, THREAD_ELEMENTS>::ELEMENTS_PER_TEX,
		VALUE_ELEMENTS_PER_TEX		= util::TexVector<ValueType, THREAD_ELEMENTS>::ELEMENTS_PER_TEX,

		// If values are also going through tex, make sure the tex vector size is the same for both the key and value texture vector types
		ELEMENTS_PER_TEX			= util::NumericTraits<ValueType>::BUILT_IN ?
											CUB_MIN(int(KEY_ELEMENTS_PER_TEX), int(VALUE_ELEMENTS_PER_TEX)) :
											KEY_ELEMENTS_PER_TEX
	};

	typedef typename util::TexVector<
		KeyType,
		ELEMENTS_PER_TEX>::VecType KeyTexType;

	// Texture binding for downsweep values
	typedef typename util::TexVector<
		ValueType,
		ELEMENTS_PER_TEX>::VecType ValueTexType;
};

#if (__CUB_CUDA_ARCH__ == 0)
} // namespace (anonymous)
#endif

} // namespace downsweep
} // namespace radix_sort
} // namespace b40c
B40C_NS_POSTFIX
