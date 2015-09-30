// Authors: Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#ifndef __COMMON_DP_H__
#define __COMMON_DP_H__

#include <stdlib.h> 
#include <stdio.h> 
#include "config_dp.h"

#ifdef SUNOS
#define bool int
#define false 0
#define true 1
#endif

/*#ifndef LINUX
#define RANDOM_MAX 2147483647
#else
#define RANDOM_MAX RAND_MAX
#endif*/

/**@name Standard Types 
 * Definition of Platform independent Types
*/
//@{

/// Type CHAR
typedef char CHAR;
typedef CHAR* P_CHAR;

/// Type BYTE 
typedef unsigned char BYTE;
typedef BYTE* P_BYTE;

/// Type SHORT is 2 bytes in size
typedef short int SHORT;
typedef SHORT* P_SHORT;

/// Type WORD is 2 bytes in size
typedef unsigned short int WORD;
typedef WORD* P_WORD;

/// Type INT is 4 bytes in size
typedef int INT;
typedef INT* P_INT;

/// Type INT is 4 bytes in size
typedef unsigned int UINT;
typedef UINT* P_UINT;

/// Type LONG is 8 bytes in size
typedef long LONG;
typedef LONG* P_LONG;

/// Type SHORTREAL is 4 bytes in size
typedef float SHORTREAL;
typedef SHORTREAL* P_SHORTREAL;

/// Type REAL is 8 bytes in size
typedef double REAL;
typedef REAL* P_REAL;

/// Type LONGREAL is 16 bytes in size
//typedef long double LONGREAL;
//typedef LONGREAL* P_LONGREAL;

#ifdef USE_SHORTREAL_KERNELCACHE
	typedef SHORTREAL KERNELCACHE_ELEM;
#else
	typedef REAL KERNELCACHE_ELEM;
#endif

typedef KERNELCACHE_ELEM P_KERNELCACHE_ELEM;

typedef LONG KERNELCACHE_IDX;

/// The io libs output [DEBUG] etc in front of every CIO::message
/// 'higher' messages filter output depending on the loglevel, i.e. CRITICAL messages
/// will print all M_CRITICAL TO M_EMERGENCY messages to
enum EMessageType_DP
{
	M_DEBUG_DP,
	M_INFO_DP,
	M_NOTICE_DP,
	M_WARN_DP,
	M_ERROR_DP,
	M_CRITICAL_DP,
	M_ALERT_DP,
	M_EMERGENCY_DP,
	M_MESSAGEONLY_DP,
	M_PROGRESS_DP
};

/// Alphabet of charfeatures/observations
enum E_ALPHABET_DP
{
	/// DNA - letters A,C,G,T,*,N,n
	DNA_DP=0,

	/// PROTEIN - letters a-z
	PROTEIN_DP=1,

	/// ALPHANUM - [0-9a-z]
	ALPHANUM_DP=2,

	/// CUBE - [1-6]
	CUBE_DP=3,

	/// NONE - type has no alphabet
	NONE_DP=4
};

//@}

#define TMP_DIR "/tmp/"
//#define TMP_DIR "/short/x46/tmp/"

#endif
