#pragma once
// Authors: Korbinian Schneeberger, Stephan Ossowski, Joerg Hagmann, Gunnar Raetsch
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "shogun/init.h"
#include "shogun/common.h"
#include "shogun/io.h"
#include "shogun/Array.h"
#include "DNAArray.h"
#include "DNAArray4.h"

//#define BinaryStream_MAP

#ifdef BinaryStream_MAP
#include "BinaryStream.h"
#endif


// ##############################################################
// ####### GLOBAL VARIABLES #####################################
// ##############################################################

#include <palmapper/Config.h>
#include <palmapper/Statistics.h>
#include <palmapper/QPalma.h>

extern Config _config;
extern Statistics _stats;
