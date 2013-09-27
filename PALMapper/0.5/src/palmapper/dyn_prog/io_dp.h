// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#ifndef __CIO_DP_H__
#define __CIO_DP_H__

#include "common_dp.h"

#include <stdio.h>
#include <stdarg.h>

class CIO_DP
{
public:
	CIO_DP();

	static void set_target(FILE* target);
	static void message(EMessageType_DP prio, const CHAR *fmt, ... );

	inline static void not_implemented() 
	{
		message(M_ERROR_DP, "Sorry, not yet implemented\n");
	}

	static void buffered_message(EMessageType_DP prio, const CHAR *fmt, ... );

	static CHAR* skip_spaces(CHAR* str);

protected:
	static void check_target();
	static void print_message_prio(EMessageType_DP prio, FILE* target);

protected:
	static FILE* target;
};
#endif
