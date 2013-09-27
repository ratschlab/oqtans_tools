// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include "io_dp.h"
#include "common_dp.h"

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

FILE* CIO_DP::target=stdout;

CIO_DP::CIO_DP()
{
}

void CIO_DP::message(EMessageType_DP prio, const CHAR *fmt, ... )
{
	check_target();
	print_message_prio(prio, target);
    va_list list;
    va_start(list,fmt);
    vfprintf(target,fmt,list);
    va_end(list);
    fflush(target);
}

void CIO_DP::buffered_message(EMessageType_DP prio, const CHAR *fmt, ... )
{
	check_target();
	print_message_prio(prio, target);
    va_list list;
    va_start(list,fmt);
    vfprintf(target,fmt,list);
    va_end(list);
}

CHAR* CIO_DP::skip_spaces(CHAR* str)
{
	INT i=0;

	if (str)
	{
		for (i=0; isspace(str[i]); i++);

		return &str[i];
	}
	else 
		return str;
}

void CIO_DP::set_target(FILE* t)
{
	target=t;
}

void CIO_DP::check_target()
{
	if (!target)
		target=stdout;
}

void CIO_DP::print_message_prio(EMessageType_DP prio, FILE* target_)
{
	switch (prio)
	{
		case M_DEBUG_DP:
			fprintf(target_, "[DEBUG] ");
			break;
		case M_PROGRESS_DP:
			fprintf(target_, "[PROGRESS] ");
			break;
		case M_INFO_DP:
			//fprintf(target_, "[INFO] ");
			break;
		case M_NOTICE_DP:
			fprintf(target_, "[NOTICE] ");
			break;
		case M_WARN_DP:
			fprintf(target_, "[WARN] ");
			break;
		case M_ERROR_DP:
			fprintf(target_, "[ERROR] ");
			break;
		case M_CRITICAL_DP:
			fprintf(target_, "[CRITICAL] ");
			break;
		case M_ALERT_DP:
			fprintf(target_, "[ALERT] ");
			break;
		case M_EMERGENCY_DP:
			fprintf(target_, "[EMERGENCY] ");
			break;
		case M_MESSAGEONLY_DP:
			break;
		default:
			break;
	}
}
