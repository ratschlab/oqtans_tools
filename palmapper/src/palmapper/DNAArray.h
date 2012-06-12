#ifndef __DNAARRAY_H__
#define __DNAARRAY_H__

#include <assert.h>

using namespace shogun ;

//#define DNAARRAY_PERFORM_SANITY_CHECKS
#include "shogun/Array.h"

//using namespace shogun ;

class CDNAArray: public CArray<unsigned char>
{
protected:
	unsigned char dna2num(char dna) const
		{
			switch (dna)
			{
			case 'A':
			case 'a': return 0 ;
			case 'c':
			case 'C': return 1 ;
			case 'g':
			case 'G': return 2 ;
			case 't':
			case 'T': return 3 ;
			case 'n':
			case 'N': return 4 ;
			case 'r':
			case 'R': return 5 ;
			case 'y':
			case 'Y': return 6 ;
			case 'm':
			case 'M': return 7 ;
			case 'k':
			case 'K': return 8 ;
			case 'w':
			case 'W': return 9 ;
			case 's':
			case 'S': return 10 ;
			case 'b':
			case 'B': return 11 ;
			case 'd':
			case 'D': return 12 ;
			case 'h':
			case 'H': return 13 ;
			case 'v':
			case 'V': return 14 ;
			} ;
			return 15 ;
		} ;

	char num2dna(unsigned char num) const 
		{
			switch (num)
			{
			case 0: return 'A' ;
			case 1: return 'C' ;
			case 2: return 'G' ;
			case 3: return 'T' ;
			case 4: return 'N' ;
			case 5: return 'R' ;
			case 6: return 'Y' ;
			case 7: return 'M' ;
			case 8: return 'K' ;
			case 9: return 'W' ;
			case 10: return 'S' ;
			case 11: return 'B' ;
			case 12: return 'D' ;
			case 13: return 'H' ;
			case 14: return 'V' ;
			} ;
			return 'N' ;
		} ;
	char m_num2dna[16] ;
	unsigned char m_dna2num[256] ;
	
public:
	CDNAArray(char* seq, size_t len): CArray<unsigned char>((len+1)/2)
		{
			for (int i=0; i<16; i++)
				m_num2dna[i]=num2dna(i) ;
			for (int i=0; i<256; i++)
				m_dna2num[i]=dna2num(i) ;

			for (size_t i=0; i<(len+1)/2; i++)
			{
				unsigned char n = m_dna2num[(size_t)seq[i*2]] ;
				//unsigned char n1=n ;
				
				n = n*16 ;
				if (i*2+1<len)
					n += m_dna2num[(size_t)seq[i*2+1]] ;
				set_element(n, i) ;
				//if (i<10)
				//	fprintf(stdout, "set elem %i: %i (%i) %c %c\n", (int)i, (int)n, (int)n1, seq[i*2], seq[i*2+1]) ;
			}
		}

	CDNAArray(size_t len): CArray<unsigned char>((len+1)/2)
		{
			for (int i=0; i<16; i++)
				m_num2dna[i]=num2dna(i) ;

			for (int i=0; i<256; i++)
				m_dna2num[i]=dna2num(i) ;

			for (size_t i=0; i<(len+1)/2; i++)
				set_element(0, i) ;
		}

	static const char*get_class_name()
		{
			return "DNAArray" ;
		}

	char get_char(size_t index) const
		{
			size_t i=index/2 ;
			unsigned char n = get_element(i) ;
			//fprintf(stdout, "elem %i: %i\n", (int)index, (int) n) ;
			
			if (index%2==0)
				n=n>>4 ;
			else
				n=n&15 ;
			return m_num2dna[n] ;
		}

	unsigned char get_elem(size_t index) const
		{
			size_t i=index/2 ;
			unsigned char n = get_element(i) ;
			
			if (index%2==0)
				n=n>>4 ;
			else
				n=n&15 ;

			return n ;
		}
	
	void set_elem(size_t index, unsigned char elem)
		{
#ifdef DNAARRAY_PERFORM_SANITY_CHECKS
			assert(elem<16) ;
#endif
			size_t i=index/2 ;
			unsigned char n = get_element(i) ;
			//unsigned char n1 = n ;

			if (index%2==0)
				n = (n & 15) + elem*16 ;
			else
				n = (n & 240) + elem ;

			set_element(n, i) ;

#ifdef DNAARRAY_PERFORM_SANITY_CHECKS
			assert(get_elem(index)==elem) ;
#endif
		}
	
} ;

#endif
