#ifndef __DNAARRAY4_H__
#define __DNAARRAY4_H__

#include <assert.h>

//#define DNAARRAY4_PERFORM_SANITY_CHECKS

class CDNAArray4: public CArray<unsigned char>
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
			} ;
			return 0 ; // convert any other letter to 'A'
		} ;

	char num2dna(unsigned char num) const 
		{
			switch (num)
			{
			case 0: return 'A' ;
			case 1: return 'C' ;
			case 2: return 'G' ;
			case 3: return 'T' ;
			} ;
			assert(0) ;
			return 'A' ; // this should not happen
		} ;
	char m_num2dna[4] ;
	unsigned char m_dna2num[256] ;
	
public:
	CDNAArray4(char* seq, size_t len): CArray<unsigned char>((len+3)/4)
		{
			for (int i=0; i<4; i++)
				m_num2dna[i]=num2dna(i) ;
			for (int i=0; i<256; i++)
				m_dna2num[i]=dna2num(i) ;

			for (size_t i=0; i<(len+3)/4; i++)
			{
				unsigned char n ;

				n = m_dna2num[(size_t)seq[i*4]]<<6 ;
				if (i*4+1<len)
					n += m_dna2num[(size_t)seq[i*4+1]]<<4 ;
				if (i*4+2<len)
					n += m_dna2num[(size_t)seq[i*4+2]]<<2 ;
				if (i*4+3<len)
					n += m_dna2num[(size_t)seq[i*4+3]] ;

				set_element(n, i) ;
			}
		}

	CDNAArray4(size_t len): CArray<unsigned char>((len+3)/4)
		{
			for (int i=0; i<4; i++)
				m_num2dna[i]=num2dna(i) ;

			for (int i=0; i<256; i++)
				m_dna2num[i]=dna2num(i) ;

			for (size_t i=0; i<(len+3)/4; i++)
				set_element(0, i) ;
		}

	static const char*get_class_name()
		{
			return "DNAArray4" ;
		}

	char get_char(size_t index) const
		{
			unsigned char n = get_elem(index) ;

			return m_num2dna[n] ;
		}

	unsigned char get_elem(size_t index) const
		{
			size_t i=index/4 ;
			unsigned char n = get_element(i) ;
			
			if (index%4==0)
				n = (n>>6) & 3;
			else if (index%4==1)
				n = (n>>4) & 3 ;
			else if (index%4==2)
				n = (n>>2) & 3 ;
			else
				n = n & 3 ;

#ifdef DNAARRAY_PERFORM_SANITY_CHECKS
			assert(n<4) ;
#endif
			return n ;
		}
	
	void set_elem(size_t index, unsigned char elem)
		{
#ifdef DNAARRAY_PERFORM_SANITY_CHECKS
			assert(elem<4) ;
#endif
			size_t i=index/4 ;
			unsigned char n = get_element(i) ;
			//unsigned char n1=n ;
			
			if (index%4==0)
				n = (n & (~(3<<6))) + (elem<<6) ;
			else if (index%4==1)
				n = (n & (~(3<<4))) + (elem<<4) ;
			else if (index%4==2)
				n = (n & (~(3<<2))) + (elem<<2) ;
			else 
				n = (n & (~3)) + elem ;

			set_element(n, i) ;

#ifdef DNAARRAY_PERFORM_SANITY_CHECKS
			assert(get_elem(index)==elem) ;
#endif
		}
} ;

#endif
