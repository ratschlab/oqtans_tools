/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 2009 Soeren Sonnenburg
 * Copyright (C) 2009 Fraunhofer Institute FIRST and Max-Planck-Society
 *
 * Taken from shogun
 */

#ifndef __BINARYSTREAM_H__
#define __BINARYSTREAM_H__

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/** @brief memory mapped emulation via binary streams (files)
*
* Implements memory mapped file emulation (\sa CMemoryMappedFile) via standard
* file operations like fseek, fread etc
*/
template <class T> class CBinaryStream
{
	public:
		/** constructor
		 *
		 * open a file for read mode
		 *
		 * @param fname name of file, zero terminated string
		 * @param flag determines read or read write mode (currently only 'r' is supported)
		 */
		CBinaryStream(const char* fname, const char* flag="r")
		{
			rw=flag;
			m_fname=strdup(fname);
			fd = fopen(fname, flag);
			if (!fd)
			{
				perror("Error opening file\n");
				exit(1);
			}

			struct stat sb;
			if (stat(fname, &sb) == -1)
			{
				perror("Error determining file size\n");
				exit(1);
			}

			length = sb.st_size;
			//fprintf(stderr, "File %s is %ld bytes in size\n", fname, length);
		}

		CBinaryStream(const CBinaryStream &bs)
		{
			rw=bs.rw;
			m_fname=strdup(bs.m_fname);

			fd = fopen(m_fname, rw);
			if (!fd)
			{
				perror("Error (re-)opening file\n");
				exit(1);
			}

			struct stat sb;
			if (stat(m_fname, &sb) == -1)
			{
				perror("Error determining file size\n");
				exit(1);
			}

			length = sb.st_size;
			assert(length==bs.length);
		}

		/** destructor */
		virtual ~CBinaryStream()
		{
			free(m_fname);
			fclose(fd);
		}

		/** get the number of objects of type T cointained in the file 
		 *
		 * @return length of file
		 */
		uint64_t get_length() const
		{
			return length/sizeof(T);
		}

		/** get the size of the file in bytes
		 *
		 * @return size of file in bytes
		 */
		uint64_t get_size() const
		{
			return length;
		}

		void read_and_forget() const
			{
				const int block_size = 32768 ;
				T buf[block_size] ;
				for (int i=0; i < get_size(); i+=32768)
				{
					fread(buf, sizeof(T), block_size, fd) ;
				}
			}
		
		void pre_buffer(T* buffer, long index, long size) const
		{
			assert(index>=0);
			assert(size>=0);

			if (size>0)
			{
				if (fseek(fd, ((long) sizeof(T))*((long) index), SEEK_SET) != 0)
				{
					fprintf(stderr, "Error seeking to %ld\n", sizeof(T)*((int64_t) index));
					exit(1);
				}

				if ( fread(buffer, sizeof(T), size, fd) != size)
				{
					fprintf(stderr, "Error calling fread index=%ld size=%ld file_length%ld\n", index, size, get_length());
					exit(1);
				}
			}
		}

		/** read next
		 *
		 * @return next element
		 */
		inline T read_next() const
		{
			T ptr;
			if ( fread(&ptr, sizeof(T), 1, fd) != 1)
			{
				fprintf(stderr, "Error calling fread (file '%s')\n", m_fname);
				exit(1);
			}
			return ptr;
		}

		/** operator overload for file read only access
		 *
		 * DOES NOT DO ANY BOUNDS CHECKING
		 *
		 * @param index index
		 * @return element at index
		 */
		inline T operator[](long index) const
		{
			if (fseek(fd, ((long) sizeof(T))*((long) index), SEEK_SET) != 0)
			{
				fprintf(stderr, "Error seeking to %ld (file '%s')\n", sizeof(T)*((int64_t) index), m_fname);
				exit(1);
			}

			T ptr;
			if ( fread(&ptr, sizeof(T), 1, fd) != 1)
			{
				fprintf(stderr, "Error calling fread (file '%s')\n", m_fname);
				exit(1);
			}
			return ptr;
		}

	protected:
		/** file descriptor */
		FILE* fd;
		/** size of file */
		uint64_t length;
		/** mode */
		const char* rw;
		/* fname */
		char* m_fname;
};
#endif // BINARY_STREAM
