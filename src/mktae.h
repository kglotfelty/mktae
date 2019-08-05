/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */

#include <dslib.h>
#include <pixlib.h>
#include <ardlib.h>
#include <stack.h>

/* 

 gcc -O3 -Wall -g tdet.c -I/vobs/ASC_BUILD/include \
  -I/vobs/ASC_BUILD/src/include -L/vobs/ASC_BUILD/lib \
  -lardlib -lds -lpix -lascdm

*/
typedef struct {
  long nasp;
  double *ra;
  double *dec;
  double *roll;
  double *duration;
} Asp_Data;


dmBlock *load_input_setup_output( char *infile,
				  char *outfile,
				  long *laxes,
				  dmDescriptor **phy_coord,
				  regRegion **reg
				  );


int load_asp_data( char * ahstfile, Asp_Data *Asp_Hist, double *sim );

int setup_pixlib_and_ardlib( char *infile, char *geomfile, 
			     char *ardlib, double *sim,
			     Ardlib_Det_BadPix_Type **bpix );




float*  process_pixels( long *laxes, 
			dmDescriptor *phy_coord,
			Asp_Data *Asp_Hist,
			Ardlib_Det_BadPix_Type **bpix,
			regRegion *reg
			);

void output_data( dmBlock *outBlock, float *outdata, long *laxes );

void cleanup( float *outdata, 
	      Asp_Data *Asp_Hist, 
	      regRegion *reg, 
	      Ardlib_Det_BadPix_Type **bpix);

int load_params( char *infile, 
		char *ahstfile,
		char *outfile,
		char *geomfile,
		char *ardlib,
		short *clobber,
		short *verbose
		 ) ;

int mktae(void);
