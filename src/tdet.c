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


#include "mktae.h"
#include <histlib.h>

/* load image size , don't need the data */

dmBlock *load_input_setup_output( char *infile,
				  char *outfile,
				  long *laxes,
				  dmDescriptor **phy_coord,
				  regRegion **reg
				  )
{
  
  dmBlock *inBlock;
  char inst[100];
  dmDescriptor *inImg;
  long naxes, *_laxes;
  char coord_name[100];

  dmBlock *outBlock;
  dmDescriptor *outds;
  dmDescriptor *sky_dss;

  if ( NULL == (inBlock = dmImageOpen( infile ))) {
    fprintf( stderr, "Image Open\n");
    return(NULL);
  }

  if ( NULL == dmKeyRead_c( inBlock, "INSTRUME", inst, 99 )) {
    err_msg("WARNING: COULDN'T FIND INSTRUME, ASSUMING ACIS");
  } else {
    if ( ds_strcmp_cis( inst, "ACIS" ) != 0 ) {
      err_msg("ERROR: MUST BE ACIS DATA");
      return(NULL);
    }
  }

  if ( NULL == ( inImg = dmImageGetDataDescriptor (inBlock ))) {
    fprintf( stderr, "Image Descn");
    return(NULL);
  }
  
  if ( 2 != (naxes = dmGetArrayDimensions( inImg, &_laxes ))) {
    fprintf( stderr, "Array Dimn\n");
    return(NULL);
  }

  laxes[0] = _laxes[0];
  laxes[1] = _laxes[1];
  dmFree( _laxes );

  if ( NULL == ( *phy_coord = dmArrayGetAxisGroup( inImg, 1 ))) {
    err_msg("ERROR: Cannot get axis group");
    return(NULL);
  }

  dmGetName( *phy_coord, coord_name, 99 );
  if ( ds_strcmp_cis( coord_name, "TDET" ) != 0 ) {
    err_msg("NOT A TDET IMAGE");
    return(NULL);
  }


  /* Create the output file */

  if ( NULL == ( outBlock = dmImageCreate( outfile, dmFLOAT, laxes, naxes ))) {
    err_msg("ERROR CREATING OUTPUT FIEL ");
    return(NULL);
  }

  outds = dmImageGetDataDescriptor( outBlock );
  dmBlockCopyWCS( inBlock, outBlock );
  dmBlockCopy( inBlock, outBlock, "HEADER");
  dmBlockCopy( inBlock, outBlock, "SUBSPACE" );
  put_param_hist_info( NULL, "mktae", NULL, 0 );

  dmBlockClose( inBlock );

  if ( NULL == ( sky_dss = dmSubspaceColOpen( outBlock, "sky" ))) {
    err_msg("WARNING NO SKY DSS FOUND, USING ALL PIXELS");
    *reg = NULL;
  } else {
    if ( NULL == (*reg = dmSubspaceColGetRegion( sky_dss ))) {
      err_msg("WARNING NO REGION FOUND ON SKY DSS, USING ALL PIXELS");
    }
  }
  
  /* Re assign to the output coord */
  *phy_coord = dmArrayGetAxisGroup( outds, 1 );
  
  return(outBlock);
  
}

/* load asp data */
int load_asp_data( char * ahstfile, 
		   Asp_Data *Asp_Hist,
		   double *sim )
{
  Stack aspStk;
  char *cur;
  
  
  if ( NULL == (aspStk = stk_build( ahstfile ))) {
    err_msg("ERROR opening aspect stack \n");
    return(-45);
  }
  
  if (( 0 == stk_count( aspStk  )) ||
      ( NULL == stk_read_num(aspStk, 1 ) ) ||
      ( 0 == strlen( stk_read_num(aspStk, 1)) )) {
    err_msg("ERROR: EMpty stack\n");
    return(-56);
  }
  
  stk_rewind(aspStk);
  while ( NULL != (cur = stk_read_next( aspStk ))) {
    dmBlock *aspBlock;
    dmDescriptor *tmpCol;
    double roll;
    double ltime;
    long ii;
    char detnam[100];
    char *ptr;
    short ccd;

    if ( NULL == ( aspBlock = dmTableOpen(cur))) {
      fprintf( stderr, "Open Table\n");
      return(-4);
    }
    
    if ( NULL == dmKeyRead_c( aspBlock, "DETNAM", detnam, 99 )) {
      fprintf( stderr, "ERROR: Needs detname");
      return(-5);
    }
    if ( 0 != ds_strncmp_cis( detnam, "ACIS-", 5 )) {
      fprintf( stderr, "MUST BE ACIS ahist\n" );
      return(-6);
    }
    ptr = strchr( detnam, '-' );
    ptr++;
    if ( strlen(ptr) > 1 ) {
      fprintf( stderr, "ASPHIST MUST BE FOR 1 chip only, detname "
	       "says otherwise");
      return(-7);
    }
    ccd = *ptr - '0';
    
    if (( ccd < 0) || (ccd>9)) {
      err_msg( "DETNAME must be ACIS-#, not ACIS-I or -S");
      return(-2341);
    }

    Asp_Hist[ccd].nasp = dmTableGetNoRows( aspBlock );
    if ( Asp_Hist[ccd].nasp < 1 ) {
      err_msg("ZERO AHIST ROWS");
      return(-3451);
    }
    
    Asp_Hist[ccd].ra=(double*)calloc(Asp_Hist[ccd].nasp, sizeof(double));
    Asp_Hist[ccd].dec=(double*)calloc(Asp_Hist[ccd].nasp, sizeof(double));
    Asp_Hist[ccd].roll=(double*)calloc(Asp_Hist[ccd].nasp, sizeof(double));
    Asp_Hist[ccd].duration=(double*)calloc(Asp_Hist[ccd].nasp, sizeof(double));
    
    if ( NULL == ( tmpCol = dmTableOpenColumn( aspBlock, "X_OFFSET"))) {
      err_msg( "XOFFSET missing");
      return(-9);
    }
    dmGetScalars_d( tmpCol, Asp_Hist[ccd].ra, 1, Asp_Hist[ccd].nasp );
    
    if ( NULL == ( tmpCol = dmTableOpenColumn( aspBlock, "Y_OFFSET"))) {
      err_msg( "YOFFSET missing");
      return(-9);
    }
    dmGetScalars_d( tmpCol, Asp_Hist[ccd].dec, 1, Asp_Hist[ccd].nasp );
    
    if ( NULL == ( tmpCol = dmTableOpenColumn( aspBlock, "ROLL_OFFSET"))) {
      err_msg( "NO ROLL");
      return(-9);
    }
    dmGetScalars_d( tmpCol, Asp_Hist[ccd].roll, 1, Asp_Hist[ccd].nasp );
    
    if ( NULL == ( tmpCol = dmTableOpenColumn( aspBlock, "DURATION" ))) {
      err_msg( "ERROR CANNOT OPEN DURATION\n");
      return(-10);
    }
    dmGetScalars_d( tmpCol, Asp_Hist[ccd].duration, 1, Asp_Hist[ccd].nasp );
    
    
    if ( NULL == dmKeyRead_d( aspBlock, "ROLL_NOM", &roll )) {
      err_msg("ERROR NO RA NOM NOR PNT\n");
      return(-11);
    }
    
    /* We want to normalize by livetime since
       duration has dtcor applied, as does livetime */
    /* Hopefully don't have to go and look for LIVTIMEx*/
    if ( NULL == dmKeyRead_d( aspBlock, "LIVETIME", &ltime )) {
      err_msg("LIVETIME\n");
      return(-224);
    }
    
    for (ii=0;ii<Asp_Hist[ccd].nasp;ii++) {
      /* Roll needed is absolute */
      Asp_Hist[ccd].roll[ii] = roll+Asp_Hist[ccd].roll[ii]; 
      //      Asp_Hist[ccd].duration[ii] /= ltime;
    }
    
    /* Could/can/should only do this once or maybe check
       that they are all the same? */
    if (( NULL == dmKeyRead_d( aspBlock, "SIM_X", sim+0 )) ||
	( NULL == dmKeyRead_d( aspBlock, "SIM_Y", sim+1 )) ||
	( NULL == dmKeyRead_d( aspBlock, "SIM_Z", sim+2 ))) {
      err_msg("ERROR: BAD SIM");
      return(-9877);
    }
    dmTableClose( aspBlock );

  } /* end while stk_read_next */

  return(0);

}

int setup_pixlib_and_ardlib( char *infile, char *geomfile, 
			     char *ardlib, double *sim,
			     Ardlib_Det_BadPix_Type **bpix )
{
  short bb;
  
  if ( PIX_GOOD != pix_init_pixlib("CHANDRA", geomfile) ) {
    return(-3564);
  }
  pix_set_detector("ACIS");
  pix_set_aimpoint_by_value( sim );
  
  double zro[3] = { 0.0, 0.0, 0.0 };
  double aln[3] = { 0.0, 0.0, 0.0 };
  double ang[3] = { 0.0, 0.0, 0.0 };
  pix_set_mirror( zro, aln, ang );
  /* +++ history */
  
  /* ----------------------------------- */


#ifdef SETUP_BADPIX
  if ( ardlib && strlen( ardlib ) && ( 0 != ds_strcmp_cis(ardlib,"none"))) {
    ardlib_set_parfile( ardlib );
    
    if ( 0 != ardlib_initialize( "CHANDRA", infile )) {
      return(-99);
    }
    
    for (bb=0;bb<10;bb++) {
      char dnam[20];
      sprintf(dnam, "ACIS""-%d", bb);
      bpix[bb] = ardlib_open_det_badpix( dnam );
    }
  } else {
    
    for (bb=0;bb<10;bb++) {
      bpix[bb] = NULL;
    }
  } 
#endif
#ifndef SETUP_BADPIX
    for (bb=0;bb<10;bb++) {
      bpix[bb] = NULL;
    }
#endif


  
  /* +++ history */

  return(0);
}	     


/*
  pixel -> tdet (physical) -> chip -> det -> sky as f(aspect)
*/
float*  process_pixels( long *laxes, 
			dmDescriptor *phy_coord,
			Asp_Data *Asp_Hist,
			Ardlib_Det_BadPix_Type **bpix,
			regRegion *reg
			)
{
  
  /* Now loop through output */
  
  float *outdata;
  long ly;
  
  outdata = ( float*)calloc( laxes[0]*laxes[1], sizeof(float)); 
  if ( ( NULL == outdata ) ||
       ( 0 == laxes[0] ) ||
       ( 0 == laxes[1] ) ) {
    err_msg("ERRROR 0 pixels or 0 memeroy\n");
    return(NULL);
  }

  for (ly=0;ly<laxes[1];ly++ ) { /* for each row */
    long index; /* output array index */
    long lx;
    index = ly * laxes[0];

    for (lx=0;lx<laxes[0];lx++,index++) { /* for each column incr both*/
      long aa;
      double pix[2];
      double tdet[2];
      double chip[2];
      short  ccd;
      double det[2];

      pix[0] = lx+1;
      pix[1] = ly+1;

      dmCoordCalc_d( phy_coord, pix, tdet );
      if ( PIX_GOOD != pix_tdet_to_chip( tdet, &ccd, chip )) {
	continue;
      }

      if ( 0 == Asp_Hist[ccd].nasp ) {
	continue;
      }
      
      if (bpix[ccd] && ardlib_det_is_badpix( bpix[ccd], chip[0], chip[1] ) ) {
	continue; 	
      }

      pix_chip_to_fpc( ccd, chip, det );
      for (aa=0;aa<Asp_Hist[ccd].nasp;aa++) { /* for each asphist bin */
	double sky[2];
	double asp[3];
	asp[0] = Asp_Hist[ccd].ra[aa];
	asp[1] = Asp_Hist[ccd].dec[aa];
	asp[2] = Asp_Hist[ccd].roll[aa];

	pix_apply_aspect( det, asp, sky );
	
	if (( NULL==reg) || ( reg && regInsideRegion( reg, sky[0], sky[1] ))) {
	  outdata[index] += Asp_Hist[ccd].duration[aa];
	}
	
      } /* end for aa*/
    } /* end for lx */
  } /* end for ly */

  return( outdata );

}


void output_data( dmBlock *outBlock, float *outdata, long *laxes )
{
  
  dmSetArray_f( dmImageGetDataDescriptor(outBlock),
		outdata, laxes[0]*laxes[1] );
  return;

}


void cleanup( float *outdata, 
	      Asp_Data *Asp_Hist, 
	      regRegion *reg, 
	      Ardlib_Det_BadPix_Type **bpix)
{
  short bz;
  short close_ardlib=0;
  free( outdata );
  regFree( reg );

  for (bz=0;bz<10;bz++) {
    if ( bpix[bz] ) {
      ardlib_close_det_badpix( bpix[bz] );
      close_ardlib=1;
    }

    if (Asp_Hist[bz].nasp ) {
      free( Asp_Hist[bz].ra );
      free( Asp_Hist[bz].dec );
      free( Asp_Hist[bz].roll );
      free( Asp_Hist[bz].duration );
    }
  }
  pix_close_pixlib();
  if ( close_ardlib)  ardlib_finalize();
  return;
}


int load_params( char *infile, 
		char *ahstfile,
		char *outfile,
		char *geomfile,
		char *ardlib,
		short *clobber,
		short *verbose
		) 
{
  clgetstr( "infile", infile, 1000 );
  clgetstr( "asphistfile", ahstfile, 1000);
  clgetstr( "outfile", outfile, 1000 );
  clgetstr( "geompar", geomfile, 1000 );
  /*   clgetstr( "ardlibpar", ardlib, 1000 );*/
  *clobber = clgetb( "clobber" );
  *verbose = clgeti( "verbose" );

  if ( *verbose ) {


  }
  return(0);


}





int mktae()
{

  char infile[1000] = "/export/wmap_in_tdet_bin1.fits[wmap]";
  char ahstfile[1000] = "/export/acisf00635_000N001_c1_asphist3.fits";
  char outfile[1000] = "!outimg.fits";
  char geomfile[1000] = "geom";
  char ardlib[1000] = "ardlib";
  short clobber = 1;
  short verbose = 0;
  
  double sim[3];
  
  long laxes[2];
  dmDescriptor *phy_coord;
  regRegion *reg;
  Asp_Data Asp_Hist[10];
  Ardlib_Det_BadPix_Type *bpix[10];
  
  dmBlock *outBlock;
  float *outdata;
  
  memset( &Asp_Hist[0], 0, 10*sizeof(Asp_Data));

  load_params( infile, ahstfile, outfile, geomfile, ardlib, 
	       &clobber, &verbose );
  
  if ( 0 != ds_clobber( outfile, clobber, NULL ) ) {
    return(-1);
  }
  if ( NULL == ( outBlock = load_input_setup_output( infile, outfile, laxes,
						     &phy_coord, &reg))) {
    return(-1);
  }
  if ( 0 != load_asp_data( ahstfile, Asp_Hist, sim )) {
    return(-1);
  }

  if ( 0 != setup_pixlib_and_ardlib( infile, geomfile, ardlib, sim, bpix )) {
    return(-3);
  }

  if ( NULL == ( outdata = process_pixels( laxes, phy_coord, Asp_Hist, 
					   bpix, reg ))) {
    return(-5);
  }
  
  output_data( outBlock, outdata, laxes );
  ds_write_pixhist_in_dm( outBlock );

  dmImageClose( outBlock);
  
  cleanup( outdata, Asp_Hist, reg, bpix );

  return(0);

}
