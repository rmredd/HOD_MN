/*   	ftwrite writes data unformatted using fortran convention --
	i.e. an integer specifying the number of bytes in the record,
	the data record, and another integer specifying the number of 
	bytes in the record.  The call is identical to the standard
	i/o library routine fwrite.
*/

#include <stdio.h>

int ftwrite(ptr,size,nitems,stream)
char *ptr ;
unsigned size, nitems ;
FILE *stream ;

{
    int nbytes, nitem1 ;
    int errno ;

    errno = 0 ;
    nbytes = size*nitems ;
    if ( fwrite(&nbytes,sizeof(int),1,stream) != 1 )  {
	errno = -10 ;
	fprintf(stderr,"write error, is the file open ? \n") ;
	}
    nitem1 = fwrite(ptr,size,nitems,stream) ;
    if ( nitem1 != nitems ) {
	errno = -20 ;
	fprintf(stderr,"write error, %d items requested, %d items written. \n",
		nitems,nitem1) ;
    }
    if ( fwrite(&nbytes,sizeof(int),1,stream) != 1 )  {
	errno = -30 ;
	fprintf(stderr,"write error on second byte label \n") ;
    }
    
    return(errno) ;

}
