#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include "tsvio.h"

/***************************************************************************
 *  Description:
 *      Read next tab-separated field
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     tsv_read_field(const char *argv[], FILE *infile,
		       char buff[], size_t buff_size)

{
    size_t  c;
    char    *p;
    
    c = 0, p = buff;
    while ( (c < buff_size) && ((*p = getc(infile)) != '\t') &&
	    (*p != '\n') && (*p != EOF) )
	++c, ++p;
    
    if ( c == buff_size )
    {
	fprintf(stderr, "%s: Buffer overflow reading field.\n", argv[0]);
	fputs(buff, stderr);
	exit(EX_DATAERR);
    }
    
    if ( *p == EOF )
	return 0;
    else
    {
	*p = '\0';
	return c;
    }
}


/***************************************************************************
 *  Description:
 *      Discard next field
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    tsv_skip_field(const char *argv[], FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != '\t') && (ch != '\n') && (ch != EOF) )
	;
    
    if ( ch == EOF )
    {
	fprintf(stderr, "%s: EOF encounterd unexpectedly.\n", argv[0]);
	exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Description:
 *      Discard the rest of the current input line.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    tsv_skip_rest_of_line(const char *argv[], FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != EOF) && (ch != '\n') )
	;
    if ( ch == EOF )
    {
	fprintf(stderr, "%s: EOF reached unexpectedly.\n", argv[0]);
	exit(EX_DATAERR);
    }
}
