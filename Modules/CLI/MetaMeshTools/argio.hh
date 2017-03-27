/*
 * author:  msturm,mstyner
 * created: 16 Apr 1997
 * changes:
 *
 * contains routines for command line parsing
 * (mostly copied from gaudi toolbox, changed to templated functions)
 *
 */

#ifndef __ARGIO_H__
#define __ARGIO_H__

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>


//Copied from itkMacro.h (7a870297f818f32a22fe720b30577cc58ce91634)
#ifndef itkMacro_h
#if __cplusplus >= 201103L
// In c++11 the override keyword allows you to explicity define that a function
// is intended to override the base-class version.  This makes the code more
// managable and fixes a set of common hard-to-find bugs.
#define ITK_OVERRIDE override
// In c++11 there is an explicit nullptr type that introduces a new keyword to
// serve as a distinguished null pointer constant: nullptr. It is of type
// nullptr_t, which is implicitly convertible and comparable to any pointer type
// or pointer-to-member type. It is not implicitly convertible or comparable to
// integral types, except for bool.
#define ITK_NULLPTR  nullptr
// In C++11 the throw-list specification has been deprecated,
// replaces with the noexcept specifier. Using this function
// specification adds the run-time check that the method does not
// throw, if it does throw then std::terminate will be called.
// Use cautiously.
#define ITK_NOEXCEPT noexcept
#else
#define ITK_OVERRIDE
#define ITK_NULLPTR  NULL
#define ITK_NOEXCEPT throw()
#endif
#endif




const int ipMAXTOKLEN = 255 ; // max. length of a token

// some workarounds for template function
typedef char* charp ;

inline float fatof( const char *str )
{
  return (float) atof( str ) ;
}

// command line parsing


template <class T>
inline T ipGetArgument( const char **argv , const char *keystr ,
                        T ( *convert ) ( const char *str ) , const T defval
                      )
{
  for( int i = 1 ; argv[ i ] ; i++ )
  {
    if( !strcmp( argv[ i ] , keystr ) )
    {
      if( argv[ i + 1 ] )
      {
        return convert( argv [ i + 1 ] ) ;
      }
      else
      {
        fprintf( stderr , "Error: ipGetArgument: argument value of option \"%s\" is missing!\n" ,
           argv[ i ] ) ;
        exit( -1 ) ;
      }
    }
  }
  return defval ;
}

inline charp ipGetStringArgument( const char **argv , const char keystr[] , const char defval[] )
{
  char *_defval = ( defval ? strdup( defval ) : ITK_NULLPTR ) ;
  for( int i = 1 ; argv[ i ] ; i++ )
  {
    if( !strcmp( argv[ i ] , keystr ) )
    {
      if( argv[ i + 1 ] )
      {
        return strdup( argv[ i + 1 ] ) ;
      }
      else
      {
        fprintf( stderr , "Error: ipGetArgument: argument value of option \"%s\" is missing!\n" ,
          argv[ i ] ) ;
        exit( -1 ) ;
      }
    }
  }
  return _defval ;
}

// reads in multiple string arguments started by keystr and ended by an '-' or
// end of args
inline int ipGetStringMultipArgument( const char **argv , const char *keystr , char **out , const int max )
{
  int i = 1 ;
  int num ;
  while( argv[ i ] )
  {
    if( !strcmp( argv[ i ] , keystr ) )
    {
      // keystr found
      if( argv[ i + 1 ] && ( argv[ i + 1 ] )[ 0 ] != '-' )
      {
        i++ ;
        num = 0 ;
        out[ 0 ] = strdup( argv[ i ] ) ;
        i++ ;
        num++ ;
        while( argv[ i ] && num < max && ( argv[ i ] )[ 0 ] != '-' )
        {
          out[ num ] = strdup( argv[ i ] ) ;
          i++ ;
          num++ ;
        }
        return num ;
        // jump back !
      }
      else
      {
        fprintf( stderr , "Error: ipGetArgument: argument value of option \"%s\" is missing!\n" , argv[ i ] ) ;
        exit( -1 ) ;
      }
    } // if keystr
    i++ ;
  } //while
  return 0 ;
}


inline int ipGetIntArgument( const char **argv , const char *keystr , const int defval )
{
  for( int i = 1 ; argv[ i ] ; i++ )
  {
    if( !strcmp( argv[ i ] , keystr ) )
    {
      if( argv[ i + 1 ] )
      {
        return atoi( argv[ i + 1 ] ) ;
      }
      else
      {
        fprintf( stderr ,"Error: ipGetArgument: argument value of option \"%s\" is missing!\n" ,
          argv[ i ] ) ;
        exit( -1 ) ;
      }
    }
  }
  return defval ;
}

inline float ipGetFloatArgument( const char **argv , const char *keystr , const float defval )
{
  for( int i = 1 ; argv[ i ] ; i++ )
  {
    if( !strcmp( argv[ i ] , keystr ) )
    {
      if( argv[ i + 1 ] )
      {
        return fatof( argv[ i + 1 ] ) ;
      }
      else
      {
        fprintf( stderr , "Error: ipGetArgument: argument value of option \"%s\" is missing!\n",
          argv[ i ] ) ;
        exit( -1 ) ;
      }
    }
  }
  return defval ;
}

inline double ipGetDoubleArgument( const char **argv , const char *keystr , const double defval )
{
  for( int i = 1 ; argv[ i ] ; i++ )
  {
    if( !strcmp( argv[ i ] , keystr ) )
    {
      if( argv[ i + 1 ] )
      {
        return atof( argv[ i + 1 ] ) ;
      }
      else
      {
        fprintf( stderr , "Error: ipGetArgument: argument value of option \"%s\" is missing!\n" ,
          argv[ i ] ) ;
        exit( -1 ) ;
      }
    }
  }
  return defval ;
}

inline int ipExistsArgument( const char **argv , const char *keystr )
{
  for( int i = 1 ; argv[ i ] ; i++ )
  {
    if( !strcmp( argv[ i ] , keystr ) )
    {
      return 1 ;
    }
  }
  return 0 ;
}


// string utilities
// appends src to resized dst and returns it
inline char *ipAppendString( char *&dst , const char *src )
{
  if( dst )
  {
    dst = (char *) realloc( dst , strlen( dst ) + strlen( src ) + 1 ) ;
    return strcat( dst , src ) ;
  }
  else
  {
    return strdup( src ) ;
  }
}


/*************************************************************
 * counts the words on one line
 ************************************************************/
inline int ipLineWordCount( const char *s )
{
  int n = 0 ;
  for( ; ; )
  {
    while( isspace( *s ) )
    {
      s++ ;
    }
    if( *s == '\0' )
    {
      return n ;
    }
    while( isgraph( *s ) )
    {
      s++ ;
    }
    n++ ;
  }
}


/*************************************************************
 * gets the generic name of a file cutting the extension
 ************************************************************/
inline char *ipGetBaseName( const char *string )
{
  unsigned int i ;
  char *ret = ITK_NULLPTR ;
  char *retp = ITK_NULLPTR ;
  if( !( ret = strdup( string ) ) )
  {
    fprintf( stderr , "Error: ipGetBaseName [%s, line %d]: strdup() failed:" ,
            __FILE__ ,  __LINE__ ) ;
    perror( "" ) ;
    exit( errno ) ;
  }
  for( i = 0 , retp = ret ; i < strlen( string ) ; i++ , retp++ )
  {
    if( *retp == '.' )
    {
      *retp = '\0' ;
      break ;
    }
  }
  return ret ;
}


/*************************************************************
 * reads a line into s (already allocated) and returns length
 * lim describes the max line length.
 * function proposed by Kernighan and Ritchie
 ************************************************************/
inline int ipfgetline( FILE* f , char* s , int lim )
{
  int c ;
  int i ;
  for( i = 0 ; ( i < lim - 1 ) && ( ( c = getc( f ) ) != EOF ) && ( c != '\n' ) ; ++i )
  {
    s[ i ] = c ;
  }
  if( c == '\n' )
  {
    s[ i ] = c ;
    ++i ;
  }
  s[ i ] = '\0' ;
  return ( i ) ;
}

// extracts string tokens from a string which are either separeted by
//    - white-space (isspace()) or
//    - punctuation (ispunct()) except '.', '+', '-', '_'
// and converts them to type <class T>.
//
// at most n tokens will be extracted, tokenval[] has to be
// allocated
// returns number of tokens found, converted tokens in tokenval[]

template <class T>
int ipExtractTokens( T *tokenval , const char *tokenstr , const int n ,
                     T (*convert) ( const char *str )
                   )
{
  char *tmp_token = new char[ ipMAXTOKLEN ] ;
  char *tmp_tokenp = tmp_token ;
  const char *tokenp = tokenstr ;
  int  i = 0 ;
  while( ( i < n ) && ( *tokenp ) )
  {
    memset( tmp_token , '\0' , ipMAXTOKLEN * sizeof( char ) ) ;
    while( isspace( *tokenp ) && *tokenp )
    {
      tokenp++ ;
    }
    tmp_tokenp = tmp_token ;
    while( (*tokenp) && ( isalnum(*tokenp)
                          || (*tokenp == '.' )
                          || (*tokenp == '-' )
                          || (*tokenp == '+' )
                          || (*tokenp == '_' )
                        )
            && ( (tmp_tokenp - tmp_token) < ipMAXTOKLEN )
         )
    {
      *( tmp_tokenp++ ) = *( tokenp++ ) ;
    }
    if( *tokenp )
    {
      tokenp++ ; // skip separator
    }
    tokenval[ i++ ] = convert( tmp_token ) ;
  }
  delete [] tmp_token ;
  return i ;
}

// extracts string tokens from a string which are either separeted by
//    - white-space (isspace()) or
//    - commas
// and converts them to type <class T>.
//
// at most n tokens will be extracted, tokenval[] has to be
// allocated
// returns number of tokens found, converted tokens in tokenval[]

template <class T>
int ipExtractSpaceSepTokens( T *tokenval , const char *tokenstr , const int n ,
          T ( *convert ) ( const char *str ) )
{
  char *tmp_token = new char [ ipMAXTOKLEN ] ;
  char *tmp_tokenp = tmp_token ;
  const char *tokenp = tokenstr ;
  int  i = 0 ;
  while( (i < n ) && ( *tokenp ) )
  {
    memset( tmp_token , '\0', ipMAXTOKLEN * sizeof(char) ) ;
    while( isspace(*tokenp) && *tokenp )
    {
      tokenp++ ;
    }
    tmp_tokenp = tmp_token ;
    while( (*tokenp) && ( isspace(*tokenp)
                          || ( *tokenp != ',' )
                        )
                     && ( (tmp_tokenp - tmp_token ) < ipMAXTOKLEN )
         )
    {
      *(tmp_tokenp++) = *(tokenp++) ;
    }
    if( *tokenp )
    {
      tokenp++ ; // skip separator
    }
    tokenval[ i++ ] = convert( tmp_token ) ;
  }
  delete [] tmp_token ;
  return i ;
}


inline int ipExtractIntTokens( int *tokenval , const char *tokenstr , const int n )
{
  char *tmp_token = new char [ ipMAXTOKLEN ] ;
  char *tmp_tokenp = tmp_token ;
  const char *tokenp = tokenstr ;
  int i = 0 ;
  while( (i < n) && (*tokenp) )
  {
    memset( tmp_token , '\0' , ipMAXTOKLEN * sizeof(char) ) ;
    while( isspace(*tokenp) && *tokenp )
    {
      tokenp++ ;
    }
    tmp_tokenp = tmp_token ;
    while( (*tokenp)
           && ( isalnum(*tokenp)
                || (*tokenp == '.')
                || (*tokenp == '-')
                || (*tokenp == '+')
                || (*tokenp == '_')
              )
           && ( ( tmp_tokenp - tmp_token ) < ipMAXTOKLEN )
         )
    {
      *(tmp_tokenp++) = *(tokenp++) ;
    }
    if( *tokenp )
    {
      tokenp++ ; // skip separator
    }
    tokenval[ i++ ] = atoi( tmp_token ) ;
  }
  delete [] tmp_token ;
  return i ;
}

inline int ipExtractFloatTokens( float *tokenval , const char *tokenstr , const int n )
{
  char *tmp_token = new char [ ipMAXTOKLEN ] ;
  char *tmp_tokenp = tmp_token ;
  const char *tokenp = tokenstr ;
  int  i = 0 ;
  while( (i < n) && (*tokenp) )
  {
    memset( tmp_token , '\0' , ipMAXTOKLEN * sizeof(char) ) ;
    while( isspace(*tokenp) && *tokenp )
    {
      tokenp++ ;
    }
    tmp_tokenp = tmp_token ;
    while( (*tokenp)
           && ( isalnum( *tokenp )
                || ( *tokenp == '.' )
                || ( *tokenp == '-' )
                || ( *tokenp == '+' )
                || ( *tokenp == '_' )
              )
           && ( (tmp_tokenp - tmp_token ) < ipMAXTOKLEN )
         )
    {
      *( tmp_tokenp++ ) = *( tokenp++ ) ;
    }
    if( *tokenp )
    {
      tokenp++ ; // skip separator
    }
    tokenval[ i++ ] = fatof( tmp_token ) ;
  }
  delete [] tmp_token ;
  return i ;
}

inline int ipExtractDoubleTokens( double *tokenval , const char *tokenstr , const int n )
{
  char *tmp_token = new char [ ipMAXTOKLEN ] ;
  char *tmp_tokenp = tmp_token ;
  const char *tokenp = tokenstr ;
  int  i = 0 ;
  while( (i < n) && (*tokenp) )
  {
    memset( tmp_token , '\0' , ipMAXTOKLEN * sizeof(char) ) ;
    while( isspace(*tokenp) && *tokenp )
    {
      tokenp++ ;
    }
    tmp_tokenp = tmp_token ;
    while( (*tokenp)
           && ( isalnum( *tokenp )
                || ( *tokenp == '.' )
                || ( *tokenp == '-' )
                || ( *tokenp == '+' )
                || ( *tokenp == '_' )
              )
           && ( ( tmp_tokenp - tmp_token ) < ipMAXTOKLEN )
         )
    {
      *( tmp_tokenp++ ) = *( tokenp++ ) ;
    }
    if( *tokenp )
    {
      tokenp++ ; // skip separator
    }
    tokenval[ i++ ] = atof( tmp_token ) ;
  }
  delete [] tmp_token ;
  return i ;
}

inline int ipExtractStringTokens( char **tokenval , const char *tokenstr , const int n )
{
  char *tmp_token = new char [ ipMAXTOKLEN ] ;
  char *tmp_tokenp = tmp_token ;
  const char *tokenp = tokenstr ;
  int i = 0 ;

  while( (i < n) && (*tokenp) )
  {
    memset( tmp_token , '\0' , ipMAXTOKLEN * sizeof(char) ) ;
    while( isspace(*tokenp) && *tokenp )
    {
      tokenp++ ;
    }
    tmp_tokenp = tmp_token ;
    while( (*tokenp)
           && ( isalnum( *tokenp )
                || ( *tokenp == '.' )
                || ( *tokenp == '-' )
                || ( *tokenp == '+' )
                || ( *tokenp == '_' )
                || ( *tokenp == '/' )
              )
           && ( (tmp_tokenp - tmp_token ) < ipMAXTOKLEN )
         )
    {
      *( tmp_tokenp++ ) = *( tokenp++ ) ;
    }
    if( *tokenp )
    {
      tokenp++ ; // skip separator
    }
    tokenval[ i++ ] = strdup( tmp_token ) ;
  }
  delete [] tmp_token ;
  return i ;
}

#endif
