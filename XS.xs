#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#define GEOCALC_STASH geocalc_stash

#if __GNUC__ >= 3
# define expect(expr,value)         __builtin_expect ((expr), (value))
# define INLINE                     static inline
#else
# define expect(expr,value)         (expr)
# define INLINE                     static
#endif

#define PI 3.14159265358979323846
#define DEG2RAD(DEG) ((DEG)*((PI)/(180.0)))
#define RAD2DEG(DEG) ((DEG)*(180.0/PI))

#define INGC HV

static HV *geocalc_stash;

typedef struct {
  /* Coords */
  long double latitude;
  long double longitude;

  /* Distance Units */
  /* 0 -> m, 1 -> km, 2 -> yards, 3 -> feet, 4 -> mile */
  int unit_conv;

  /* Earth Radius */
  long double radius;
} GCX;

typedef struct {
    long double lat;
    long double lon;
    long double final_bearing;
} DESTINATION;

INLINE void
geocalc_init( GCX *gcx, HV * options )
{
  Zero( gcx, 1, GCX );
  gcx->latitude  = (long double) SvNV( *hv_fetch(options, "lat", 3, 0) );
  gcx->longitude = (long double) SvNV( *hv_fetch(options, "lon", 3, 0) );

  SV ** sv = hv_fetch(options, "units", 5, 0);
  if( sv == (SV**)NULL ) {
    gcx->unit_conv = 0;    // Default to m
  } else {
    if( strEQ( SvPV_nolen( *sv ), "m" ) ) {
      gcx->unit_conv = 0;
    } else if( strEQ( SvPV_nolen( *sv ), "k-m" ) ) {
      gcx->unit_conv = 1;
    } else if( strEQ( SvPV_nolen( *sv ), "yd" ) ) {
      gcx->unit_conv = 2;
    } else if( strEQ( SvPV_nolen( *sv ), "ft" ) ) {
      gcx->unit_conv = 3;
    } else if( strEQ( SvPV_nolen( *sv ), "mi" ) ) {
      gcx->unit_conv = 4;
    } else if( strEQ( SvPV_nolen( *sv ), "" ) ) {
      gcx->unit_conv = 0;
    } else {
      warn("Unrecognised unit (defaulting to m)");
      gcx->unit_conv = 0;
    }
  }
  /*
  sv = hv_fetch(options, "radius", 6, 0);
  if( sv == (SV**)NULL ) {
    gcx->radius = 6371;    // Default to KM
  } else {
    gcx->radius = 6371;    // Default to KM
  }
  */
  gcx->radius    = 6371; // Earth Radius
}

INLINE long double
convert_km( long double input, int to_unit )
{
  double output = 0;

  switch( to_unit )
  {
    case 0  : output = input * 1000;
              break;
    case 1  : output = input;               // kilometer
              break;
    case 2  : output = input * 1093.6133;   // yard
              break;
    case 3  : output = input * 3280.8399;   // feet
              break;
    case 4  : output = input * 0.621371192; // mile
              break;
    default : output = input;               // kilometer
              break;
  }

  return output;
}

INLINE long double
convert_to_m( long double input, int from_unit )
{
  double output = 0;

  switch( from_unit )
  {
    case 0  : output = input;            // m
              break;
    case 1  : output = input * 1000;     // km
              break;
    case 2  : output = input * 0.9144;   // yard
              break;
    case 3  : output = input * 0.3048;   // yard
              break;
    case 4  : output = input * 1609.344; // mile
              break;
    default : output = input * 1000;     // kilometer
              break;
  }

  return output;
}

INLINE long double
cint( double x ) {
    double f;
    if ( modf( x, &f ) >= 0.5 )
        return ( x >= 0 ) ? ceil( x ) : floor( x );
    else
        return ( x < 0 ) ? ceil( x ) : floor( x );
}

INLINE long double
_precision_ld( long double r, int places ) {
    long double off = pow( 10, -places );
    return cint( r * off ) / off;
}

INLINE long double
_ib_precision( long double brng, int precision, int mul ) {
    return _precision_ld( (long double)fmod( mul * ( RAD2DEG( brng ) ) + 360, 360 ), precision );
}

INLINE long double
_fb_precision( long double brng, int precision ) {
    return _precision_ld( (long double)fmod( ( RAD2DEG( brng ) ) + 180, 360 ), precision );
}


INLINE int
number_of_decimals( SV * sv )
{
    int precision = -6;
    if( SvOK( sv ) )
        precision = SvIV( sv );
    else
        precision = -6;

    return precision;
}

INLINE void
_destination_point( GCX *self, double bearing, double s, int precision, DESTINATION *dest )
{
    s = convert_to_m( s, self->unit_conv );

    long double r_major    = 6378137;           // Equatorial Radius, WGS84
    long double r_minor    = 6356752.314245179; // defined as constant
    long double f          = 1/298.257223563;   // 1/f = ( $r_major - $r_minor ) / $r_major

    long double alpha1     = DEG2RAD( bearing );
    long double sinAlpha1  = sin( alpha1 );
    long double cosAlpha1  = cos( alpha1 );

    long double tanU1      = ( 1 - f ) * tan( DEG2RAD( self->latitude ) );
    long double cosU1      = 1 / sqrt( ( 1 + ( tanU1 * tanU1 ) ) );
    long double sinU1      = tanU1 * cosU1;
    long double sigma1     = atan2( tanU1, cosAlpha1 );
    long double sinAlpha   = cosU1 * sinAlpha1;
    long double cosSqAlpha = 1 - ( sinAlpha * sinAlpha );

    long double uSq        = cosSqAlpha * ( ( r_major * r_major ) - ( r_minor * r_minor ) ) / ( r_minor * r_minor );
    long double A          = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    long double B          = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

    long double sigma      = s / ( r_minor * A );
    long double sigmaP     = PI * 2;

    long double cos2SigmaM = cos(2*sigma1 + sigma);
    long double sinSigma   = sin(sigma);
    long double cosSigma   = cos(sigma);

    while ( abs(sigma - sigmaP) > 1 * pow( 10, -12 ) ) {
        cos2SigmaM = cos(2*sigma1 + sigma);
        sinSigma   = sin(sigma);
        cosSigma   = cos(sigma);

        long double deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
                  B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
        sigmaP                 = sigma;
        sigma                  = s / (r_minor*A) + deltaSigma;
    }

    long double tmp    = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1;
    long double lat2   = atan2( sinU1*cosSigma + cosU1*sinSigma*cosAlpha1, (1-f)*sqrt(sinAlpha*sinAlpha + tmp*tmp) );

    long double lambda = atan2(sinSigma*sinAlpha1, cosU1*cosSigma - sinU1*sinSigma*cosAlpha1);
    long double C      = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
    long double L      = lambda - (1-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));

    long double lon2   = fmod( DEG2RAD( self->longitude ) + L + ( PI * 3 ), PI * 2 ) - PI;
    long double revAz  = atan2( sinAlpha, -tmp ); // final bearing, if required

    dest->lat = _precision_ld( RAD2DEG( lat2 ), precision );
    dest->lon = _precision_ld( RAD2DEG( lon2 ), precision );
    dest->final_bearing = _precision_ld( RAD2DEG( revAz ), precision );
}


MODULE = Geo::Calc::XS         PACKAGE = Geo::Calc::XS
BOOT:
{
    geocalc_stash = gv_stashpv( "Geo::Calc::XS", 1 );
}

PROTOTYPES: DISABLE

void CLONE (...)
    CODE:
        geocalc_stash = 0;

void new ( char *klass, ... )
    PREINIT:
        int add_count = items - 1;
    PPCODE:
    {
        SV *pv = NEWSV ( 0, sizeof( GCX ) );
        HV *options = newHV();
        int i = 0;
        SvPOK_only( pv );
        if( add_count % 2 != 0 )
            croak( "Please check your parameters while initiating the module\n" );

        for( i = 0; i < add_count; i = i + 2 ) {
            hv_store_ent( options, ST( i + 1 ), newSVsv(ST( i + 2)), 0);
        }

        geocalc_init( (GCX *)SvPVX( pv ), options );

        SvREFCNT_dec((SV *) options);

        XPUSHs( sv_2mortal( sv_bless(
           newRV_noinc( pv ),
           strEQ( klass, "Geo::Calc::XS" ) ? GEOCALC_STASH : gv_stashpv( klass, 1 )
        ) ) );
    }

void get_lat ( GCX *self, ... )
    PPCODE:
    {
        XPUSHs( sv_2mortal( newSVnv( self->latitude ) ) );
    }

void get_lon ( GCX *self, ... )
    PPCODE:
    {
        XPUSHs( sv_2mortal( newSVnv( self->longitude ) ) );
    }

void get_units ( GCX *self, ... )
    PPCODE:
    {
        /* Distance Units */
        /* 0 -> m, 1 -> km, 2 -> yards, 3 -> feet, 4 -> mile */
        switch ( self->unit_conv)
        {
            case 0:
                XPUSHs( sv_2mortal( newSVpv( "m", 0 ) ) );
                break;
            case 1:
                XPUSHs( sv_2mortal( newSVpv( "k-m", 0 ) ) );
                break;
            case 2:
                XPUSHs( sv_2mortal( newSVpv( "yd", 0 ) ) );
                break;
            case 3:
                XPUSHs( sv_2mortal( newSVpv( "ft", 0 ) ) );
                break;
            default:
                XPUSHs( sv_2mortal( newSVpv( "mi", 0 ) ) );
        }
    }

void get_radius ( GCX *self, ... )
    PPCODE:
    {
        XPUSHs( sv_2mortal( newSVnv( self->radius ) ) );
    }

void distance_to( GCX *self, INGC *to_latlon, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 2 ) );

        long double lat1 = DEG2RAD( self->latitude );
        long double lon1 = DEG2RAD( self->longitude );
        long double lat2 = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) );
        long double lon2 = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lon", 3, 0) ) );

        long double t = pow( sin( ( lat2 - lat1 ) / 2 ), 2 ) + ( pow( cos( lat1 ), 2 ) * pow( sin( ( lon2 - lon1 )/2 ), 2 ) );
        long double d = convert_km( self->radius * ( 2 * atan2( sqrt(t), sqrt(1-t) ) ), self->unit_conv );

        XPUSHs( sv_2mortal( newSVnv( _precision_ld( d, precision ) ) ) );
    }

void destination_point( GCX *self, double bearing, double s, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 3 ) );
        DESTINATION dest;
        _destination_point( self, bearing, s, precision, &dest );

        HV *retval = newHV();
        hv_store( retval, "lat", 3, newSVnv( dest.lat ), 0 );
        hv_store( retval, "lon", 3, newSVnv( dest.lon ), 0 );
        hv_store( retval, "final_bearing", 13, newSVnv( dest.final_bearing ), 0 );
        XPUSHs( sv_2mortal( newRV_noinc( (SV*)retval ) ) );
    }

void boundry_box( GCX *self, double width, ... )
    PPCODE:
    {
        int precision = -6;
        double height = 0;
        DESTINATION dest;

        if( SvOK( ST( 2 ) ) ) {
            height = SvNV( ST( 2 ) );
        } else {
            width = width * 2;
            height = width;
        }

        if( SvOK( ST( 3 ) ) )
            precision = number_of_decimals( ST( 3 ) );

        HV *retval = newHV();

        _destination_point( self, 180, height / 2, precision, &dest );
        hv_store( retval, "lat_min", 7, newSVnv( dest.lat ), 0 );

        _destination_point( self, 270, width  / 2, precision, &dest );
        hv_store( retval, "lon_min", 7, newSVnv( dest.lon ), 0 );

        _destination_point( self,   0, height / 2, precision, &dest );
        hv_store( retval, "lat_max", 7, newSVnv( dest.lat ), 0 );

        _destination_point( self,  90, width  / 2, precision, &dest );
        hv_store( retval, "lon_max", 7, newSVnv( dest.lon ), 0 );

        XPUSHs( sv_2mortal( newRV_noinc( (SV*)retval ) ) );
    }

void midpoint_to( GCX *self, INGC *to_latlon, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 2 ) );

        long double lat1 = DEG2RAD( self->latitude );
        long double lon1 = DEG2RAD( self->longitude );

        long double lat2 = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) );
        long double dlon = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lon", 3, 0) ) - self->longitude );

        long double bx = cos( lat2 ) * cos( dlon );
        long double by = cos( lat2 ) * sin( dlon );

        long double lat3 = atan2( sin( lat1 ) + sin ( lat2 ), sqrt( ( pow( ( cos( lat1 ) + bx ), 2 ) ) + pow( by, 2 ) ) );
        long double lon3 = fmod( lon1 + atan2( by, cos( lat1 ) + bx ) + ( PI * 3 ), PI * 2 ) - PI;

        HV *retval = newHV();
        hv_store( retval, "lat", 3, newSVnv( _precision_ld( RAD2DEG( lat3 ), precision ) ), 0 );
        hv_store( retval, "lon", 3, newSVnv( _precision_ld( RAD2DEG( lon3 ), precision ) ), 0 );

        XPUSHs( sv_2mortal( newRV_noinc( (SV*)retval ) ) );
    }

void intersection( GCX *self, double brng1, INGC *to_latlon, double brng2, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 4 ) );

        long double lat1   = DEG2RAD( self->latitude );
        long double lon1   = DEG2RAD( self->longitude );
        long double lat2   = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) );
        long double lon2   = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lon", 3, 0) ) );
        long double brng13 = DEG2RAD( brng1 );
        long double brng23 = DEG2RAD( brng2 );

        long double dlat   = lat2 - lat1;
        long double dlon   = lon2 - lon1;

        long double dist12 = 2 * asin( sqrt( pow( sin( dlat/2 ), 2 ) + cos( lat1 ) * cos( lat2 ) * pow( sin( dlon/2 ), 2 ) ) );
        if( dist12 == 0 ) {
            XPUSHs( &PL_sv_undef );
            return;
        }

        // initial/final bearings between points
        long double brnga, brngb;
        if( sin( dist12 ) * cos( lat1 ) > 0 ) {
            brnga = acos( ( sin( lat2 ) - sin( lat1 ) * cos( dist12 ) ) / ( sin( dist12 ) * cos( lat1 ) ) );
        } else {
            brnga = 0;
        }

        if( sin( dist12 ) * cos( lat2 ) > 0 ) {
            brngb = acos( ( sin( lat1 ) - sin( lat2 ) * cos( dist12 ) ) / ( sin( dist12 ) * cos( lat2 ) ) );
        } else {
            brngb = 0;
        }

        long double brng12, brng21;
        if( sin( dlon ) > 0 ) {
            brng12 = brnga;
            brng21 = PI*2 - brngb;
        } else {
            brng12 = PI*2 - brnga;
            brng21 = brngb;
        }

        long double alpha1 = fmod( brng13 - brng12 + ( PI * 3 ), PI * 2 ) - PI;
        long double alpha2 = fmod( brng21 - brng23 + ( PI * 3 ), PI * 2 ) - PI;

        if( ( sin( alpha1 ) == 0 ) && ( sin( alpha2 ) == 0 ) ) { // infinite intersections
            XPUSHs( &PL_sv_undef );
            return;
        }

        if( sin( alpha1 ) * sin( alpha2 ) < 0 ) { // ambiguous intersection
            XPUSHs( &PL_sv_undef );
            return;
        }

        long double alpha3 = acos( -cos( alpha1 ) * cos( alpha2 ) + sin( alpha1 ) * sin( alpha2 ) * cos( dist12 ) );
        long double dist13 = atan2( sin( dist12 ) * sin( alpha1 ) * sin( alpha2 ), cos( alpha2 ) + cos( alpha1 ) * cos( alpha3 ) );
        long double lat3 = asin( sin( lat1 ) * cos( dist13 ) + cos( lat1 ) * sin( dist13 ) * cos( brng13 ) );
        long double dlon13 = atan2( sin( brng13 ) * sin( dist13 ) * cos( lat1 ), cos( dist13 ) - sin( lat1 ) * sin( lat3 ) );
        long double lon3 = fmod( lon1 + dlon13 + ( PI * 3 ), PI * 2 ) - PI;

        HV *retval = newHV();
        hv_store( retval, "lat", 3, newSVnv( _precision_ld( RAD2DEG( lat3 ), precision ) ), 0 );
        hv_store( retval, "lon", 3, newSVnv( _precision_ld( RAD2DEG( lon3 ), precision ) ), 0 );

        XPUSHs( sv_2mortal( newRV_noinc( (SV*)retval ) ) );
    }

void distance_at( GCX *self, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 1 ) );

        long double lat = DEG2RAD( self->latitude );

        // Set up "Constants"
        double m1 = 111132.92; // latitude calculation term 1
        double m2 = -559.82;   // latitude calculation term 2
        double m3 = 1.175;     // latitude calculation term 3
        double m4 = -0.0023;   // latitude calculation term 4
        double p1 = 111412.84; // longitude calculation term 1
        double p2 = -93.5;     // longitude calculation term 2
        double p3 = 0.118;     // longitude calculation term 3.

        HV *retval = newHV();
        hv_store( retval, "m_lat", 5, newSVnv( _precision_ld( m1 + (m2 * cos(2 * lat)) + (m3 * cos(4 * lat)) + ( m4 * cos(6 * lat) ), precision ) ), 0 );
        hv_store( retval, "m_lon", 5, newSVnv( _precision_ld( ( p1 * cos(lat)) + (p2 * cos(3 * lat)) + (p3 * cos(5 * lat) ), precision ) ), 0 );

        XPUSHs( sv_2mortal( newRV_noinc( (SV*)retval ) ) );
    }

void bearing_to( GCX *self, INGC *to_latlon, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 2 ) );

        long double lat1 = DEG2RAD( self->latitude );
        long double lat2 = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) );
        long double dlon = DEG2RAD( self->longitude - (long double) SvNV( *hv_fetch(to_latlon, "lon", 3, 0) ) );

        long double brng = atan2( sin( dlon ) * cos( lat2 ), ( cos( lat1 ) * sin( lat2 ) ) - ( sin( lat1 ) * cos( lat2 ) * cos( dlon ) ) );

        XPUSHs( sv_2mortal( newSVnv( _ib_precision( brng, precision, -1 ) ) ) );
    }

void final_bearing_to( GCX *self, INGC *to_latlon, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 2 ) );

        long double lat1 = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) );
        long double lat2 = DEG2RAD( self->latitude );
        long double dlon = -DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lon", 3, 0) ) - self->longitude );

        long double brng = atan2( sin( dlon ) * cos( lat2 ), ( cos( lat1 ) * sin( lat2 ) ) - ( sin( lat1 ) * cos( lat2 ) * cos( dlon ) ) );

        XPUSHs( sv_2mortal( newSVnv( _fb_precision( brng, precision ) ) ) );
    }

void rhumb_distance_to( GCX *self, INGC *to_latlon, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 2 ) );

        long double lat1 = DEG2RAD( self->latitude );
        long double lat2 = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) );
        long double dlat = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) - self->latitude );
        long double dlon = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lon", 3, 0) ) - self->longitude );
        if( dlon < 0 )
            dlon = -dlon;

        long double dphi = log( tan( lat2/2 + PI/4 ) / tan( lat1/2 + PI/4 ) );

        long double q;

        if ( dphi != 0 ) {
            q = dlat/dphi;
        } else {
            q = cos(lat1); // E-W line gives dPhi=0
        }

        if( dlon > PI )
            dlon = PI*2 - dlon;

        long double dist = sqrt( pow( dlat, 2 ) + pow( q, 2 ) * pow( dlon, 2 ) ) * self->radius;

        XPUSHs( sv_2mortal( newSVnv( _precision_ld( convert_km( dist, self->unit_conv ), precision ) ) ) );
}

void rhumb_bearing_to( GCX *self, INGC *to_latlon, ... )
    PPCODE:
    {
        int precision = number_of_decimals( ST( 2 ) );

        long double lat1 = DEG2RAD( self->latitude );
        long double lat2 = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lat", 3, 0) ) );
        long double dlon = DEG2RAD( (long double) SvNV( *hv_fetch(to_latlon, "lon", 3, 0) ) - self->longitude );

        long double dphi = log( tan( lat2/2 + PI/4 ) / tan( lat1/2 + PI/4 ) );
        long double abs_dphi = dphi;
        if( abs_dphi < 0 ) {
            abs_dphi = -abs_dphi;
        }

        if( abs_dphi > PI ) {
            dlon = ( dlon > 0 ) ? -( PI*2 - dlon ) : ( PI*2 + dlon );
        }

        XPUSHs( sv_2mortal( newSVnv( _ib_precision( atan2( dlon, dphi ), precision, 1 ) ) ) );
    }

void rhumb_destination_point( GCX *self, double brng, double s, ... )
    PPCODE:
    {
        int precision    = number_of_decimals( ST( 3 ) );
        long double d    = ( convert_to_m( s, self->unit_conv ) / 1000 ) / self->radius;

        long double lat1 = DEG2RAD( self->latitude );
        long double lon1 = DEG2RAD( self->longitude );
        brng             = DEG2RAD( brng );

        long double lat2 = lat1 + ( d * cos( brng ) );
        long double dlat = lat2 - lat1;
        long double dphi = log( tan( lat2/2 + PI/4 ) / tan( lat1/2 + PI/4 ) );
        long double q    = ( dphi != 0 ) ? dlat/dphi : cos(lat1); //# E-W line gives dPhi=0
        long double dlon = d * sin( brng ) / q;
        // check for some daft bugger going past the pole
        long double lat2_abs = lat2;
        if( lat2_abs < 0 )
            lat2_abs = -lat2_abs;

        if( lat2_abs > PI/2 )
            lat2 = ( lat2 ) > 0 ? PI-lat2 : -(PI-lat2);

        long double lon2 = fmod( lon1 + dlon + ( PI * 3 ), PI * 2 ) - PI;

        HV *retval = newHV();
        hv_store( retval, "lat", 3, newSVnv( _precision_ld( RAD2DEG( lat2 ), precision ) ), 0 );
        hv_store( retval, "lon", 3, newSVnv( _precision_ld( RAD2DEG( lon2 ), precision ) ), 0 );

        XPUSHs( sv_2mortal( newRV_noinc( (SV*)retval ) ) );
    }
