package Geo::Calc::XS;

require 5.4.0;

use strict;
use warnings;
use utf8;

use Exporter;
use XSLoader;

our @ISA = qw( Exporter DynaLoader );

our %EXPORT_TAGS = ( 'all' => [ 'new', 'distance_to' ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = ();
our $VERSION = '0.25';

XSLoader::load 'Geo::Calc::XS', $VERSION;

# Copyrights 2011-2012 by Sorin Alexandru Pop.
# For other contributors see ChangeLog.
# See the manual pages for details on the licensing terms.

=head1 NAME

Geo::Calc::XS - simple geo calculator for points and distances

=head1 SYNOPSIS

 use Geo::Calc::XS;

 my $gc            = Geo::Calc::XS->new( lat => 40.417875, lon => -3.710205 );
 my $distance      = $gc->distance_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $brng          = $gc->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $f_brng        = $gc->final_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $midpoint      = $gc->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $destination   = $gc->destination_point( 90, 1, -6 );
 my $bbox          = $gc->boundry_box( 3, 4, -6 );
 my $r_distance    = $gc->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $r_brng        = $gc->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $r_destination = $gc->rhumb_destination_point( 30, 1, -6 );
 my $point         = $gc->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180, -6 );

=head1 DESCRIPTION

C<Geo::Calc::XS> implements a variety of calculations for latitude/longitude points

All these formulare are for calculations on the basis of a spherical earth
(ignoring ellipsoidal effects) which is accurate enough* for most purposes.

[ In fact, the earth is very slightly ellipsoidal; using a spherical model
gives errors typically up to 0.3% ].

Benchmarnking this module and Geo::Calc I found out that this module is sometime
more than 8000 times faster.

=head1 Geo::Calc::XS->new()

 $gc = Geo::Calc::XS->new( lat => 40.417875, lon => -3.710205 ); # Somewhere in Madrid
 $gc = Geo::Calc::XS->new( lat => 51.503269, lon => 0, units => 'k-m' ); # The O2 Arena in London

Creates a new Geo::Calc::XS object from a latitude and longitude. The default
deciaml precision is -6 for all functions => meaning by default it always
returns the results with 6 deciamls.

The default unit distance is 'm' (meter), but you cand define another unit using 'units'.
Accepted values are: 'm' (meters), 'k-m' (kilometers), 'yd' (yards), 'ft' (feet) and 'mi' (miles)

Returns ref to a Geo::Calc::XS object.

=head2 Parameters

=over 4

=item lat

C<>=> latitude of the point ( required )

=item lon

C<>=> longitude of the point ( required )

=item radius

C<>=> earth radius in km ( defaults to 6371 )

=back

=cut

=head1 METHODS

=head2 distance_to

 $gc->distance_to( $point[, $precision] )
 $gc->distance_to( { lat => 40.422371, lon => -3.704298 } )
 $gc->distance_to( Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298 ) )

This uses the "haversine" formula to calculate great-circle distances between
the two points - that is, the shortest distance over the earth's surface -
giving an `as-the-crow-flies` distance between the points (ignoring any hills!)

The haversine formula `remains particularly well-conditioned for numerical
computation even at small distances` - unlike calculations based on the spherical
law of cosines. It was published by R W Sinnott in Sky and Telescope, 1984,
though known about for much longer by navigators. (For the curious, c is the
angular distance in radians, and a is the square of half the chord length between
the points).

Returns with the distance using the precision defined or -6
( -6 = 6 decimals ( eg 4.000001 ) )

=cut

=head2 bearing_to

 $gc->bearing_to( $point[, $precision] );
 $gc->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 $gc->bearing_to( Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298 ), -6 );

In general, your current heading will vary as you follow a great circle path
(orthodrome); the final heading will differ from the initial heading by varying
degrees according to distance and latitude (if you were to go from say 35N,45E
(Baghdad) to 35N,135E (Osaka), you would start on a heading of 60 and end up on
a heading of 120!).

This formula is for the initial bearing (sometimes referred to as forward
azimuth) which if followed in a straight line along a great-circle arc will take
you from the start point to the end point

Returns the (initial) bearing from this point to the supplied point, in degrees
with the specified pricision

see http://williams.best.vwh.net/avform.htm#Crs

=cut

=head2 final_bearing_to

 my $f_brng = $gc->final_bearing_to( $point[, $precision] );
 my $f_brng = $gc->final_bearing_to( { lat => 40.422371, lon => -3.704298 } );
 my $f_brng = $gc->final_bearing_to( Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298 ) );

Returns final bearing arriving at supplied destination point from this point;
the final bearing will differ from the initial bearing by varying degrees
according to distance and latitude

=cut

=head2 midpoint_to

 $gc->midpoint_to( $point[, $precision] );
 $gc->midpoint_to( { lat => 40.422371, lon => -3.704298 } );
 $gc->midpoint_to( Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298 ) );

Returns the midpoint along a great circle path between the initial point and
the supplied point.

see http://mathforum.org/library/drmath/view/51822.html for derivation

=cut

=head2 destination_point

 $gc->destination_point( $bearing, $distance[, $precision] );
 $gc->destination_point( 90, 1 );

Returns the destination point and the final bearing using Vincenty inverse
formula for ellipsoids.

=cut

=head2 destination_point_hs

 $gc->destination_point_hs( $bearing, $distance[, $precision] );
 $gc->destination_point_hs( 90, 1 );

Returns the destination point from this point having travelled the given
distance on the given initial bearing (bearing may vary before destination is
reached)

see http://williams.best.vwh.net/avform.htm#LL

=cut

=head2 boundry_box

 $gc->boundry_box( $width[, $height[, $precision]] ); # in km
 $gc->boundry_box( 3, 4 ); # will generate a 3x4m box around the point
 $gc->boundry_box( 1 ); # will generate a 2x2m box around the point (radius)

Returns the boundry box min/max having the initial point defined as the center
of the boundry box, given the widht and height

=cut

=head2 rhumb_distance_to

 $gc->rhumb_distance_to( $point[, $precision] );
 $gc->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 } );
 $gc->rhumb_distance_to( Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298 ) );

Returns the distance from this point to the supplied point, in km, travelling
along a rhumb line.

A 'rhumb line' (or loxodrome) is a path of constant bearing, which crosses all
meridians at the same angle.

Sailors used to (and sometimes still) navigate along rhumb lines since it is
easier to follow a constant compass bearing than to be continually adjusting
the bearing, as is needed to follow a great circle. Rhumb lines are straight
lines on a Mercator Projection map (also helpful for navigation).

Rhumb lines are generally longer than great-circle (orthodrome) routes. For
instance, London to New York is 4% longer along a rhumb line than along a
great circle . important for aviation fuel, but not particularly to sailing
vessels. New York to Beijing . close to the most extreme example possible
(though not sailable!) . is 30% longer along a rhumb line.

see http://williams.best.vwh.net/avform.htm#Rhumb

=cut

=head2 rhumb_bearing_to

 $gc->rhumb_bearing_to( $point[, $precision] );
 $gc->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 } );
 $gc->rhumb_bearing_to( Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298 ) );

Returns the bearing from this point to the supplied point along a rhumb line,
in degrees

=cut

=head2 rhumb_destination_point

 $gc->rhumb_destination_point( $brng, $distance[, $precision] );
 $gc->rhumb_destination_point( 30, 1 );

Returns the destination point from this point having travelled the given distance
(in km) on the given bearing along a rhumb line.

=cut

=head2 intersection

 $gc->intersection( $brng1, $point, $brng2[, $precision] );
 $gc->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180 );
 $gc->intersection( 90, Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298 ), 180 );

Returns the point of intersection of two paths defined by point and bearing

see http://williams.best.vwh.net/avform.htm#Intersection

=cut

=head2 distance_at

Returns the distance in meters for 1deg of latitude and longitude at the
specified latitude

 my $m_distance = $self->distance_at([$precision]);
 my $m_distance = $self->distance_at();
 # at lat 2 with precision -6 returns { m_lat => 110575.625009, m_lon => 111252.098718 }

=cut

=head1 BUGS

All complex software has bugs lurking in it, and this module is no
exception.

Please report any bugs through the web interface at L<http://rt.cpan.org>.

=head1 AUTHOR

Sorin Alexandru Pop C<< <asp@cpan.org> >>

=head1 THANKS

Marius Crisan C<< <crisan.marius@gmail.com> >>

=head1 LICENSE

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut

__END__

1;