#!/usr/bin/perl
#
# Test for Geo::Calc::XS
#

use strict;
use warnings;
use utf8;

use Test::More;
use Test::Deep;
use Time::HiRes qw( gettimeofday tv_interval );

use_ok( 'Geo::Calc::XS' );

my $gc_m  = Geo::Calc::XS->new( lat => 40.417875, lon => -3.710205, units => 'm' );
my $gc_m0 = Geo::Calc::XS->new( lat => 40.422371, lon => -3.704298, units => 'm' );
my $gc_km = Geo::Calc::XS->new( lat => 40.417875, lon => -3.710205, units => 'k-m' );
my $gc_mi = Geo::Calc::XS->new( lat => 34.159545, lon => -118.243103, units => 'mi' );
my $gc_ft = Geo::Calc::XS->new( lat => 34.159545, lon => -118.243103, units => 'ft' );
{
    ## boundry_box
    # 1x1 km
    is_deeply( $gc_m->boundry_box( 1000, 1000, -6 ), { 'lat_max' => 40.422377, 'lon_max' => -3.704309, 'lon_min' => -3.7161, 'lat_min' => 40.413373 }, 'boundry box 1/1 km' );
    is_deeply( $gc_m->boundry_box( 2000, 2000, -6 ), { 'lat_max' => 40.426878, 'lon_max' => -3.698414, 'lon_min' => -3.721995, 'lat_min' => 40.408872 }, 'boundry box 2/2 km' );
    is_deeply( $gc_m->boundry_box( 2000, 2000, -12 ), { 'lat_max' => 40.426878090262, 'lon_max' => -3.698414102952, 'lon_min' => -3.721995897047, 'lat_min' => 40.408871900362 }, 'boundry box 2/2 km with 12 digits' );

    # 2x2 km == 2000, 2000
    is_deeply( $gc_m->boundry_box( 1000, undef, -3 ), { 'lat_max' => 40.427, 'lon_max' => -3.698, 'lon_min' => -3.721, 'lat_min' => 40.409 }, 'boundry box 2/2 km with 1 args ( 1000, undef, -3 )' );
    is_deeply( $gc_m->boundry_box( 1000, undef ), { 'lat_max' => 40.426878, 'lon_max' => -3.698414, 'lon_min' => -3.721995, 'lat_min' => 40.408872 }, 'boundry box 2/2 km with args ( 1000, undef )' );
    is_deeply( $gc_m->boundry_box( 1000 ), { 'lat_max' => 40.426878, 'lon_max' => -3.698414, 'lon_min' => -3.721995, 'lat_min' => 40.408872 }, 'boundry box 2/2 km with args ( 1000 )' );

    # 6x8 km == 6000x8000 meters
    is_deeply( $gc_m->boundry_box( 6000, 8000 ), { 'lat_max' => 40.453887, 'lon_max' => -3.674832, 'lon_min' => -3.745577, 'lat_min' => 40.381863 }, 'boundry box 6/8 km' );

    ## destination_point
    is_deeply( $gc_m->destination_point( 44.3, 1000, -6 ), { 'lat' => 40.424321, 'lon' => -3.701973, 'final_bearing' => 44.305337 }, 'destination point 1' );
    is_deeply( $gc_m->destination_point( 13.443, 1000, -6 ), $gc_km->destination_point( 13.443, 1, -6 ), 'destination point 2' );

    ## distance_to
    is( $gc_m->distance_to( { lat => 40.422371, lon => -3.704298 } ), 707.106482, 'distance_to 40.422371/-3.704298' );
    is( $gc_m->distance_to( $gc_m0 ), 707.106482, 'distance_to 40.422371/-3.704298 ( Geo::Calc::XS Object )' );
    is( $gc_m->distance_to( { lat => 51.500795, lon => -0.142264 }, -3 ), 1269060.915, 'distance_to buckingham palace' );

    ## midpoint
    is_deeply( $gc_m->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 ) , { 'lat' => 40.420123, 'lon' => -3.707251 }, 'midpoint' );
    is_deeply( $gc_m->midpoint_to( $gc_m0, -6 ) , { 'lat' => 40.420123, 'lon' => -3.707251 }, 'midpoint ( Geo::Calc::XS Object )' );
    is_deeply( $gc_m->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 ) , $gc_km->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 ), 'midpoints are the same' );
    is_deeply( $gc_m->midpoint_to( { lat => 48.149367, lon => 11.748848 }, -6 ) , { 'lat' => 44.543903, 'lon' => 3.506819 }, 'midpoint' );

    ## intersection
    is_deeply( $gc_m->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180, -6 ), { 'lat' => 40.417875, 'lon' => -3.704298 }, 'intersection' );
    is_deeply( $gc_m->intersection( 90, $gc_m0, 180, -6 ), { 'lat' => 40.417875, 'lon' => -3.704298 }, 'intersection ( Geo::Calc::XS Object )' );
    is_deeply( $gc_m->intersection( 43, { lat => 40.729828, lon => -73.883743 }, 12, -6 ), { 'lat' => 54.967178, 'lon' => 85.065586 }, 'over intersection' );

    ## distance_at
    is_deeply( $gc_m->distance_at(), { m_lon => 84871.014948, m_lat => 111042.645811 }, 'distance at latitude' );

    ## initial bearing
    is( $gc_m->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 ), 45.004851, 'initial bearing 1' );
    is( $gc_m->bearing_to( $gc_m0, -6 ), 45.004851, 'initial bearing 1 ( Geo::Calc::XS Object )' );
    is( $gc_m->bearing_to( { lat => 12, lon => -85 }, -6 ), 273.683864, 'initial bearing 2' );
    is( $gc_m->bearing_to( { lat => 1, lon => 10 }, -6 ), 158.973869, 'initial bearing 3' );
    is( $gc_m->bearing_to( { lat => 46, lon => 5 }, -6 ), 45.753222, 'initial bearing 4' );

    ## final bearing
    is( $gc_m->final_bearing_to( { lat => 40.422371, lon => -3.704298 } ), 45.008681, 'final bearing' );
    is( $gc_m->final_bearing_to( $gc_m0 ), 45.008681, 'final bearing ( Geo::Calc::XS Object )' );
    is( $gc_m->final_bearing_to( { lat => 12, lon => -85 } ), 230.962738, 'final bearing' );
    is( $gc_m->final_bearing_to( { lat => 1, lon => 10 } ), 164.144976, 'final bearing' );
    is( $gc_m->final_bearing_to( { lat => 46, lon => 5 } ), 51.729940, 'final bearing' );

    ## using rhumb
    is( $gc_m->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 }, -6 ), 707.094665, 'rhumb distance' );
    is( $gc_m->rhumb_distance_to( $gc_m0, -6 ), 707.094665, 'rhumb distance ( Geo::Calc::XS Object )' );

    is( $gc_m->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 } ), 45.006766, 'rhumb bearing 1' );
    is( $gc_m->rhumb_bearing_to( $gc_m0 ), 45.006766, 'rhumb bearing 1 ( Geo::Calc::XS Object )' );
    is( $gc_m->rhumb_bearing_to( { lat => 12, lon => -85 } ), 248.409098, 'rhumb bearing 2' );
    is( $gc_m->rhumb_bearing_to( { lat => 10, lon => -6 } ), 183.829572, 'rhumb bearing 3' );
    is( $gc_m->rhumb_bearing_to( { lat => 46, lon => 5 } ), 48.644443, 'rhumb bearing 4' );

    is_deeply( $gc_m->rhumb_destination_point( 30, 1000, -6 ), { 'lat' => 40.425663, 'lon' => -3.704298 }, 'rhumb destination point' );
    is_deeply( $gc_m->rhumb_destination_point( 30, 1000, -6 ), $gc_km->rhumb_destination_point( 30, 1, -6 ), 'rhumb destination point' );
}

{
    is_deeply( $gc_m->boundry_box( 1000, 1000, -6 ), $gc_km->boundry_box( 1, 1, -6 ), 'boundry box 1/1 km' );
    is_deeply( $gc_m->boundry_box( 5500, 1200, -6 ), $gc_km->boundry_box( 5.5, 1.2, -6 ), 'boundry box 5.5/1.2 km' );
}

{
    is_deeply( $gc_mi->boundry_box( 1, 1, -12 ), $gc_ft->boundry_box( 5280, 5280, -12 ), 'boundry box 1/1 mi' );
    is( $gc_mi->distance_to( { lat => 40.422371, lon => -3.704298 } ), 6119.644619, 'distance_to 40.422371/-3.704298 in miles' );
    is( $gc_mi->distance_to( $gc_m0 ), 6119.644619, 'distance_to 40.422371/-3.704298 in miles ( Geo::Calc::XS Object )' );
    is( $gc_mi->distance_to( { lat => 40.853293, lon => -73.987427 }, -5 ), 2554.82647, 'distance_to NY -> LA' );
}

done_testing();
