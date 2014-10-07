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

{
    my $t0 = [gettimeofday];
    map{ $gc_m->distance_to( { lat => 40.422371, lon => -3.704298 } ) } ( 1..100000 );
    ok( 1, sprintf( 'Took %s to finish 100000 distance calculations', tv_interval ( $t0 ) ) );

    $t0 = [gettimeofday];
    map{ $gc_m->distance_to( $gc_m0 ) } ( 1..100000 );
    ok( 1, sprintf( 'Took %s to finish 100000 distance calculations to anther Geo::Calc::XS object', tv_interval ( $t0 ) ) );
}

done_testing();
