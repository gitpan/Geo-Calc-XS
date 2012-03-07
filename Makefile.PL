#!/usr/bin/env perl

use 5.4.0;
use strict;
use warnings;
use utf8;

use ExtUtils::MakeMaker qw(:DEFAULT);
use File::Spec;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'Geo::Calc::XS',
    VERSION_FROM => File::Spec->catfile( 'lib', 'Geo', 'Calc', 'XS.pm' ), # finds $VERSION
    PREREQ_PM    => {
            'File::Spec'  => 0,
            'Test::More'  => 0,
            'Test::Deep'  => 0,
        },
    LIBS         => [''], # e.g., '-lm'
    DEFINE       => '', # e.g., '-DBYPASS'
    INC          => '', # e.g., '-I/usr/include/other'
    SIGN         => 1,
);
