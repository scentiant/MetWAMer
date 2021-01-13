#!/usr/bin/env bash

# PLEASE be certain that the libxml2 library is available
# on your system BEFORE running this build script!
# On Ubuntu, issue ``sudo apt install libxml2-dev".

# build_it.sh - Script to fully automate the build process
#               for the MetWAMer package.  A much more elegant
#               solution would be to use GNU Autotools, but
#               many configuration options were buried in the
#               original Makefiles to assist with testing
#               differing variants of the methods, and I
#               don't have time to re-tool the package.
#               That being said, the current approach
#               should work perfectly fine on most systems.

# Michael E Sparks (mespar1@gmail.com)
# Last modified: 20 December 2020

# Copyright (c) 2010,2013 Michael E Sparks
# All rights reserved.

# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

# Build the IMMpractical library
cd IMMpractical/src/
make
cd -

# Build the MetWAMer package
cd src
make
cd -

# Build the BSSM4GSQ package (for training examples)
cd ./data/Bmori/splice_site_models/BSSM4GSQ_forBombyx/src
make
cd -

# Look at what's been built
echo -e "\nExecutables are available in ./bin/\n"

# Hint for cleaning the build
#cd IMMpractical/src/ ; make clean ; cd - ; cd src/ ; make cleaner ; cd -
#rm ./data/Bmori/splice_site_models/BSSM4GSQ_forBombyx/bin/*
