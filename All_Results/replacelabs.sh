#!/bin/bash
echo "Replacing labels of .paramnames files"

find ./ -iname '*.paramnames' -exec sed -i -e 's/\\Omega_m/\\Omega_\{\{\\rm m\},0\}/g' {} \;
find ./ -iname '*.paramnames' -exec sed -i -e 's/\\Omega_b/\\Omega_\{\{\\rm b\},0\}/g' {} \;
find ./ -iname '*.paramnames' -exec sed -i -e 's/\\Omega_{DE}/\\Omega_\{\{\\rm DE\},0\}/g' {} \;
#find ./ -iname '*.paramnames' -exec sed -i -e 's/\\sigma_{8}/\\sigma_\{8,0\}/g' {} \;
find ./ -iname '*.paramnames' -exec sed -i -e 's/\\sigma_{8,0}/\\sigma_{8}/g' {} \;
find ./ -iname '*.paramnames' -exec sed -i -e 's/n_s/n_\{\\rm s\}/g' {} \;

echo "Replacement finished"
