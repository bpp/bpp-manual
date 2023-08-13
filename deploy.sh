#!/bin/bash

if ! [[ -d ./site ]]
then
    mkdir site
fi
if ! [[ -d ./site/docs ]]
then
    mkdir ./site/docs
fi

cp mkdocs.yml ./site/ &&
cp bpp-4-manual.md ./site/docs/ &&
cp quickstart.md ./site/docs/index.md &&
cp mathjax.js ./site/docs/ &&
echo "copied files to site" &&
(cd ./site; mkdocs gh-deploy)
