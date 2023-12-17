---
hide:
  - navigation
  - toc
---

# BPP Quick Start
Here we provide a quick start for users interested in installing and running BPP on an example dataset for frogs. Once you have succeeded in installing and running the program you can read the documentation to learn how to prepare input files and run the program on your own data.

## Linux
```
wget -c https://github.com/bpp/bpp/releases/download/v4.7.0/bpp-4.7.0-linux-x86_64.tar.gz -O - | tar -xz
export PATH=$PATH:$PWD/bpp-4.7.0-linux-x86_64/bin
cd bpp-4.7.0-linux-x86_64/examples/frogs/
bpp --cfile A00.bpp.ctl
```

## Mac
```
wget -c https://github.com/bpp/bpp/releases/download/v4.7.0/bpp-4.7.0-macos-x86_64.tar.gz -O - | tar -xz
export PATH=$PATH:$PWD/bpp-4.7.0-macos-x86_64/bin
cd bpp-4.7.0-macos-x86_64/examples/frogs/
bpp --cfile A00.bpp.ctl
```

