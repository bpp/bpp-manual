---
hide:
  - navigation
  - toc
---

# BPP Quick Start
Here we provide a quick start for users interested in installing and running BPP on an example dataset for frogs. Once you have succeeded in installing and running the program you can read the documentation to learn how to prepare input files and run the program on your own data.

## Linux
```bash
wget -c https://github.com/bpp/bpp/releases/download/v4.8.0/bpp-4.8.0-linux-x86_64.tar.gz -O - | tar -xz
export PATH=$PATH:$PWD/bpp-4.8.0-linux-x86_64/bin
cd bpp-4.8.0-linux-x86_64/examples/frogs/
bpp --cfile A00.bpp.ctl
```

## Mac x86_64 (Intel/AMD chips)
```bash
wget -c https://github.com/bpp/bpp/releases/download/v4.8.0/bpp-4.8.0-macos-x86_64.tar.gz -O - | tar -xz
export PATH=$PATH:$PWD/bpp-4.8.0-macos-x86_64/bin
cd bpp-4.8.0-macos-x86_64/examples/frogs/
bpp --cfile A00.bpp.ctl
```
## Mac aarch64 (Apple Silicon M1/M2/M3)
```bash
wget -c https://github.com/bpp/bpp/releases/download/v4.8.0/bpp-4.8.0-macos-aarch64.tar.gz -O - | tar -xz
export PATH=$PATH:$PWD/bpp-4.8.0-macos-aarch64/bin
cd bpp-4.8.0-macos-aarch64/examples/frogs/
bpp --cfile A00.bpp.ctl
```

## Windows
```bash
wget -c https://github.com/bpp/bpp/releases/download/v4.8.0/bpp-4.8.0-win-x86_64.zip
```
Extract (unzip) the contents of this file. You will now have the binary distribution in a folder called bpp-4.8.0-win-x86_64. The bpp executable is called bpp.exe.
