#!/bin/bash
head -n 5 supercell.CONFIG > Sed_supercell
sed -n -e '/Si/{N;p;}'  supercell.CONFIG >> Sed_supercell
sed -n -e '/O /{N;p;}'  supercell.CONFIG >> Sed_supercell

