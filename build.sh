#!/bin/sh


if [ -f unfft ]
then
  rm unfft
fi

goimports -w .
go fmt
go build

if [ -f unfft ]
then
  ./unfft
fi

feh -d -g+1000+0 fft_c.png
# feh -d -g+1000+0 u.png
# feh -d -g+1000+0 s.png
# feh -d -g+1000+0 un_s.png


#convert s.png n.30-70.png -compose Mathematics -define compose:args='0,1,1,-0.5' -composite c.png
#feh -d -g+1000+0 c.png
