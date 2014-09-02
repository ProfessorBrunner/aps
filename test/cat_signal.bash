#!/bin/bash
cat $1/signal*.dat > $1/x_signal.dat
#rm $1/signal*.dat
mv $1/x_signal.dat $1/signal.dat

cat $1/kl_signal*.dat > $1/x_kl_signal.dat
#rm $1/kl_signal*.dat
mv $1/x_kl_signal.dat $1/kl_signal.dat