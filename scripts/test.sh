#!/bin/bash

for k in `cat file.txt`
do
tail -n 4 $k.log > check502.txt
grep "/soft/apps/Paid/g09/l502.exe" check502.txt > exist502.txt
done
if [ -s exist502.txt ]
then
grep -B5 "/soft/apps/Paid/g09/l502.exe" check502.txt > sub502.txt
grep "==>" sub502.txt > sub502name.txt
else
echo "*****************There's no l502 task.********************"
fi

